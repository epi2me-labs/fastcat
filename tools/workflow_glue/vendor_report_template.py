"""Convert a wf-template report into embedded fastcat C report assets."""

import argparse
import gzip
import json
from pathlib import Path
import re


DOCS_JSON_RE = re.compile(
    r"const docs_json = '(.*?)';\n        const render_items = ",
    flags=re.DOTALL,
)
RENDER_ITEMS_RE = re.compile(
    r"const render_items = (\[.*?\]);\n        root\.Bokeh\.embed\.embed_items",
    flags=re.DOTALL,
)
DOCS_JSON_PLACEHOLDER = "__FASTCAT_BOKEH_DOCS_JSON__"
BOUND_VALUE_SENTINEL = "__FASTCAT_DEADBEEF__"
DOM_UUID_RE = re.compile(r"[0-9a-f]{32}")

MIN_EXPECTED_FIGURES = 17


def deterministic_gzip(data):
    """Compress data without Python-version-specific gzip header metadata."""
    compressed = bytearray(gzip.compress(data, compresslevel=9, mtime=0))
    compressed[9] = 255  # gzip OS field: unknown, rather than the host platform
    return bytes(compressed)


def argparser():
    parser = argparse.ArgumentParser(
        description="patch wf-template HTML and emit embedded C assets",
        add_help=False)
    parser.add_argument("--input", required=True, help="example wf-template report HTML")
    parser.add_argument("--template-output", help="optional patched runtime template HTML")
    parser.add_argument("--header", required=True, help="output header path")
    parser.add_argument("--source", required=True, help="output C source path")
    parser.add_argument(
        "--verify-source",
        help="verify this existing C source contains the generated template")
    parser.add_argument("--symbol", default="fc_report_template_html", help="symbol base name")
    return parser


def bokeh_figures(docs):
    doc = next(iter(docs.values()))
    return [
        root for root in doc.get("roots", [])
        if root.get("type") == "object" and root.get("name") == "Figure"
    ]


def axis_label(fig, side):
    axes = fig.get("attributes", {}).get(side) or []
    for item in axes:
        if item.get("type") == "object" and item.get("name", "").endswith("Axis"):
            return item.get("attributes", {}).get("axis_label")
    return None


def figure_by_axis_labels(figures, x_label, y_label, renderer_count=None, source_count=None):
    matches = []
    for fig in figures:
        if axis_label(fig, "below") != x_label:
            continue
        if axis_label(fig, "left") != y_label:
            continue
        if renderer_count is not None:
            renderers = fig.get("attributes", {}).get("renderers", [])
            if len(renderers) != renderer_count:
                continue
        if source_count is not None and len(embedded_sources(fig)) != source_count:
            continue
        matches.append(fig)
    if len(matches) != 1:
        raise RuntimeError(
            f"expected exactly one figure for x={x_label!r}, y={y_label!r}, "
            f"renderers={renderer_count!r}, sources={source_count!r}; found {len(matches)}")
    return matches[0]


def figure_by_subtitle(figures, subtitle, source_count=None):
    matches = []
    for fig in figures:
        try:
            text = subtitle_model(fig).get("attributes", {}).get("text")
        except RuntimeError:
            continue
        if text != subtitle:
            continue
        if source_count is not None and len(embedded_sources(fig)) != source_count:
            continue
        matches.append(fig)
    if len(matches) != 1:
        raise RuntimeError(
            f"expected exactly one figure for subtitle={subtitle!r}, sources={source_count!r}; "
            f"found {len(matches)}")
    return matches[0]


def embedded_sources(fig):
    sources = []
    for renderer in fig.get("attributes", {}).get("renderers", []):
        source = renderer.get("attributes", {}).get("data_source")
        if (
                isinstance(source, dict)
                and source.get("type") == "object"
                and source.get("name") == "ColumnDataSource"):
            sources.append(source)
    return sources


def first_embedded_source(fig):
    sources = embedded_sources(fig)
    if sources:
        return sources[0]
    raise RuntimeError(f"could not find embedded ColumnDataSource for {fig.get('id')}")


def remove_lineplot_marker_renderers(fig):
    """Drop ezcharts' unconditional lineplot scatter markers from the template."""
    renderers = fig.get("attributes", {}).get("renderers", [])
    fig.get("attributes", {})["renderers"] = [
        renderer for renderer in renderers
        if renderer.get("attributes", {}).get("glyph", {}).get("name") != "Scatter"
    ]


def add_binding(bindings, name, model, attr):
    bindings.append({
        "name": name,
        "model_id": model["id"],
        "attr": attr,
    })


def subtitle_model(fig):
    above = fig.get("attributes", {}).get("above") or []
    for item in above:
        if item.get("type") == "object" and item.get("name") == "Title":
            return item
    raise RuntimeError(f"could not find subtitle Title for {fig.get('id')}")


def named_model(value, name):
    if isinstance(value, dict):
        attrs = value.get("attributes")
        if isinstance(attrs, dict) and attrs.get("name") == name and "id" in value:
            return value
        for item in value.values():
            hit = named_model(item, name)
            if hit is not None:
                return hit
    elif isinstance(value, list):
        for item in value:
            hit = named_model(item, name)
            if hit is not None:
                return hit
    return None


def numeric_range_model(fig, axis):
    range_obj = fig.get("attributes", {}).get(f"{axis}_range")
    if not isinstance(range_obj, dict):
        raise RuntimeError(f"could not find {axis}_range for {fig.get('id')}")
    return range_obj


def factor_range_model(fig):
    range_obj = fig.get("attributes", {}).get("x_range")
    if not isinstance(range_obj, dict) or range_obj.get("name") != "FactorRange":
        raise RuntimeError(f"could not find FactorRange for {fig.get('id')}")
    return range_obj


def heatmap_color_mapper_model(fig):
    renderer = fig.get("attributes", {}).get("renderers", [])[0]
    glyph = renderer.get("attributes", {}).get("glyph", {})
    attrs = glyph.get("attributes", {})
    mapper = attrs.get("color_mapper")
    if not isinstance(mapper, dict):
        fill = attrs.get("fill_color", {})
        mapper = fill.get("transform") if isinstance(fill, dict) else None
    if not isinstance(mapper, dict):
        raise RuntimeError(f"could not find heatmap color mapper for {fig.get('id')}")
    return mapper


def colorbar_title(fig):
    right = fig.get("attributes", {}).get("right") or []
    for item in right:
        if item.get("type") == "object" and item.get("name") == "ColorBar":
            return item.get("attributes", {}).get("title")
    return None


def add_source_binding(bindings, name, fig):
    add_binding(bindings, f"{name}.source.data", first_embedded_source(fig), "data")


def add_embedded_source_binding(bindings, name, fig, index, field):
    sources = embedded_sources(fig)
    if index >= len(sources):
        raise RuntimeError(
            f"could not find embedded ColumnDataSource[{index}] for {fig.get('id')}")
    add_binding(bindings, f"{name}.{field}.data", sources[index], "data")


def add_subtitle_binding(bindings, name, fig, required=True):
    try:
        model = subtitle_model(fig)
    except RuntimeError:
        if not required:
            return
        raise
    add_binding(bindings, f"{name}.subtitle.text", model, "text")


def add_numeric_range_bindings(bindings, name, fig, axis):
    range_obj = numeric_range_model(fig, axis)
    add_binding(bindings, f"{name}.{axis}_range.start", range_obj, "start")
    add_binding(bindings, f"{name}.{axis}_range.end", range_obj, "end")
    add_binding(bindings, f"{name}.{axis}_range.bounds", range_obj, "bounds")


def add_factor_range_binding(bindings, name, fig):
    add_binding(bindings, f"{name}.x_range.factors", factor_range_model(fig), "factors")


def add_heatmap_color_bindings(bindings, name, fig):
    mapper = heatmap_color_mapper_model(fig)
    add_binding(bindings, f"{name}.color_mapper.high", mapper, "high")


def add_figure_attr_binding(bindings, name, fig, attr):
    add_binding(bindings, f"{name}.{attr}", fig, attr)


def add_null_state_bindings(bindings, name, fig, colorbar=False):
    label = named_model(fig, "fastcat_null_label")
    if label is None:
        raise RuntimeError(f"could not find null label for {fig.get('id')}")
    x_axis = named_model(fig, "fastcat_null_xaxis")
    y_axis = named_model(fig, "fastcat_null_yaxis")
    x_grid = named_model(fig, "fastcat_null_xgrid")
    y_grid = named_model(fig, "fastcat_null_ygrid")
    if x_axis is None or y_axis is None or x_grid is None or y_grid is None:
        raise RuntimeError(f"could not find null-state chrome models for {fig.get('id')}")

    add_binding(bindings, f"{name}.null.text", label, "text")
    add_binding(bindings, f"{name}.null.visible", label, "visible")
    add_binding(bindings, f"{name}.x_axis.visible", x_axis, "visible")
    add_binding(bindings, f"{name}.y_axis.visible", y_axis, "visible")
    add_binding(bindings, f"{name}.x_grid.visible", x_grid, "visible")
    add_binding(bindings, f"{name}.y_grid.visible", y_grid, "visible")
    add_binding(bindings, f"{name}.toolbar_location", fig, "toolbar_location")
    if colorbar:
        color_bar = named_model(fig, "fastcat_null_colorbar")
        if color_bar is None:
            raise RuntimeError(f"could not find null-state colorbar for {fig.get('id')}")
        add_binding(bindings, f"{name}.colorbar.visible", color_bar, "visible")


def scrub_bound_values(docs, bindings):
    """Make runtime-filled Bokeh values obvious if C misses a replacement."""
    by_id = {}

    def walk(value):
        if isinstance(value, dict):
            model_id = value.get("id")
            if model_id is not None:
                # Bokeh may include lightweight references and full model objects
                # with the same id; only the full objects carry attributes
                if model_id not in by_id or "attributes" in value:
                    by_id[model_id] = value
            for item in value.values():
                walk(item)
        elif isinstance(value, list):
            for item in value:
                walk(item)

    walk(docs)
    for binding in bindings:
        model = by_id.get(binding["model_id"])
        if model is None:
            raise RuntimeError(f"could not find model {binding['model_id']}")
        attrs = model.get("attributes")
        if not isinstance(attrs, dict):
            raise RuntimeError(f"model {binding['model_id']} has no attributes")
        attrs[binding["attr"]] = BOUND_VALUE_SENTINEL


def canonicalize_bokeh_ids(docs, bindings, html):
    """Replace Bokeh's generated document and model IDs with stable IDs.

    Model IDs are present in both ``docs_json`` and the adjacent
    ``render_items`` JavaScript. Canonicalising only the document would leave
    those structures referring to different models at runtime.
    """
    id_map = {}
    canonical_docs = {}

    for index, (doc_id, document) in enumerate(docs.items()):
        canonical_id = f"d{index}"
        id_map[doc_id] = canonical_id
        canonical_docs[canonical_id] = document

    def collect(value):
        if isinstance(value, dict):
            model_id = value.get("id")
            if isinstance(model_id, str) and model_id not in id_map:
                id_map[model_id] = f"m{len(id_map) - len(canonical_docs)}"
            for item in value.values():
                collect(item)
        elif isinstance(value, list):
            for item in value:
                collect(item)

    def replace(value):
        if isinstance(value, dict):
            model_id = value.get("id")
            if model_id in id_map:
                value["id"] = id_map[model_id]
            for item in value.values():
                replace(item)
        elif isinstance(value, list):
            for item in value:
                replace(item)

    collect(canonical_docs)
    replace(canonical_docs)
    for binding in bindings:
        binding["model_id"] = id_map[binding["model_id"]]

    render_items_match = RENDER_ITEMS_RE.search(html)
    if not render_items_match:
        raise RuntimeError("could not find Bokeh render_items")
    render_items = json.loads(render_items_match.group(1))
    root_count = 0
    for item in render_items:
        for root_id in item.get("roots", {}).values():
            if root_id not in id_map:
                id_map[root_id] = f"r{root_count}"
                root_count += 1

    # The docs_json assignment is replaced below, but render_items remains in
    # the HTML and must use the same canonical document/model/root IDs. Root
    # IDs also occur as HTML div IDs, so replacing quoted values covers both.
    for old, new in sorted(id_map.items(), key=lambda item: len(item[0]), reverse=True):
        html = re.sub(rf"(?<=[\"']){re.escape(old)}(?=[\"'])", new, html)
    return canonical_docs, html


def canonicalize_generated_dom_ids(html):
    """Replace wf-template UUIDs in DOM IDs and all of their references."""
    replacements = {}
    for dom_id in re.findall(r'\bid="([^"]+)"', html):
        for token in DOM_UUID_RE.findall(dom_id):
            if token not in replacements:
                replacements[token] = f"dom{len(replacements)}"
    for old, new in replacements.items():
        html = html.replace(old, new)
    return html


def extract_bokeh_document(html):
    match = DOCS_JSON_RE.search(html)
    if not match:
        raise RuntimeError("could not find Bokeh docs_json")

    docs_json = match.group(1)
    # bokeh embeds docs_json as a single-quoted JS string, so string-valued
    # model attributes carry one extra level of escaping compared with raw JSON
    docs_json = docs_json.replace("\\'", "'")
    docs_json = docs_json.replace("\\\\", "\\")
    docs = json.loads(docs_json)
    figures = bokeh_figures(docs)
    if len(figures) < MIN_EXPECTED_FIGURES:
        raise RuntimeError(
            f"expected at least {MIN_EXPECTED_FIGURES} Bokeh figures, found {len(figures)}")

    bindings = []

    fig_quality = figure_by_axis_labels(figures, "Quality", "Number of reads", renderer_count=1, source_count=1)
    fig_length = figure_by_axis_labels(figures, "Read length / kb", "Number of reads", renderer_count=1, source_count=1)
    fig_yield = figure_by_axis_labels(figures, "Read length / kb", "Yield above length / Gbases", renderer_count=2, source_count=2)
    fig_alignment_accuracy = figure_by_axis_labels(figures, "Alignment Concordance / %", "Number of reads", renderer_count=1, source_count=1)
    fig_alignment_query_coverage = figure_by_axis_labels(figures, "Read alignment coverage / %", "Number of reads", renderer_count=1, source_count=1)
    fig_coverage_depth = figure_by_axis_labels(
        figures,
        "Minimum Depth",
        "Reference Positions Covered / %",
        renderer_count=2)
    fig_time_yield = figure_by_axis_labels(figures, "Elapsed time / hours", "Cumulative base yield / Gbp", renderer_count=2, source_count=2)
    fig_time_yield_bin = figure_by_axis_labels(figures, "Elapsed time / hours", "Base yield / Gbp", renderer_count=2, source_count=2)
    fig_time_reads = figure_by_axis_labels(figures, "Elapsed time / hours", "Cumulative read count", renderer_count=2, source_count=2)
    fig_time_reads_bin = figure_by_axis_labels(figures, "Elapsed time / hours", "Read count", renderer_count=2, source_count=2)
    fig_time_length = figure_by_axis_labels(figures, "Elapsed time / hours", "Read length / kb", renderer_count=8, source_count=2)
    fig_time_qscore = figure_by_axis_labels(figures, "Elapsed time / hours", "Quality", renderer_count=8, source_count=2)
    fig_time_speed = figure_by_axis_labels(figures, "Elapsed time / hours", "Translocation speed / bp/s", renderer_count=8, source_count=2)
    fig_time_accuracy = figure_by_axis_labels(figures, "Elapsed time / hours", "Alignment Concordance / %", renderer_count=8, source_count=2)
    fig_speed_qscore = figure_by_axis_labels(figures, "Translocation speed / bp/s", "Quality", renderer_count=1, source_count=1)
    fig_length_qscore = figure_by_axis_labels(figures, "Read length / kb", "Quality", renderer_count=1, source_count=1)
    fig_accuracy_qscore = named_model(docs, "fastcat_accuracy_qscore")
    if (
            fig_accuracy_qscore is None
            or fig_accuracy_qscore.get("name") != "Figure"
            or len(embedded_sources(fig_accuracy_qscore)) != 1):
        raise RuntimeError("could not identify alignment concordance heatmap")
    fig_channel_reads = None
    fig_channel_mean_length = None
    fig_channel_mean_quality = None
    fig_channel_mean_concordance = None
    for fig in figures:
        if axis_label(fig, "below") != "x" or axis_label(fig, "left") != "y":
            continue
        if len(embedded_sources(fig)) != 2:
            continue
        title = colorbar_title(fig)
        if title == "Read count":
            fig_channel_reads = fig
        elif title == "Mean read length / kb":
            fig_channel_mean_length = fig
        elif title == "Mean quality":
            fig_channel_mean_quality = fig
        elif title == "Mean alignment concordance / Q":
            fig_channel_mean_concordance = fig
    if (
            fig_channel_reads is None
            or fig_channel_mean_length is None
            or fig_channel_mean_quality is None
            or fig_channel_mean_concordance is None):
        raise RuntimeError("could not identify all channel occupancy figures")

    add_source_binding(bindings, "quality", fig_quality)
    add_subtitle_binding(bindings, "quality", fig_quality)
    add_numeric_range_bindings(bindings, "quality", fig_quality, "x")
    add_numeric_range_bindings(bindings, "quality", fig_quality, "y")

    add_source_binding(bindings, "length", fig_length)
    add_subtitle_binding(bindings, "length", fig_length)
    add_numeric_range_bindings(bindings, "length", fig_length, "x")
    add_numeric_range_bindings(bindings, "length", fig_length, "y")

    add_source_binding(bindings, "yield", fig_yield)
    add_embedded_source_binding(bindings, "yield", fig_yield, 1, "hover_source")
    add_subtitle_binding(bindings, "yield", fig_yield)
    add_numeric_range_bindings(bindings, "yield", fig_yield, "x")
    add_numeric_range_bindings(bindings, "yield", fig_yield, "y")

    add_source_binding(bindings, "alignment_accuracy", fig_alignment_accuracy)
    add_subtitle_binding(bindings, "alignment_accuracy", fig_alignment_accuracy)
    add_numeric_range_bindings(bindings, "alignment_accuracy", fig_alignment_accuracy, "x")
    add_numeric_range_bindings(bindings, "alignment_accuracy", fig_alignment_accuracy, "y")

    add_source_binding(bindings, "alignment_query_coverage", fig_alignment_query_coverage)
    add_subtitle_binding(bindings, "alignment_query_coverage", fig_alignment_query_coverage)
    add_numeric_range_bindings(bindings, "alignment_query_coverage", fig_alignment_query_coverage, "x")
    add_numeric_range_bindings(bindings, "alignment_query_coverage", fig_alignment_query_coverage, "y")

    coverage_source = named_model(docs, "fastcat_coverage_source")
    coverage_hover_source = named_model(docs, "fastcat_coverage_hover_source")
    coverage_subtitle = named_model(docs, "fastcat_coverage_subtitle")
    coverage_targets = named_model(docs, "fastcat_coverage_targets")
    if (
            coverage_source is None or coverage_hover_source is None or coverage_subtitle is None
            or coverage_targets is None):
        raise RuntimeError("could not identify coverage target models")
    add_binding(bindings, "coverage_depth.source.data", coverage_source, "data")
    add_binding(bindings, "coverage_depth.hover_source.data", coverage_hover_source, "data")
    add_binding(bindings, "coverage_depth.subtitle.text", coverage_subtitle, "text")
    add_numeric_range_bindings(bindings, "coverage_depth", fig_coverage_depth, "x")
    add_numeric_range_bindings(bindings, "coverage_depth", fig_coverage_depth, "y")
    add_binding(bindings, "coverage_depth.targets.data", coverage_targets, "data")

    add_source_binding(bindings, "length_qscore", fig_length_qscore)
    add_subtitle_binding(bindings, "length_qscore", fig_length_qscore, required=False)
    add_heatmap_color_bindings(bindings, "length_qscore", fig_length_qscore)
    add_numeric_range_bindings(bindings, "length_qscore", fig_length_qscore, "x")
    add_numeric_range_bindings(bindings, "length_qscore", fig_length_qscore, "y")
    add_null_state_bindings(bindings, "length_qscore", fig_length_qscore, colorbar=True)

    add_source_binding(bindings, "speed_qscore", fig_speed_qscore)
    add_subtitle_binding(bindings, "speed_qscore", fig_speed_qscore, required=False)
    add_heatmap_color_bindings(bindings, "speed_qscore", fig_speed_qscore)
    add_numeric_range_bindings(bindings, "speed_qscore", fig_speed_qscore, "x")
    add_numeric_range_bindings(bindings, "speed_qscore", fig_speed_qscore, "y")
    add_null_state_bindings(bindings, "speed_qscore", fig_speed_qscore, colorbar=True)

    add_source_binding(bindings, "accuracy_qscore", fig_accuracy_qscore)
    add_heatmap_color_bindings(bindings, "accuracy_qscore", fig_accuracy_qscore)
    add_figure_attr_binding(bindings, "accuracy_qscore", fig_accuracy_qscore, "width")
    add_numeric_range_bindings(bindings, "accuracy_qscore", fig_accuracy_qscore, "x")
    add_numeric_range_bindings(bindings, "accuracy_qscore", fig_accuracy_qscore, "y")
    add_null_state_bindings(bindings, "accuracy_qscore", fig_accuracy_qscore, colorbar=True)

    remove_lineplot_marker_renderers(fig_time_yield)
    remove_lineplot_marker_renderers(fig_time_yield_bin)
    remove_lineplot_marker_renderers(fig_time_reads)
    remove_lineplot_marker_renderers(fig_time_reads_bin)

    add_source_binding(bindings, "time_yield", fig_time_yield)
    add_embedded_source_binding(bindings, "time_yield", fig_time_yield, 1, "hover_source")
    add_numeric_range_bindings(bindings, "time_yield", fig_time_yield, "x")
    add_numeric_range_bindings(bindings, "time_yield", fig_time_yield, "y")
    add_null_state_bindings(bindings, "time_yield", fig_time_yield)

    add_source_binding(bindings, "time_yield_bin", fig_time_yield_bin)
    add_embedded_source_binding(bindings, "time_yield_bin", fig_time_yield_bin, 1, "hover_source")
    add_numeric_range_bindings(bindings, "time_yield_bin", fig_time_yield_bin, "x")
    add_numeric_range_bindings(bindings, "time_yield_bin", fig_time_yield_bin, "y")
    add_null_state_bindings(bindings, "time_yield_bin", fig_time_yield_bin)

    add_source_binding(bindings, "time_reads", fig_time_reads)
    add_embedded_source_binding(bindings, "time_reads", fig_time_reads, 1, "hover_source")
    add_numeric_range_bindings(bindings, "time_reads", fig_time_reads, "x")
    add_numeric_range_bindings(bindings, "time_reads", fig_time_reads, "y")
    add_null_state_bindings(bindings, "time_reads", fig_time_reads)

    add_source_binding(bindings, "time_reads_bin", fig_time_reads_bin)
    add_embedded_source_binding(bindings, "time_reads_bin", fig_time_reads_bin, 1, "hover_source")
    add_numeric_range_bindings(bindings, "time_reads_bin", fig_time_reads_bin, "x")
    add_numeric_range_bindings(bindings, "time_reads_bin", fig_time_reads_bin, "y")
    add_null_state_bindings(bindings, "time_reads_bin", fig_time_reads_bin)

    add_source_binding(bindings, "time_length", fig_time_length)
    add_embedded_source_binding(bindings, "time_length", fig_time_length, 1, "hover_source")
    add_numeric_range_bindings(bindings, "time_length", fig_time_length, "x")
    add_numeric_range_bindings(bindings, "time_length", fig_time_length, "y")
    add_null_state_bindings(bindings, "time_length", fig_time_length)

    add_source_binding(bindings, "time_qscore", fig_time_qscore)
    add_embedded_source_binding(bindings, "time_qscore", fig_time_qscore, 1, "hover_source")
    add_numeric_range_bindings(bindings, "time_qscore", fig_time_qscore, "x")
    add_numeric_range_bindings(bindings, "time_qscore", fig_time_qscore, "y")
    add_null_state_bindings(bindings, "time_qscore", fig_time_qscore)

    add_source_binding(bindings, "time_speed", fig_time_speed)
    add_embedded_source_binding(bindings, "time_speed", fig_time_speed, 1, "hover_source")
    add_numeric_range_bindings(bindings, "time_speed", fig_time_speed, "x")
    add_numeric_range_bindings(bindings, "time_speed", fig_time_speed, "y")
    add_null_state_bindings(bindings, "time_speed", fig_time_speed)

    add_source_binding(bindings, "time_accuracy", fig_time_accuracy)
    add_embedded_source_binding(bindings, "time_accuracy", fig_time_accuracy, 1, "hover_source")
    add_numeric_range_bindings(bindings, "time_accuracy", fig_time_accuracy, "x")
    add_numeric_range_bindings(bindings, "time_accuracy", fig_time_accuracy, "y")
    add_null_state_bindings(bindings, "time_accuracy", fig_time_accuracy)

    add_source_binding(bindings, "channel_reads", fig_channel_reads)
    add_embedded_source_binding(bindings, "channel_reads", fig_channel_reads, 1, "hover_source")
    add_subtitle_binding(bindings, "channel_reads", fig_channel_reads)
    add_heatmap_color_bindings(bindings, "channel_reads", fig_channel_reads)
    add_figure_attr_binding(bindings, "channel_reads", fig_channel_reads, "width")
    add_numeric_range_bindings(bindings, "channel_reads", fig_channel_reads, "x")
    add_numeric_range_bindings(bindings, "channel_reads", fig_channel_reads, "y")
    add_null_state_bindings(bindings, "channel_reads", fig_channel_reads, colorbar=True)

    add_source_binding(bindings, "channel_mean_length", fig_channel_mean_length)
    add_embedded_source_binding(
        bindings, "channel_mean_length", fig_channel_mean_length, 1, "hover_source")
    add_subtitle_binding(bindings, "channel_mean_length", fig_channel_mean_length)
    add_heatmap_color_bindings(bindings, "channel_mean_length", fig_channel_mean_length)
    add_figure_attr_binding(bindings, "channel_mean_length", fig_channel_mean_length, "width")
    add_numeric_range_bindings(bindings, "channel_mean_length", fig_channel_mean_length, "x")
    add_numeric_range_bindings(bindings, "channel_mean_length", fig_channel_mean_length, "y")
    add_null_state_bindings(bindings, "channel_mean_length", fig_channel_mean_length, colorbar=True)

    add_source_binding(bindings, "channel_mean_quality", fig_channel_mean_quality)
    add_embedded_source_binding(
        bindings, "channel_mean_quality", fig_channel_mean_quality, 1, "hover_source")
    add_subtitle_binding(bindings, "channel_mean_quality", fig_channel_mean_quality)
    add_heatmap_color_bindings(bindings, "channel_mean_quality", fig_channel_mean_quality)
    add_figure_attr_binding(bindings, "channel_mean_quality", fig_channel_mean_quality, "width")
    add_numeric_range_bindings(bindings, "channel_mean_quality", fig_channel_mean_quality, "x")
    add_numeric_range_bindings(bindings, "channel_mean_quality", fig_channel_mean_quality, "y")
    add_null_state_bindings(bindings, "channel_mean_quality", fig_channel_mean_quality, colorbar=True)

    add_source_binding(bindings, "channel_mean_concordance", fig_channel_mean_concordance)
    add_embedded_source_binding(
        bindings, "channel_mean_concordance", fig_channel_mean_concordance, 1, "hover_source")
    add_subtitle_binding(bindings, "channel_mean_concordance", fig_channel_mean_concordance)
    add_heatmap_color_bindings(bindings, "channel_mean_concordance", fig_channel_mean_concordance)
    add_figure_attr_binding(bindings, "channel_mean_concordance", fig_channel_mean_concordance, "width")
    add_numeric_range_bindings(bindings, "channel_mean_concordance", fig_channel_mean_concordance, "x")
    add_numeric_range_bindings(bindings, "channel_mean_concordance", fig_channel_mean_concordance, "y")
    add_null_state_bindings(bindings, "channel_mean_concordance", fig_channel_mean_concordance, colorbar=True)

    scrub_bound_values(docs, bindings)
    assignment = f"const docs_json = {DOCS_JSON_PLACEHOLDER};\n        const render_items = "
    html = html[:match.start()] + assignment + html[match.end():]
    docs, html = canonicalize_bokeh_ids(docs, bindings, html)
    docs_json = json.dumps(docs, separators=(",", ":"))
    return html, docs_json, bindings


def remove_tag_block_containing(html, tag, marker):
    marker_pos = html.find(marker)
    if marker_pos < 0:
        return html

    tag_open = f"<{tag}"
    tag_close = f"</{tag}>"
    start = html.rfind(tag_open, 0, marker_pos)
    if start < 0:
        raise RuntimeError(f"could not find <{tag}> containing {marker!r}")

    pos = start
    depth = 0
    while True:
        next_open = html.find(tag_open, pos)
        next_close = html.find(tag_close, pos)
        if next_close < 0:
            raise RuntimeError(f"could not find closing </{tag}> for {marker!r}")
        if next_open >= 0 and next_open < next_close:
            depth += 1
            pos = next_open + len(tag_open)
            continue
        depth -= 1
        pos = next_close + len(tag_close)
        if depth == 0:
            return html[:start] + html[pos:]


def remove_versions_panel(html):
    html = html.replace(
        '<div class=""><a class="dropdown-item" href="#versions">Versions</a></div>',
        "",
    )
    return remove_tag_block_containing(html, "section", 'id="versions"')


def replace_params_table(html):
    pattern = re.compile(
        r'(<section[^>]*id="parameters"[^>]*>.*?<h2[^>]*>)Workflow parameters'
        r'(</h2>.*?<tbody>\s*).*?(\s*</tbody>)',
        flags=re.DOTALL,
    )
    html, count = pattern.subn(
        r"\1Run configuration\2__FASTCAT_RUN_CONFIG_ROWS__\3",
        html,
        count=1,
    )
    if count != 1:
        raise RuntimeError("could not find workflow parameters table")
    return html


def patch_conditional_section(html, heading, section_id, kind, display_placeholder):
    """Give a report section and its menu link stable conditional markers."""
    heading_marker = f">{heading}</h2>"
    heading_pos = html.find(heading_marker)
    if heading_pos < 0:
        raise RuntimeError(f"could not find {heading!r} report section")
    section_start = html.rfind("<section", 0, heading_pos)
    section_end = html.find(">", section_start)
    if section_start < 0 or section_end < 0:
        raise RuntimeError(f"could not find section wrapper for {heading!r}")
    opening = html[section_start:section_end + 1]
    id_match = re.search(r'\bid="([^"]+)"', opening)
    if not id_match:
        raise RuntimeError(f"could not find section id for {heading!r}")
    old_id = id_match.group(1)
    opening = opening.replace(f'id="{old_id}"', f'id="{section_id}"', 1)
    opening = opening[:-1] + (
        f' data-fastcat-section="{kind}" '
        f'style="{display_placeholder}">')
    html = html[:section_start] + opening + html[section_end + 1:]

    nav_pattern = re.compile(
        rf'<a(?P<attrs>[^>]*href="#{re.escape(old_id)}"[^>]*)>'
        rf'{re.escape(heading)}</a>')
    nav_match = nav_pattern.search(html)
    if not nav_match:
        raise RuntimeError(f"could not find navigation link for {heading!r}")
    attrs = nav_match.group("attrs").replace(
        f'href="#{old_id}"', f'href="#{section_id}"')
    replacement = (
        f'<a{attrs} data-fastcat-section="{kind}" '
        f'style="{display_placeholder}">{heading}</a>')
    return html[:nav_match.start()] + replacement + html[nav_match.end():]


def replace_alignment_values(html):
    """Replace scaffold alignment values with runtime HTML placeholders."""
    placeholders = {
        "total-reads": "__FASTCAT_ALIGNMENT_TOTAL_READS_HTML__",
        "mapped-reads": "__FASTCAT_ALIGNMENT_MAPPED_READS_HTML__",
        "unmapped-reads": "__FASTCAT_ALIGNMENT_UNMAPPED_READS_HTML__",
        "mapping-rate": "__FASTCAT_ALIGNMENT_MAPPING_RATE_HTML__",
    }
    for key, placeholder in placeholders.items():
        pattern = re.compile(
            rf'(<td[^>]*data-fastcat-alignment-value="{re.escape(key)}"[^>]*>)'
            r'.*?(</td>)',
            flags=re.DOTALL)
        html, count = pattern.subn(rf"\1{placeholder}\2", html, count=1)
        if count != 1:
            raise RuntimeError(f"could not replace alignment value {key!r}")

    details_pattern = re.compile(
        r'(<tbody[^>]*data-fastcat-alignment-details=""[^>]*>)'
        r'.*?(</tbody>)',
        flags=re.DOTALL)
    html, count = details_pattern.subn(
        r"\1__FASTCAT_ALIGNMENT_DETAIL_ROWS__\2", html, count=1)
    if count != 1:
        raise RuntimeError("could not replace alignment detail rows")
    return html


def replace_coverage_tab_buttons(html):
    """Replace the scaffold's placeholder tab with runtime-generated tabs."""
    pattern = re.compile(
        r'<li[^>]*>\s*<button[^>]*>__FASTCAT_COVERAGE_TAB_BUTTONS__</button>\s*</li>',
        flags=re.DOTALL)
    html, count = pattern.subn("__FASTCAT_COVERAGE_TAB_BUTTONS__", html, count=1)
    if count != 1:
        raise RuntimeError("could not replace coverage tab buttons")
    return html


def patch_accuracy_tabs(html):
    """Add runtime visibility markers to optional alignment-derived tabs."""
    tabs = (
        ("Alignment Concordance", "__FASTCAT_TIME_ACCURACY_DISPLAY__"),
        ("Alignment concordance", "__FASTCAT_ACCURACY_QSCORE_DISPLAY__"),
        ("Mean alignment concordance", "__FASTCAT_CHANNEL_CONCORDANCE_DISPLAY__"),
    )
    for index, (label, placeholder) in enumerate(tabs):
        pattern = re.compile(rf'(<button[^>]*)>{re.escape(label)}(</button>)')
        replacement = (
            rf'\1 data-fastcat-accuracy-tab="{index}" style="{placeholder}">'
            rf'{label}\2')
        html, count = pattern.subn(replacement, html, count=1)
        if count != 1:
            raise RuntimeError(f"could not find {label!r} accuracy tab")
    return html


def apply_fastcat_template_wording(html):
    html = html.replace("fastcat sequencing report", "Fastcat Sequencing QC Report")
    html = html.replace(
        "Results generated through the fastcat workflow provided by Oxford Nanopore Technologies.",
        "Sequencing QC generated by fastcat.",
    )
    html = remove_versions_panel(html)
    html = replace_params_table(html)
    html = replace_alignment_values(html)
    html = replace_coverage_tab_buttons(html)
    html = patch_accuracy_tabs(html)
    html = patch_conditional_section(
        html,
        "Alignment statistics",
        "fastcat-alignment-section",
        "alignment",
        "__FASTCAT_ALIGNMENT_DISPLAY__")
    html = patch_conditional_section(
        html,
        "Depth of coverage",
        "fastcat-coverage-section",
        "coverage",
        "__FASTCAT_COVERAGE_DISPLAY__")
    html = html.replace(" nextflow workflow (", " (")
    # The template's section menu is custom markup, not a Bootstrap dropdown.
    # Bootstrap ScrollSpy treats its items as dropdown entries and crashes when
    # trying to activate their nonexistent .dropdown parent.
    html = html.replace(' data-bs-spy="scroll"', "")
    return html


def bytes_to_c_array(data, width=16):
    parts = [f"0x{b:02x}" for b in data]
    lines = []
    for i in range(0, len(parts), width):
        lines.append("    " + ", ".join(parts[i : i + width]))
    return ",\n".join(lines)


def c_string(text):
    return json.dumps(text)


def emit_bindings(bindings, symbol):
    lines = [
        f"const fc_report_template_binding {symbol}_bindings[] = {{",
    ]
    for binding in bindings:
        lines.append(
            "    { "
            f"{c_string(binding['name'])}, "
            f"{c_string(binding['model_id'])}, "
            f"{c_string(binding['attr'])} "
            "},")
    lines.append("};")
    lines.append("")
    lines.append(f"const size_t {symbol}_binding_count = {len(bindings)};")
    return "\n".join(lines)


def c_array_bytes(source_text, symbol, name):
    match = re.search(
        rf"const unsigned char {re.escape(symbol)}_{re.escape(name)}\[\] = \{{\n(.*?)\n\}};",
        source_text,
        flags=re.DOTALL,
    )
    if not match:
        raise RuntimeError(f"could not find {symbol}_{name} in generated C source")
    return bytes(int(value, 16) for value in re.findall(r"0x[0-9a-f]{2}", match.group(1)))


def verify_c_assets(html, docs_json, bindings, source_path, symbol):
    """Verify generated asset contents without depending on zlib byte output."""
    source_text = source_path.read_text(encoding="utf-8")
    raw = html.encode("utf-8")
    docs_raw = docs_json.encode("utf-8")
    gz = c_array_bytes(source_text, symbol, "gz")
    docs_gz = c_array_bytes(source_text, symbol, "docs_json_gz")

    if gzip.decompress(gz) != raw:
        raise RuntimeError("embedded report template does not match generated HTML")
    if gzip.decompress(docs_gz) != docs_raw:
        raise RuntimeError("embedded Bokeh document does not match generated document")

    expected_bindings = emit_bindings(bindings, symbol)
    bindings_match = re.search(
        rf"(const fc_report_template_binding {re.escape(symbol)}_bindings\[\] = \{{.*?"
        rf"const size_t {re.escape(symbol)}_binding_count = \d+;)",
        source_text,
        flags=re.DOTALL,
    )
    if not bindings_match or bindings_match.group(1) != expected_bindings:
        raise RuntimeError("embedded Bokeh bindings do not match generated bindings")

    expected_sizes = {
        "gz_len": len(gz),
        "len": len(raw),
        "docs_json_gz_len": len(docs_gz),
        "docs_json_len": len(docs_raw),
    }
    for name, value in expected_sizes.items():
        declaration = f"const size_t {symbol}_{name} = {value};"
        if declaration not in source_text:
            raise RuntimeError(f"embedded report template has an incorrect {name}")


def emit_c_assets(html, docs_json, bindings, header_path, source_path, symbol):
    raw = html.encode("utf-8")
    compressed = deterministic_gzip(raw)
    array_payload = bytes_to_c_array(compressed)

    docs_raw = docs_json.encode("utf-8")
    docs_compressed = deterministic_gzip(docs_raw)
    docs_array_payload = bytes_to_c_array(docs_compressed)

    guard = f"{symbol.upper()}_DATA_H"
    header_text = f"""#ifndef {guard}
#define {guard}

#include <stddef.h>

typedef struct {{
    const char *name;
    const char *model_id;
    const char *attr;
}} fc_report_template_binding;

extern const unsigned char {symbol}_gz[];
extern const size_t {symbol}_gz_len;
extern const size_t {symbol}_len;

extern const unsigned char {symbol}_docs_json_gz[];
extern const size_t {symbol}_docs_json_gz_len;
extern const size_t {symbol}_docs_json_len;

extern const fc_report_template_binding {symbol}_bindings[];
extern const size_t {symbol}_binding_count;

#endif
"""

    source_text = f"""#include "{header_path.name}"

/*
 * this template is emitted as gzip-compressed bytes rather than a C string so
 * arbitrary HTML, JavaScript, CSS, quotes, backslashes, and embedded data do
 * not need fragile C string escaping, and the executable does not carry the
 * full uncompressed report template in .rodata
 */

const size_t {symbol}_gz_len = {len(compressed)};
const size_t {symbol}_len = {len(raw)};
const size_t {symbol}_docs_json_gz_len = {len(docs_compressed)};
const size_t {symbol}_docs_json_len = {len(docs_raw)};

{emit_bindings(bindings, symbol)}

const unsigned char {symbol}_gz[] = {{
{array_payload}
}};

const unsigned char {symbol}_docs_json_gz[] = {{
{docs_array_payload}
}};
"""

    header_path.parent.mkdir(parents=True, exist_ok=True)
    source_path.parent.mkdir(parents=True, exist_ok=True)
    header_path.write_text(header_text, encoding="utf-8")
    source_path.write_text(source_text, encoding="utf-8")
    print(f"wrote {header_path}")
    print(f"wrote {source_path}")


def main(args):
    html = Path(args.input).read_text(encoding="utf-8")
    html = apply_fastcat_template_wording(html)
    html = canonicalize_generated_dom_ids(html)
    html, docs_json, bindings = extract_bokeh_document(html)

    if args.template_output:
        template_path = Path(args.template_output)
        template_path.parent.mkdir(parents=True, exist_ok=True)
        template_path.write_text(html, encoding="utf-8")

    emit_c_assets(html, docs_json, bindings, Path(args.header), Path(args.source), args.symbol)
    if args.verify_source:
        verify_c_assets(html, docs_json, bindings, Path(args.verify_source), args.symbol)
    return 0
