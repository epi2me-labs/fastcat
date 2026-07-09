"""Convert a wf-template report into embedded fastcat C report assets."""

import argparse
import gzip
import json
from pathlib import Path
import re


DOCS_JSON_RE = re.compile(
    r"const docs_json = '(.*?)';\n        const render_items",
    flags=re.DOTALL,
)
DOCS_JSON_PLACEHOLDER = "__FASTCAT_BOKEH_DOCS_JSON__"
BOUND_VALUE_SENTINEL = "__FASTCAT_DEADBEEF__"

MIN_EXPECTED_FIGURES = 14


def argparser():
    parser = argparse.ArgumentParser(
        description="patch wf-template HTML and emit embedded C assets",
        add_help=False)
    parser.add_argument("--input", required=True, help="example wf-template report HTML")
    parser.add_argument("--template-output", help="optional patched runtime template HTML")
    parser.add_argument("--header", required=True, help="output header path")
    parser.add_argument("--source", required=True, help="output C source path")
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
    fig_yield = figure_by_axis_labels(figures, "Read length / kb", "Yield above length / Gbases", renderer_count=1, source_count=1)
    fig_time_yield = figure_by_axis_labels(figures, "Elapsed time / hours", "Cumulative base yield / Gbp", renderer_count=1, source_count=1)
    fig_time_yield_bin = figure_by_axis_labels(figures, "Elapsed time / hours", "Base yield / Gbp", renderer_count=1, source_count=1)
    fig_time_reads = figure_by_axis_labels(figures, "Elapsed time / hours", "Cumulative read count", renderer_count=1, source_count=1)
    fig_time_reads_bin = figure_by_axis_labels(figures, "Elapsed time / hours", "Read count", renderer_count=1, source_count=1)
    fig_time_length = figure_by_axis_labels(figures, "Elapsed time / hours", "Read length / kb", renderer_count=8, source_count=1)
    fig_time_qscore = figure_by_axis_labels(figures, "Elapsed time / hours", "Quality", renderer_count=8, source_count=1)
    fig_time_speed = figure_by_axis_labels(figures, "Elapsed time / hours", "Translocation speed / bp/s", renderer_count=8, source_count=1)
    fig_speed_qscore = figure_by_axis_labels(figures, "Translocation speed / bp/s", "Quality", renderer_count=1, source_count=1)
    fig_length_qscore = figure_by_axis_labels(figures, "Read length / kb", "Quality", renderer_count=1, source_count=1)
    fig_channel_reads = None
    fig_channel_bases = None
    for fig in figures:
        if axis_label(fig, "below") != "x" or axis_label(fig, "left") != "y":
            continue
        if len(embedded_sources(fig)) != 2:
            continue
        title = colorbar_title(fig)
        if title == "read count":
            fig_channel_reads = fig
        elif title == "base count / 1 M":
            fig_channel_bases = fig
    if fig_channel_reads is None or fig_channel_bases is None:
        raise RuntimeError("could not identify both channel occupancy figures")

    add_source_binding(bindings, "quality", fig_quality)
    add_subtitle_binding(bindings, "quality", fig_quality)
    add_numeric_range_bindings(bindings, "quality", fig_quality, "x")
    add_numeric_range_bindings(bindings, "quality", fig_quality, "y")

    add_source_binding(bindings, "length", fig_length)
    add_subtitle_binding(bindings, "length", fig_length)
    add_numeric_range_bindings(bindings, "length", fig_length, "x")
    add_numeric_range_bindings(bindings, "length", fig_length, "y")

    add_source_binding(bindings, "yield", fig_yield)
    add_subtitle_binding(bindings, "yield", fig_yield)
    add_numeric_range_bindings(bindings, "yield", fig_yield, "x")
    add_numeric_range_bindings(bindings, "yield", fig_yield, "y")

    add_source_binding(bindings, "length_qscore", fig_length_qscore)
    add_subtitle_binding(bindings, "length_qscore", fig_length_qscore, required=False)
    add_heatmap_color_bindings(bindings, "length_qscore", fig_length_qscore)
    add_numeric_range_bindings(bindings, "length_qscore", fig_length_qscore, "x")
    add_numeric_range_bindings(bindings, "length_qscore", fig_length_qscore, "y")

    add_source_binding(bindings, "speed_qscore", fig_speed_qscore)
    add_subtitle_binding(bindings, "speed_qscore", fig_speed_qscore, required=False)
    add_heatmap_color_bindings(bindings, "speed_qscore", fig_speed_qscore)
    add_numeric_range_bindings(bindings, "speed_qscore", fig_speed_qscore, "x")
    add_numeric_range_bindings(bindings, "speed_qscore", fig_speed_qscore, "y")

    remove_lineplot_marker_renderers(fig_time_yield)
    remove_lineplot_marker_renderers(fig_time_yield_bin)
    remove_lineplot_marker_renderers(fig_time_reads)
    remove_lineplot_marker_renderers(fig_time_reads_bin)

    add_source_binding(bindings, "time_yield", fig_time_yield)
    add_numeric_range_bindings(bindings, "time_yield", fig_time_yield, "x")
    add_numeric_range_bindings(bindings, "time_yield", fig_time_yield, "y")

    add_source_binding(bindings, "time_yield_bin", fig_time_yield_bin)
    add_numeric_range_bindings(bindings, "time_yield_bin", fig_time_yield_bin, "x")
    add_numeric_range_bindings(bindings, "time_yield_bin", fig_time_yield_bin, "y")

    add_source_binding(bindings, "time_reads", fig_time_reads)
    add_numeric_range_bindings(bindings, "time_reads", fig_time_reads, "x")
    add_numeric_range_bindings(bindings, "time_reads", fig_time_reads, "y")

    add_source_binding(bindings, "time_reads_bin", fig_time_reads_bin)
    add_numeric_range_bindings(bindings, "time_reads_bin", fig_time_reads_bin, "x")
    add_numeric_range_bindings(bindings, "time_reads_bin", fig_time_reads_bin, "y")

    add_source_binding(bindings, "time_length", fig_time_length)
    add_numeric_range_bindings(bindings, "time_length", fig_time_length, "x")
    add_numeric_range_bindings(bindings, "time_length", fig_time_length, "y")

    add_source_binding(bindings, "time_qscore", fig_time_qscore)
    add_numeric_range_bindings(bindings, "time_qscore", fig_time_qscore, "x")
    add_numeric_range_bindings(bindings, "time_qscore", fig_time_qscore, "y")

    add_source_binding(bindings, "time_speed", fig_time_speed)
    add_numeric_range_bindings(bindings, "time_speed", fig_time_speed, "x")
    add_numeric_range_bindings(bindings, "time_speed", fig_time_speed, "y")

    add_source_binding(bindings, "channel_reads", fig_channel_reads)
    add_embedded_source_binding(bindings, "channel_reads", fig_channel_reads, 1, "hover_source")
    add_subtitle_binding(bindings, "channel_reads", fig_channel_reads)
    add_heatmap_color_bindings(bindings, "channel_reads", fig_channel_reads)
    add_figure_attr_binding(bindings, "channel_reads", fig_channel_reads, "width")
    add_numeric_range_bindings(bindings, "channel_reads", fig_channel_reads, "x")
    add_numeric_range_bindings(bindings, "channel_reads", fig_channel_reads, "y")

    add_source_binding(bindings, "channel_bases", fig_channel_bases)
    add_embedded_source_binding(bindings, "channel_bases", fig_channel_bases, 1, "hover_source")
    add_subtitle_binding(bindings, "channel_bases", fig_channel_bases)
    add_heatmap_color_bindings(bindings, "channel_bases", fig_channel_bases)
    add_figure_attr_binding(bindings, "channel_bases", fig_channel_bases, "width")
    add_numeric_range_bindings(bindings, "channel_bases", fig_channel_bases, "x")
    add_numeric_range_bindings(bindings, "channel_bases", fig_channel_bases, "y")

    scrub_bound_values(docs, bindings)
    docs_json = json.dumps(docs, separators=(",", ":"))
    assignment = f"const docs_json = {DOCS_JSON_PLACEHOLDER};\n        const render_items"
    html = html[:match.start()] + assignment + html[match.end():]
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


def apply_fastcat_template_wording(html):
    html = html.replace("fastcat sequencing report", "Fastcat Sequencing QC Report")
    html = html.replace(
        "Results generated through the fastcat workflow provided by Oxford Nanopore Technologies.",
        "Sequencing QC generated by fastcat.",
    )
    html = remove_versions_panel(html)
    html = replace_params_table(html)
    html = html.replace(" nextflow workflow (", " (")
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


def emit_c_assets(html, docs_json, bindings, header_path, source_path, symbol):
    raw = html.encode("utf-8")
    compressed = gzip.compress(raw, compresslevel=9, mtime=0)
    array_payload = bytes_to_c_array(compressed)

    docs_raw = docs_json.encode("utf-8")
    docs_compressed = gzip.compress(docs_raw, compresslevel=9, mtime=0)
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
    html, docs_json, bindings = extract_bokeh_document(html)

    if args.template_output:
        template_path = Path(args.template_output)
        template_path.parent.mkdir(parents=True, exist_ok=True)
        template_path.write_text(html, encoding="utf-8")

    emit_c_assets(html, docs_json, bindings, Path(args.header), Path(args.source), args.symbol)
    return 0
