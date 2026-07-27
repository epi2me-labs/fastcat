"""Create workflow report."""
import json
from pathlib import Path

from bokeh.models import (
    ColumnDataSource, CustomJS, HoverTool, LinearColorMapper, Range1d, Title)
from dominate.tags import b, p, script, table, tbody, td, th, thead, tr
from dominate.util import raw
import ezcharts as ezc
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Grid, Tabs
from ezcharts.layout.snippets.table import DataTable
import numpy as np
import pandas as pd

from . import bokeh_relational
from . import bokeh_matrix
from . import bokeh_null
from . import bokeh_ribbon
from . import bokeh_seqsummary
from .bokeh_style import patch_bokeh_plot_class
from .util import get_main_logger, get_named_logger, wf_parser  # noqa: ABS101


CHANNEL_MAP_ROOT = Path(__file__).resolve().parent / "data"
FASTCAT_BLUE = "#0173b2"
CHANNEL_MAP_PATHS = {
    "flongle": CHANNEL_MAP_ROOT / "asicb512fx_flongle_channel_map.csv",
    "minion": CHANNEL_MAP_ROOT / "asicb512fx_channel_map.csv",
    "promethion": CHANNEL_MAP_ROOT / "3000_opx_channel_map.csv",
}


patch_bokeh_plot_class(ezc.plots.BokehPlot)


COVERAGE_TAB_SCRIPT = """
document.addEventListener("click", (event) => {
    const tab = event.target.closest("[data-fastcat-coverage-target]");
    if (tab == null) {
        return;
    }

    const index = Number.parseInt(tab.dataset.fastcatCoverageTarget, 10);
    const doc = Bokeh.documents[0];
    const source = doc?.get_model_by_name("fastcat_coverage_source");
    const hoverSource = doc?.get_model_by_name("fastcat_coverage_hover_source");
    const targets = doc?.get_model_by_name("fastcat_coverage_targets");
    const subtitle = doc?.get_model_by_name("fastcat_coverage_subtitle");
    const xRange = doc?.get_model_by_name("fastcat_coverage_x_range");
    const yRange = doc?.get_model_by_name("fastcat_coverage_y_range");
    if (!Number.isInteger(index)
            || [source, hoverSource, targets, subtitle, xRange, yRange].some((model) => model == null)
            || index < 0 || index >= targets.data.x.length) {
        return;
    }

    source.data = {
        x: Array.from(targets.data.x[index]),
        y: Array.from(targets.data.y[index]),
    };
    hoverSource.data = {
        x: Array.from(targets.data.hover_x[index]),
        y: Array.from(targets.data.hover_y[index]),
        width: Array.from(targets.data.hover_width[index]),
        height: Array.from(targets.data.hover_height[index]),
    };
    subtitle.text = targets.data.subtitle[index];
    xRange.start = 1;
    xRange.end = targets.data.x_end[index];
    xRange.bounds = [0, targets.data.x_end[index]];
    yRange.start = 0;
    yRange.end = targets.data.y_end[index];
    yRange.bounds = [0, targets.data.y_end[index]];
    source.change.emit();
    hoverSource.change.emit();

    document.querySelectorAll("[data-fastcat-coverage-target]").forEach((button) => {
        const active = button === tab;
        button.classList.toggle("active", active);
        button.setAttribute("aria-selected", active ? "true" : "false");
    });
});
"""


def load_dense_heatmap(
        df, x, y, value, subtitle, null_message, match_aspect=False,
        hover_data=None, hover_items=None):
    """Load a dense rectangular TSV into a Bokeh heatmap matrix."""
    if df.empty:
        return bokeh_null.null_plot(null_message, x, y)

    data = df.pivot(
        index=x,
        columns=y,
        values=value)
    data = data.rename_axis(x)
    data = data.rename_axis(y, axis=1)
    finite = np.isfinite(data.to_numpy(dtype=float))
    if not finite.any():
        return bokeh_null.null_plot(null_message, x, y)

    return bokeh_matrix.heatmap(
        data,
        annot=False,
        vmax_quantile=0.995,
        cbar_kws={"title": "count"},
        subtitle=subtitle,
        match_aspect=match_aspect,
        hover_data=hover_data,
        hover_items=hover_items)


def aspect_matched_chart_styles(plot):
    """Size an aspect-matched heatmap without stretching its data frame."""
    x_span = float(plot._fig.x_range.end - plot._fig.x_range.start)
    y_span = float(plot._fig.y_range.end - plot._fig.y_range.start)
    aspect = x_span / max(y_span, 1e-12)
    # Allow for the colorbar and surrounding plot chrome.
    max_width = int(round(500.0 * aspect + 120.0))
    return {
        "width": "100%",
        "max-width": f"{max_width}px",
        "margin-left": "auto",
        "margin-right": "auto",
    }


def simple_kv_table(rows):
    """Render a plain key/value table."""
    with table(cls="table table-sm"):
        with thead():
            with tr():
                th("Property")
                th("Value")
        with tbody():
            for key, value in rows:
                with tr():
                    td(str(key))
                    td("" if value is None else str(value))


def load_alignment_counts(path):
    """Load shared flagstat output into report-level alignment counts."""
    df = pd.read_csv(path, sep="\t")
    count_columns = (
        "total", "primary", "secondary", "supplementary", "unmapped",
        "qcfail", "duplicate", "duplex", "duplex_forming")
    for column in count_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce").fillna(0).astype(int)

    mapped = df[df["ref"] != "*"]
    totals = {column: int(df[column].sum()) for column in count_columns}
    mapped_records = int(mapped["total"].sum())
    mapped_primary = int(mapped["primary"].sum())
    unmapped_reads = totals["unmapped"]
    total_reads = mapped_primary + unmapped_reads
    mapping_rate = 100.0 * mapped_primary / total_reads if total_reads else 0.0
    return {
        **totals,
        "mapped_records": mapped_records,
        "mapped_primary": mapped_primary,
        "total_reads": total_reads,
        "mapping_rate": mapping_rate,
    }


def add_alignment_summary(path):
    """Render alignment headline values and overlapping record categories."""
    counts = load_alignment_counts(path)
    headline = (
        ("Total reads", counts["total_reads"], "total-reads"),
        ("Reads mapped", counts["mapped_primary"], "mapped-reads"),
        ("Unmapped reads", counts["unmapped"], "unmapped-reads"),
        ("Mapping rate", f'{counts["mapping_rate"]:.2f}%', "mapping-rate"),
    )
    with table(cls="table table-sm"):
        with thead():
            with tr():
                for label, _, _ in headline:
                    th(label)
        with tbody():
            with tr():
                for _, value, key in headline:
                    td(str(value), **{"data-fastcat-alignment-value": key})

    p("Record categories may overlap; their counts are not expected to sum to the total.")
    detail_rows = (
        ("Mapped records", counts["mapped_records"]),
        ("Primary mappings", counts["primary"]),
        ("Secondary mappings", counts["secondary"]),
        ("Supplementary mappings", counts["supplementary"]),
        ("Total records", counts["total"]),
        ("QC-fail", counts["qcfail"]),
        ("Duplicate", counts["duplicate"]),
        ("Duplex", counts["duplex"]),
        ("Duplex-forming", counts["duplex_forming"]),
    )
    with table(cls="table table-sm"):
        with thead():
            with tr():
                th("Record category")
                th("Count")
        with tbody(**{"data-fastcat-alignment-details": ""}):
            for label, value in detail_rows:
                with tr():
                    td(label)
                    td(str(value))


def add_alignment_plots(hist_dir):
    """Render mapped-primary alignment concordance and query-span histograms."""
    accuracy = bokeh_seqsummary.load_histogram(hist_dir, "accuracy")
    query_coverage = bokeh_seqsummary.load_histogram(hist_dir, "coverage")
    with Grid(columns=2):
        EZChart(bokeh_seqsummary.alignment_histogram_plot(
            accuracy,
            80.0,
            100.0,
            "Alignment Concordance / %",
            color=FASTCAT_BLUE))
        EZChart(bokeh_seqsummary.alignment_histogram_plot(
            query_coverage,
            85.0,
            100.0,
            "Read alignment coverage / %",
            color=FASTCAT_BLUE))


def load_coverage_targets(summary_paths, distribution_paths):
    """Load aggregate legacy coverage tables for the scaffold CDF selector."""
    if len(summary_paths) != len(distribution_paths):
        raise ValueError(
            "coverage summary and distribution path counts must match")
    targets = []
    for summary_path, distribution_path in zip(summary_paths, distribution_paths):
        summary = pd.read_csv(summary_path, sep="\t")
        distribution = pd.read_csv(
            distribution_path,
            sep="\t",
            names=["chrom", "depth", "fraction"])
        name = Path(distribution_path).name.removesuffix(".dist.txt")
        mean_depth = float(summary.iloc[-1]["mean"]) if not summary.empty else 0.0
        distribution = distribution[distribution["depth"] > 0]
        x = distribution["depth"].astype(float).tolist()
        y = (100.0 * distribution["fraction"].astype(float)).tolist()
        width = np.diff(x + [x[-1] + 1.0]).tolist() if x else []
        targets.append({
            "name": name,
            "x": x,
            "y": y,
            "subtitle": f"Mean depth: {mean_depth:.2f}×",
            "x_end": max(x, default=1.0) or 1.0,
            "y_end": 100.0,
            "hover_x": x,
            "hover_y": y,
            "hover_width": width,
            "hover_height": [200.0] * len(x),
        })
    return targets


def add_coverage_plot(summary_paths, distribution_paths):
    """Render one aggregate CDF figure controlled by a target selector."""
    targets = load_coverage_targets(summary_paths, distribution_paths)
    if not targets:
        return

    first = targets[0]
    source = ColumnDataSource(
        data={"x": first["x"], "y": first["y"]},
        name="fastcat_coverage_source")
    hover_source = ColumnDataSource(
        data={"x": first["hover_x"], "y": first["hover_y"],
              "width": first["hover_width"], "height": first["hover_height"]},
        name="fastcat_coverage_hover_source")
    target_source = ColumnDataSource(
        data={key: [target[key] for target in targets] for key in (
            "name", "x", "y", "subtitle", "x_end", "y_end",
            "hover_x", "hover_y", "hover_width", "hover_height")},
        name="fastcat_coverage_targets")
    # Retain the target data source in the serialised Bokeh document; the
    # report's HTML tabs use it to update the visible source at runtime.
    source.js_on_change("data", CustomJS(args={"targets": target_source}, code=""))
    plot = ezc.plots.BokehPlot(
        x_range=Range1d(
            1.0, first["x_end"], name="fastcat_coverage_x_range"),
        y_range=Range1d(
            0.0, first["y_end"], name="fastcat_coverage_y_range"))
    plot._fig.xaxis.axis_label = "Minimum Depth"
    plot._fig.yaxis.axis_label = "Reference Positions Covered / %"
    renderer = plot._fig.line(
        x="x", y="y", source=source, line_width=2, color=FASTCAT_BLUE)
    hover = plot._fig.select_one(HoverTool)
    if hover is None:
        hover = HoverTool()
        plot._fig.add_tools(hover)
    hover.name = "fastcat_coverage_hover"
    bokeh_relational.add_vertical_hover_bands(
        plot, x="x", y="y", width="width", height="height",
        source=hover_source)
    hover.tooltips = [
        ("Minimum Depth", "@x{0}"),
        ("Reference Positions Covered / %", "@y{0.0}"),
    ]
    subtitle = Title(
        text=first["subtitle"],
        name="fastcat_coverage_subtitle")
    plot._fig.add_layout(subtitle, "above")

    tabs = Tabs()
    with tabs.add_tab("__FASTCAT_COVERAGE_TAB_BUTTONS__"):
        EZChart(plot, height="500px")
    script(raw(COVERAGE_TAB_SCRIPT))


def _summary_display_text(value, formatter=str, missing_text="(not set)"):
    if value is None:
        return missing_text
    text = str(value)
    if text == "" or text.lower() == "nan":
        return missing_text
    if text.startswith("__FASTCAT_"):
        return text
    return formatter(value)


def _format_duration_seconds(value):
    total = int(round(float(value)))
    if total < 0:
        return "(not meaningful)"
    hours = total // 3600
    minutes = (total % 3600) // 60
    seconds = total % 60
    return f"{hours}h {minutes}m {seconds}s"


def _format_yield_bases(value):
    value = float(value)
    if value >= 1e9:
        return f"{value / 1e9:.3g} Gbp"
    if value >= 1e6:
        return f"{value / 1e6:.3g} Mbp"
    if value >= 1e3:
        return f"{value / 1e3:.3g} kbp"
    return f"{int(round(value))} bp"


def _format_n50(value):
    value = float(value)
    if value < 1000.0:
        return f"{int(round(value))}"
    return f"{value / 1000.0:.3g} kbp"


def _format_n_range(row, center, upper, lower):
    """Format an N statistic with its high-to-low N10--N90 range."""
    return (
        f"{_summary_display_text(row.get(center), _format_n50)} "
        f"({_summary_display_text(row.get(upper), _format_n50)} - "
        f"{_summary_display_text(row.get(lower), _format_n50)})")


def _format_l_range(row, center, lower, upper):
    """Format an L statistic with its low-to-high L10--L90 range."""
    format_count = lambda value: f"{int(float(value))}"
    return (
        f"{_summary_display_text(row.get(center), format_count)} "
        f"({_summary_display_text(row.get(lower), format_count)} - "
        f"{_summary_display_text(row.get(upper), format_count)})")


def load_run_summary(path):
    """Load fastcat run_summary.tsv and format the top summary table."""
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if df.empty:
        return []
    row = df.iloc[0].to_dict()
    discovered = lambda key, formatter=str: _summary_display_text(
        row.get(key), formatter, "(not found)")
    rows = [
        ("Sample ID", _summary_display_text(row.get("sample_id"))),
        ("Input type", _summary_display_text(row.get("input_type"))),
        ("Barcode ID", _summary_display_text(row.get("barcode_id"))),
        ("Barcode alias", discovered("barcode_alias")),
        ("Run ID", discovered("run_id")),
        ("Basecaller model", discovered("basecaller_model")),
        ("Modbase models", discovered("modbase_models")),
        ("Flowcell ID", discovered("flow_cell_id")),
        ("Run date", discovered("run_date")),
        ("Run duration", discovered("run_duration_s", _format_duration_seconds)),
        ("Yield", _summary_display_text(row.get("base_count"), _format_yield_bases)),
        ("Read count", _summary_display_text(row.get("read_count"), lambda v: f"{int(float(v))}")),
        ("N50 (N10-N90)", _format_n_range(row, "n50_bp", "n10_bp", "n90_bp")),
        ("L50 (L10-L90)", _format_l_range(row, "l50", "l10", "l90")),
    ]
    return rows


def load_length_qscore_heatmap(path):
    """Load fastcat length/qscore TSV into a Bokeh heatmap matrix."""
    df = pd.read_csv(path, sep="\t")
    df["Read length / kb"] = (df["length_bin_start"] + df["length_bin_end"]) * 0.0005
    df["Quality"] = (df["qscore_bin_start"] + df["qscore_bin_end"]) * 0.5
    hover_data = pd.DataFrame({
        "x": df["Read length / kb"],
        "y": df["Quality"],
        "width": (df["length_bin_end"] - df["length_bin_start"]) * 0.001,
        "height": df["qscore_bin_end"] - df["qscore_bin_start"],
        "metric": df["Read length / kb"],
        "quality": df["Quality"],
        "count": df["count"],
    })
    plot = load_dense_heatmap(
        df,
        "Read length / kb",
        "Quality",
        "count",
        "",
        "No read length vs quality data available",
        hover_data=hover_data,
        hover_items=[
            ("Read length / kb", "@metric{0.000}"),
            ("Quality", "@quality{0.0}"),
            ("count", "@count{0.0}"),
        ])
    plot._fig.y_range.end = qscore_heatmap_upper(df)
    return plot


def load_speed_qscore_heatmap(path):
    """Load fastcat qscore/translocation-speed TSV into a Bokeh heatmap matrix."""
    df = pd.read_csv(path, sep="\t")
    df["Translocation speed / bp/s"] = (
        df["speed_bin_start"] + df["speed_bin_end"]) * 0.5
    df["Quality"] = (df["qscore_bin_start"] + df["qscore_bin_end"]) * 0.5
    hover_data = pd.DataFrame({
        "x": df["Translocation speed / bp/s"],
        "y": df["Quality"],
        "width": df["speed_bin_end"] - df["speed_bin_start"],
        "height": df["qscore_bin_end"] - df["qscore_bin_start"],
        "metric": df["Translocation speed / bp/s"],
        "quality": df["Quality"],
        "count": df["count"],
    })
    plot = load_dense_heatmap(
        df,
        "Translocation speed / bp/s",
        "Quality",
        "count",
        "",
        "No translocation speed vs quality data available",
        hover_data=hover_data,
        hover_items=[
            ("Translocation speed / bp/s", "@metric{0.0}"),
            ("Quality", "@quality{0.0}"),
            ("count", "@count{0.0}"),
        ])
    plot._fig.y_range.end = qscore_heatmap_upper(df)
    return plot


def load_accuracy_qscore_heatmap(path):
    """Load fastcat alignment-concordance/qscore TSV into a Bokeh heatmap matrix."""
    df = pd.read_csv(path, sep="\t")
    df["Alignment Concordance / Q"] = (
        df["concordance_qscore_bin_start"]
        + df["concordance_qscore_bin_end"]) * 0.5
    df["Quality"] = (df["qscore_bin_start"] + df["qscore_bin_end"]) * 0.5
    data = df.pivot(
        index="Alignment Concordance / Q",
        columns="Quality",
        values="count")
    hover_data = pd.DataFrame({
        "x": df["Alignment Concordance / Q"],
        "y": df["Quality"],
        "width": df["concordance_qscore_bin_end"]
            - df["concordance_qscore_bin_start"],
        "height": df["qscore_bin_end"] - df["qscore_bin_start"],
        "metric": df["Alignment Concordance / Q"],
        "quality": df["Quality"],
        "count": df["count"],
    })
    plot = bokeh_matrix.heatmap(
        data,
        annot=False,
        vmax_quantile=0.995,
        cbar_kws={"title": "read count"},
        match_aspect=True,
        hide_toolbar=True,
        hover_data=hover_data,
        hover_items=[
            ("Alignment concordance / Q", "@metric{0.0}"),
            ("Quality", "@quality{0.0}"),
            ("count", "@count{0.0}"),
        ])
    plot._fig.name = "fastcat_accuracy_qscore"
    qscore_end = qscore_heatmap_upper(df)
    plot._fig.x_range.end = qscore_end
    plot._fig.y_range.end = qscore_end
    return plot


def qscore_heatmap_upper(df):
    """Choose a Q-score ceiling only when its high-Q tail is material."""
    counts = pd.to_numeric(df["count"], errors="coerce").clip(lower=0.0).fillna(0.0)
    starts = pd.to_numeric(df["qscore_bin_start"], errors="coerce")
    total = counts.sum()
    if total <= 0.0:
        return 40.0
    if counts[starts >= 50.0].sum() / total >= 0.005:
        return 60.0
    if counts[starts >= 40.0].sum() / total >= 0.005:
        return 50.0
    return 40.0


def accuracy_qscore_chart_styles(plot):
    """Apply the channel heatmap's fixed colour-bar allowance to a square matrix."""
    frame_size = 500
    colourbar_width = 120
    plot._fig.width = frame_size + colourbar_width
    plot._fig.height = frame_size
    return {
        "width": "100%",
        "max-width": f"{frame_size + 2 * colourbar_width}px",
        "margin-left": "auto",
        "margin-right": "auto",
    }


def format_hours(value):
    """Format a time in hours for compact categorical bin labels."""
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.2f}".rstrip("0").rstrip(".")


def load_time_bin_summary(path):
    """Load fastcat time_bin_summary.tsv and add plotting columns."""
    df = pd.read_csv(path, sep="\t")
    numeric_cols = (
        "time_bin_start_s", "time_bin_end_s", "read_count", "base_count",
        "cumulative_read_count", "cumulative_base_count",
        "length_mean", "length_q05", "length_q25", "length_q50",
        "length_q75", "length_q95",
        "qscore_mean", "qscore_q05", "qscore_q25", "qscore_q50",
        "qscore_q75", "qscore_q95",
        "translocation_speed_mean",
        "translocation_speed_q05", "translocation_speed_q25", "translocation_speed_q50",
        "translocation_speed_q75", "translocation_speed_q95",
        "accuracy_read_count", "accuracy_mean",
        "accuracy_q05", "accuracy_q25", "accuracy_q50",
        "accuracy_q75", "accuracy_q95")
    for col in numeric_cols:
        if col in df:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.sort_values("time_bin_start_s").reset_index(drop=True)
    df["time_mid_h"] = (
        df["time_bin_start_s"] + df["time_bin_end_s"]) / 7200.0
    df["time_end_h"] = df["time_bin_end_s"] / 3600.0
    labels = [
        f"{format_hours(row.time_bin_start_s / 3600.0)}-"
        f"{format_hours(row.time_bin_end_s / 3600.0)}"
        for row in df.itertuples()
    ]
    df["time_bin"] = pd.Categorical(labels, categories=labels, ordered=True)
    df["base_yield_gbp"] = df["base_count"] / 1e9
    df["cumulative_base_yield_gbp"] = df["cumulative_base_count"] / 1e9
    df["length_mean_kb"] = df["length_mean"] / 1000.0
    df["speed_mean"] = df["translocation_speed_mean"]
    for col in ("q05", "q25", "q50", "q75", "q95"):
        df[f"length_plot_{col}"] = df[f"length_{col}"] / 1000.0
        df[f"qscore_plot_{col}"] = df[f"qscore_{col}"]
        df[f"speed_plot_{col}"] = df[f"translocation_speed_{col}"]
        df[f"accuracy_plot_{col}"] = df[f"accuracy_{col}"]
    return df


def channel_map_name(max_channel):
    if max_channel <= 128:
        return "flongle"
    if max_channel <= 512:
        return "minion"
    return "promethion"


def orient_channel_cells(device, cells):
    if device == "minion":
        cells["x"] = cells["y_raw"]
        cells["y"] = -cells["x_raw"]
    else:
        cells["x"] = cells["x_raw"]
        cells["y"] = -cells["y_raw"]
    return cells


def load_channel_map(name):
    df = pd.read_csv(CHANNEL_MAP_PATHS[name], skipinitialspace=True)
    cells = []
    for row in df.to_dict("records"):
        channel = int(row["name"])
        for mux in range(1, 5):
            x_col = f"mux_{mux}_phys_x"
            y_col = f"mux_{mux}_phys_y"
            if x_col not in df.columns or y_col not in df.columns:
                x_col = "mux_1_phys_x"
                y_col = "mux_1_phys_y"
            cells.append({
                "channel": channel,
                "mux": mux,
                "x_raw": float(row[x_col]),
                "y_raw": float(row[y_col]),
            })
    return orient_channel_cells(name, pd.DataFrame(cells))


def load_channel_summary_heatmap(path, metric):
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return bokeh_null.null_plot("No channel occupancy data available")
    else:
        df["channel"] = pd.to_numeric(df["channel"], errors="coerce").fillna(0).astype(int)
        df["mux"] = pd.to_numeric(df["mux"], errors="coerce").fillna(0).astype(int)
        df["read_count"] = pd.to_numeric(df["read_count"], errors="coerce").fillna(0.0)
        df["base_count"] = pd.to_numeric(df["base_count"], errors="coerce").fillna(0.0)
        if "mean_quality" in df:
            df["mean_quality"] = pd.to_numeric(df["mean_quality"], errors="coerce")
        if "mean_concordance_qscore" in df:
            df["mean_concordance_qscore"] = pd.to_numeric(
                df["mean_concordance_qscore"], errors="coerce")
        max_channel = int(df["channel"].max()) if not df.empty else 1
        unknown = df[df["mux"] == 0].copy()
        known = df[df["mux"] != 0].copy()
        if not unknown.empty:
            spread = []
            for mux in range(1, 5):
                chunk = unknown.copy()
                chunk["mux"] = mux
                chunk["read_count"] = chunk["read_count"] / 4.0
                chunk["base_count"] = chunk["base_count"] / 4.0
                spread.append(chunk)
            df = pd.concat([known] + spread, ignore_index=True)
        else:
            df = known

    cells = load_channel_map(channel_map_name(max_channel))
    merged = cells.merge(df, on=["channel", "mux"], how="left")
    merged["read_count"] = merged["read_count"].fillna(0.0)
    merged["base_count"] = merged["base_count"].fillna(0.0)
    color_title = "read count"
    value_format = "@count{0}"
    if metric == "mean_read_length":
        merged["count"] = np.divide(
            merged["base_count"],
            merged["read_count"],
            out=np.zeros(len(merged), dtype=float),
            where=merged["read_count"] > 0.0) / 1000.0
        color_title = "mean read length / kb"
        value_format = "@count{0.0}"
    elif metric == "mean_quality":
        merged["count"] = merged["mean_quality"].fillna(0.0)
        color_title = "mean quality"
        value_format = "@count{0.0}"
    elif metric == "mean_concordance_qscore":
        merged["count"] = merged["mean_concordance_qscore"].fillna(0.0)
        color_title = "mean alignment concordance / Q"
        value_format = "@count{0.0}"
    else:
        merged["count"] = merged[metric]
    merged["width"] = 1.0
    merged["height"] = 1.0
    data = merged.pivot(index="x", columns="y", values="count")
    data = data.reindex(
        index=np.arange(data.index.min(), data.index.max() + 1.0, 1.0),
        columns=np.arange(data.columns.min(), data.columns.max() + 1.0, 1.0))
    data = data.rename_axis(None)
    data = data.rename_axis(None, axis=1)

    occupied = ((merged["read_count"] > 0.0) | (merged["base_count"] > 0.0)).sum()
    if occupied == 0:
        return bokeh_null.null_plot("No channel occupancy data available")
    subtitle = f"{int(occupied)} sequencing wells"

    return bokeh_matrix.heatmap(
        data,
        annot=False,
        vmax_quantile=0.995,
        cbar_kws={"title": color_title},
        subtitle=subtitle,
        hide_axes=True,
        match_aspect=True,
        hide_toolbar=True,
        hover_data=merged[["x", "y", "width", "height", "channel", "mux", "count"]],
        hover_items=[
            ("channel", "@channel{0}"),
            ("mux", "@mux{0}"),
            (color_title, value_format),
        ])


def time_lineplot(df, y, ylabel, color=FASTCAT_BLUE, x="time_mid_h", seed_origin=False):
    """Create a line plot against elapsed run time in hours."""
    null_messages = {
        "Cumulative base yield / Gbp": "No cumulative base yield through time data available",
        "Base yield / Gbp": "No base yield through time data available",
        "Cumulative read count": "No cumulative read count through time data available",
    }
    plot_df = df[[x, y]].copy()
    plot_df[x] = pd.to_numeric(plot_df[x], errors="coerce")
    plot_df[y] = pd.to_numeric(plot_df[y], errors="coerce")
    plot_df = plot_df.replace([np.inf, -np.inf], np.nan).dropna(subset=[x, y])
    if plot_df.empty:
        return bokeh_null.null_plot(
            null_messages.get(ylabel, f"No {ylabel.lower()} through time data available"),
            "Elapsed time / hours",
            ylabel)
    if seed_origin and not plot_df.empty and plot_df.iloc[0][x] > 0.0:
        first = plot_df.iloc[[0]].copy()
        first[x] = 0.0
        first[y] = 0.0
        plot_df = pd.concat([first, plot_df], ignore_index=True)

    plot = bokeh_relational.lineplot(
        data=plot_df,
        x=x,
        y=y,
        estimator=None,
        errorbar=None,
        color=color,
        hover_items=[
            ("Elapsed time / hours", f"@{{{x}}}{{0.00}}"),
            (ylabel, "@y{0}")])
    plot._fig.xaxis.axis_label = "Elapsed time / hours"
    plot._fig.yaxis.axis_label = ylabel
    x_values = plot_df[x].to_numpy(dtype=float)
    if len(x_values) > 1:
        widths = np.diff(x_values, append=x_values[-1])
        widths[-1] = widths[-2]
    else:
        widths = np.array([1.0])
    hover_df = plot_df.copy()
    hover_df["x"] = x_values
    hover_df["y"] = plot_df[y].to_numpy(dtype=float)
    hover_df["width"] = widths
    hover_df["height"] = max(float(plot_df[y].max()), 1.0) * 2.0
    bokeh_relational.add_vertical_hover_bands(
        plot, hover_df, x="x", y="y", width="width", height="height")
    if ylabel == "Cumulative read count":
        bokeh_seqsummary.apply_scientific_count_formatter(plot._fig.yaxis[0])
    bokeh_null.attach_null_state(plot)
    return plot


def time_stepplot(df, y, ylabel, color=FASTCAT_BLUE):
    """Create a stepped per-bin plot over elapsed run time in hours."""
    null_messages = {
        "Base yield / Gbp": "No base yield through time data available",
        "Read count": "No read count through time data available",
    }
    if df.empty:
        return bokeh_null.null_plot(
            null_messages.get(ylabel, f"No {ylabel.lower()} through time data available"),
            "Elapsed time / hours",
            ylabel)
    plot_df = pd.DataFrame({
        "elapsed time (hours)": (df["time_bin_start_s"] / 3600.0).astype(float),
        ylabel: pd.to_numeric(df[y], errors="coerce").fillna(0.0).astype(float),
    })
    if not df.empty:
        plot_df.loc[len(plot_df)] = {
            "elapsed time (hours)": float(df["time_bin_end_s"].iloc[-1]) / 3600.0,
            ylabel: float(plot_df[ylabel].iloc[-1]),
        }

    plot = bokeh_relational.lineplot(
        data=plot_df,
        x="elapsed time (hours)",
        y=ylabel,
        estimator=None,
        errorbar=None,
        color=color,
        drawstyle="steps-post",
        hover_items=[
            ("Elapsed time / hours", "@{elapsed time (hours)}{0.00}"),
            (ylabel, f"@{{{ylabel}}}{{0.00}}"
             if ylabel == "Base yield / Gbp" else f"@{{{ylabel}}}{{0}}")])
    plot._fig.xaxis.axis_label = "Elapsed time / hours"
    plot._fig.yaxis.axis_label = ylabel
    hover_df = pd.DataFrame({
        "x": (df["time_bin_start_s"] + df["time_bin_end_s"]) / 7200.0,
        "y": plot_df[ylabel].iloc[:len(df)].to_numpy(),
        "width": (df["time_bin_end_s"] - df["time_bin_start_s"]) / 3600.0,
        "height": max(plot_df[ylabel].max(), 1.0) * 2.0,
        "elapsed time (hours)": df["time_bin_end_s"] / 3600.0,
        ylabel: plot_df[ylabel].iloc[:len(df)].to_numpy(),
    })
    bokeh_relational.add_vertical_hover_bands(
        plot, hover_df, x="x", y="y", width="width", height="height")
    if ylabel == "Read count":
        bokeh_seqsummary.apply_scientific_count_formatter(plot._fig.yaxis[0])
    elif ylabel == "Base yield / Gbp":
        bokeh_seqsummary.apply_fixed_precision_formatter(plot._fig.yaxis[0], "0,0.00")
    bokeh_null.attach_null_state(plot)
    return plot


def quantile_ribbonplot(df, prefix, ylabel):
    """Create a ribbon plot from Q05-Q95 values."""
    null_messages = {
        "length": "No read length through time data available",
        "qscore": "No quality through time data available",
        "speed": "No translocation speed through time data available",
        "accuracy": "No alignment concordance through time data available",
    }
    x = "elapsed time (hours)"
    y = ylabel
    cols = {
        "time_end_h": x,
        "read_count": "read_count",
        f"{prefix}_plot_q05": "q05",
        f"{prefix}_plot_q25": "q25",
        f"{prefix}_plot_q50": y,
        f"{prefix}_plot_q75": "q75",
        f"{prefix}_plot_q95": "q95",
    }
    plot_df = df[list(cols)].rename(columns=cols)
    plot_df = plot_df.replace([np.inf, -np.inf], np.nan).dropna(
        subset=[x, "q05", "q25", y, "q75", "q95"])
    if plot_df.empty:
        return bokeh_null.null_plot(
            null_messages.get(prefix, f"No {ylabel.lower()} available"),
            "Elapsed time / hours",
            ylabel)
    if not plot_df.empty and plot_df.iloc[0][x] > 0.0:
        first = plot_df.iloc[[0]].copy()
        first[x] = 0.0
        plot_df = pd.concat([first, plot_df], ignore_index=True)
    plot = bokeh_ribbon.ribbonplot(
        plot_df,
        x=x,
        y=y,
        lower="q05",
        q1="q25",
        q2=y,
        q3="q75",
        upper="q95",
        linecolor="#222222",
        color=FASTCAT_BLUE,
        hover_items=[
            ("Elapsed time / hours", f"@{{{x}}}{{0.00}}"),
            ("Read count", "@read_count{0}"),
            ("Q05", "@q05{0.00}"),
            ("Q25", "@q25{0.00}"),
            ("Q50", f"@{{{y}}}{{0.00}}"),
            ("Q75", "@q75{0.00}"),
            ("Q95", "@q95{0.00}"),
        ])
    if plot is not None:
        plot._fig.xaxis.axis_label = "Elapsed time / hours"
        plot._fig.yaxis.axis_label = ylabel
    return plot


def add_time_bin_plots(path):
    """Add time-bin plots for one sample to the current report section."""
    df = load_time_bin_summary(path)

    tabs = Tabs()
    with tabs.add_tab("Yield"):
        EZChart(
            time_lineplot(
                df, "cumulative_base_yield_gbp",
                "Cumulative base yield / Gbp",
                x="time_end_h",
                seed_origin=True),
            height="360px")
        EZChart(time_stepplot(df, "base_yield_gbp", "Base yield / Gbp"), height="360px")
    with tabs.add_tab("Reads"):
        EZChart(
            time_lineplot(
                df, "cumulative_read_count",
                "Cumulative read count",
                color=FASTCAT_BLUE,
                x="time_end_h",
                seed_origin=True),
            height="360px")
        EZChart(time_stepplot(df, "read_count", "Read count"), height="360px")
    with tabs.add_tab("Read length"):
        length_plot = quantile_ribbonplot(df, "length", "Read length / kb")
        EZChart(length_plot, height="420px")
    with tabs.add_tab("Quality"):
        qscore_plot = quantile_ribbonplot(df, "qscore", "Quality")
        EZChart(qscore_plot, height="420px")
    with tabs.add_tab("Translocation speed"):
        speed_plot = quantile_ribbonplot(
            df, "speed", "Translocation speed / bp/s")
        EZChart(speed_plot, height="420px")
    with tabs.add_tab("Alignment Concordance"):
        accuracy_plot = quantile_ribbonplot(
            df, "accuracy", "Alignment Concordance / %")
        accuracy_plot._fig.y_range.start = 80
        accuracy_plot._fig.y_range.end = 100
        EZChart(accuracy_plot, height="420px")


def add_channel_summary_plots(path):
    """Add channel occupancy plots for one sample."""
    reads_plot = load_channel_summary_heatmap(path, "read_count")
    mean_length_plot = load_channel_summary_heatmap(path, "mean_read_length")
    channel_summary = pd.read_csv(path, sep="\t")
    has_quality = "mean_quality" in channel_summary.columns
    has_concordance = (
        "mean_concordance_qscore" in channel_summary.columns
        and pd.to_numeric(
            channel_summary["mean_concordance_qscore"], errors="coerce").notna().any())
    tabs = Tabs()
    with tabs.add_tab("Read counts"):
        EZChart(
            reads_plot,
            height="auto",
            additional_styles=aspect_matched_chart_styles(reads_plot))
    with tabs.add_tab("Mean read length"):
        EZChart(
            mean_length_plot,
            height="auto",
            additional_styles=aspect_matched_chart_styles(mean_length_plot))
    if has_quality:
        mean_quality_plot = load_channel_summary_heatmap(path, "mean_quality")
    else:
        mean_quality_plot = None
    if has_concordance:
        mean_concordance_plot = load_channel_summary_heatmap(path, "mean_concordance_qscore")
    else:
        mean_concordance_plot = None
    mean_metric_plots = [
        plot for plot in (mean_quality_plot, mean_concordance_plot)
        if plot is not None]
    if mean_metric_plots:
        mappers = [plot._fig.select_one(LinearColorMapper) for plot in mean_metric_plots]
        shared_high = np.ceil(max(mapper.high for mapper in mappers) / 5.0) * 5.0
        for mapper in mappers:
            mapper.high = shared_high
    if mean_quality_plot is not None:
        with tabs.add_tab("Mean quality"):
            EZChart(
                mean_quality_plot,
                height="auto",
                additional_styles=aspect_matched_chart_styles(mean_quality_plot))
    if mean_concordance_plot is not None:
        with tabs.add_tab("Mean alignment concordance"):
            EZChart(
                mean_concordance_plot,
                height="auto",
                additional_styles=aspect_matched_chart_styles(mean_concordance_plot))


def add_quality_diagnostics_plots(length_qscore_path, speed_qscore_path, accuracy_qscore_path=None):
    """Add quality diagnostic heatmaps to a shared tab set."""
    tabs = Tabs()
    with tabs.add_tab("Translocation speed"):
        EZChart(load_speed_qscore_heatmap(speed_qscore_path), height="650px")
    with tabs.add_tab("Read length"):
        EZChart(load_length_qscore_heatmap(length_qscore_path), height="650px")
    if accuracy_qscore_path:
        with tabs.add_tab("Alignment concordance"):
            plot = load_accuracy_qscore_heatmap(accuracy_qscore_path)
            EZChart(
                plot,
                height="auto",
                additional_styles=accuracy_qscore_chart_styles(plot))


def main(args):
    """Run the entry point."""
    get_main_logger("workflow-glue")
    logger = get_named_logger("Report")
    grouped = args.analysis_group is not None
    logger.info(
        f"Creating report for analysis group '{args.analysis_group}'."
        if grouped
        else "Creating report for all samples."
    )
    report_title = "fastcat sequencing report"
    if grouped:
        report_title += f" for analysis group '{args.analysis_group}'"
    report = labs.LabsReport(
        report_title, "fastcat",
        args.params, args.versions, args.wf_version,
        created_date="__FASTCAT_REPORT_DATE_HTML__")

    with open(args.metadata, "r") as f:
        metadata = json.load(f)
        sample_details = sorted([{
            'sample': d['alias'],
            'barcode': d['barcode']
        } for d in metadata], key=lambda x: x["sample"])

    sample_names = [d['sample'] for d in sample_details]
    samples_with_stats = tuple(d["alias"] for d in metadata if d["has_stats"])

    if grouped:
        with report.add_section(
            f"Analysis group report: '{args.analysis_group}'", "Intro"
        ):
            p(
                "This report contains all samples of the ",
                b(args.analysis_group),
                " analysis group: ",
                ", ".join(sample_names),
                "."
            )

    client_fields = None
    if args.client_fields:
        with open(args.client_fields) as f:
            try:
                client_fields = json.load(f)
            except json.decoder.JSONDecodeError:
                error = "ERROR: Client info is not correctly formatted"

        with report.add_section("Workflow Metadata", "Workflow Metadata"):
            if client_fields:
                df = pd.DataFrame.from_dict(
                    client_fields, orient="index", columns=["Value"])
                df.index.name = "Key"

                # Examples from the client had lists as values so join lists
                # for better display
                df['Value'] = df.Value.apply(
                    lambda x: ', '.join(
                        [str(i) for i in x]) if isinstance(x, list) else x)

                DataTable.from_pandas(df)
            else:
                p(error)

    with report.add_section("Run summary", "Run summary"):
        if args.run_summary:
            simple_kv_table(load_run_summary(args.run_summary[0]))
        else:
            detail = sample_details[0]
            simple_kv_table([
                ("Sample ID", detail["sample"]),
                ("Barcode ID", detail["barcode"]),
            ])

    if args.stats:
        with report.add_section("Run statistics", "Run statistics"):
            stats = tuple(args.stats)
            if len(stats) == 1:
                stats = stats[0]
                samples_with_stats = samples_with_stats[0]
            bokeh_seqsummary.SeqSummary(stats, sample_names=samples_with_stats)

    if args.alignment_stats and args.flagstats:
        with report.add_section("Alignment statistics", "Alignment statistics"):
            p(
                "Read-level mapping statistics. Alignment distributions include "
                "primary QC-passing mapped reads.")
            add_alignment_summary(args.flagstats)
            add_alignment_plots(args.alignment_stats)

    if args.coverage_summary and args.coverage_distribution:
        with report.add_section("Depth of coverage", "Depth of coverage"):
            p(
                "The cumulative depth distribution shows the percentage of reference "
                "positions covered at or above each depth for the selected target set.")
            add_coverage_plot(args.coverage_summary, args.coverage_distribution)

    if args.time_bin_summary:
        time_samples = samples_with_stats
        if isinstance(time_samples, str):
            time_samples = (time_samples,)

        with report.add_section("Sequencing over time", "Sequencing over time"):
            if len(args.time_bin_summary) == 1:
                add_time_bin_plots(args.time_bin_summary[0])
            else:
                tabs = Tabs()
                for sample, path in zip(time_samples, args.time_bin_summary):
                    with tabs.add_tab(sample):
                        add_time_bin_plots(path)

    if args.length_qscore_2d and args.speed_qscore_2d:
        heatmap_samples = samples_with_stats
        if isinstance(heatmap_samples, str):
            heatmap_samples = (heatmap_samples,)

        with report.add_section("Quality diagnostics", "Quality diagnostics"):
            if len(args.length_qscore_2d) == 1 and len(args.speed_qscore_2d) == 1:
                add_quality_diagnostics_plots(
                    args.length_qscore_2d[0],
                    args.speed_qscore_2d[0],
                    args.accuracy_qscore_2d[0] if args.accuracy_qscore_2d else None)
            else:
                tabs = Tabs()
                for index, (sample, length_path, speed_path) in enumerate(zip(
                        heatmap_samples,
                        args.length_qscore_2d,
                        args.speed_qscore_2d)):
                    with tabs.add_tab(sample):
                        accuracy_path = (
                            args.accuracy_qscore_2d[index]
                            if args.accuracy_qscore_2d else None)
                        add_quality_diagnostics_plots(length_path, speed_path, accuracy_path)

    if args.channel_summary:
        channel_samples = samples_with_stats
        if isinstance(channel_samples, str):
            channel_samples = (channel_samples,)

        with report.add_section("Channel occupancy", "Channel occupancy"):
            if len(args.channel_summary) == 1:
                add_channel_summary_plots(args.channel_summary[0])
            else:
                tabs = Tabs()
                for sample, path in zip(channel_samples, args.channel_summary):
                    with tabs.add_tab(sample):
                        add_channel_summary_plots(path)

    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument("--analysis_group", help="Analysis group name.")
    parser.add_argument(
        "--stats", nargs='+',
        help="Fastcat histogram directory or per-read stats, ordered as per entries in --metadata.")
    parser.add_argument(
        "--metadata", required=True,
        help="sample metadata")
    parser.add_argument(
        "--run-summary", nargs="+",
        help="Fastcat run_summary.tsv files, ordered as per entries in --metadata.")
    parser.add_argument(
        "--length-qscore-2d", nargs="+",
        help="Fastcat length_qscore_2d.tsv files, ordered as per entries in --metadata.")
    parser.add_argument(
        "--speed-qscore-2d", nargs="+",
        help="Fastcat speed_qscore_2d.tsv files, ordered as per entries in --metadata.")
    parser.add_argument(
        "--accuracy-qscore-2d", nargs="+",
        help="Fastcat accuracy_qscore_2d.tsv files, ordered as per entries in --metadata.")
    parser.add_argument(
        "--time-bin-summary", nargs="+",
        help="Fastcat time_bin_summary.tsv files, ordered as per entries in --metadata.")
    parser.add_argument(
        "--channel-summary", nargs="+",
        help="Fastcat channel_summary.tsv files, ordered as per entries in --metadata.")
    parser.add_argument(
        "--alignment-stats",
        help="Fastcat BAM histogram directory containing accuracy and coverage histograms.")
    parser.add_argument(
        "--flagstats",
        help="Fastcat flagstats.tsv used to scaffold alignment headline values.")
    parser.add_argument(
        "--coverage-summary", nargs="+",
        help="Coverage summary tables ordered to match --coverage-distribution.")
    parser.add_argument(
        "--coverage-distribution", nargs="+",
        help="Aggregate coverage distribution tables ordered to match --coverage-summary.")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument(
        "--client_fields",
        help="A JSON file containing useful key/values for display")
    parser.add_argument(
        "--wf_version", default='unknown',
        help="version of the executed workflow")
    return parser
