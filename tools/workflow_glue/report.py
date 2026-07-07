"""Create workflow report."""
import json
from pathlib import Path

from dominate.tags import b, p
import ezcharts as ezc
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Tabs
from ezcharts.layout.snippets.table import DataTable
import numpy as np
import pandas as pd

from . import bokeh_matrix
from . import bokeh_ribbon
from .util import get_main_logger, get_named_logger, wf_parser  # noqa: ABS101


CHANNEL_MAP_ROOT = Path(__file__).resolve().parent / "data"
CHANNEL_MAP_PATHS = {
    "flongle": CHANNEL_MAP_ROOT / "asicb512fx_flongle_channel_map.csv",
    "minion": CHANNEL_MAP_ROOT / "asicb512fx_channel_map.csv",
    "promethion": CHANNEL_MAP_ROOT / "3000_opx_channel_map.csv",
}


def load_dense_heatmap(df, x, y, value, subtitle):
    """Load a dense rectangular TSV into a Bokeh heatmap matrix."""
    data = df.pivot(
        index=x,
        columns=y,
        values=value)
    data = data.rename_axis(x)
    data = data.rename_axis(y, axis=1)

    return bokeh_matrix.heatmap(
        data,
        annot=False,
        vmax_quantile=0.995,
        cbar_kws={"title": "count"},
        subtitle=subtitle)


def load_length_qscore_heatmap(path):
    """Load fastcat length/qscore TSV into a Bokeh heatmap matrix."""
    df = pd.read_csv(path, sep="\t")
    df["read length (kb)"] = (df["length_bin_start"] + df["length_bin_end"]) * 0.0005
    df["qscore"] = (df["qscore_bin_start"] + df["qscore_bin_end"]) * 0.5
    return load_dense_heatmap(
        df,
        "read length (kb)",
        "qscore",
        "count",
        "Estimated counts")


def load_speed_qscore_heatmap(path):
    """Load fastcat qscore/translocation-speed TSV into a Bokeh heatmap matrix."""
    df = pd.read_csv(path, sep="\t")
    df["translocation speed (nt/s)"] = (
        df["speed_bin_start"] + df["speed_bin_end"]) * 0.5
    df["qscore"] = (df["qscore_bin_start"] + df["qscore_bin_end"]) * 0.5
    return load_dense_heatmap(
        df,
        "translocation speed (nt/s)",
        "qscore",
        "count",
        "Estimated counts")


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
        "translocation_speed_q75", "translocation_speed_q95")
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
    df["base_yield_mbp"] = df["base_count"] / 1e6
    df["cumulative_base_yield_gbp"] = df["cumulative_base_count"] / 1e9
    df["length_mean_kb"] = df["length_mean"] / 1000.0
    df["speed_mean"] = df["translocation_speed_mean"]
    for col in ("q05", "q25", "q50", "q75", "q95"):
        df[f"length_plot_{col}"] = df[f"length_{col}"] / 1000.0
        df[f"qscore_plot_{col}"] = df[f"qscore_{col}"]
        df[f"speed_plot_{col}"] = df[f"translocation_speed_{col}"]
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
        max_channel = 1
    else:
        df["channel"] = pd.to_numeric(df["channel"], errors="coerce").fillna(0).astype(int)
        df["mux"] = pd.to_numeric(df["mux"], errors="coerce").fillna(0).astype(int)
        df["read_count"] = pd.to_numeric(df["read_count"], errors="coerce").fillna(0.0)
        df["base_count"] = pd.to_numeric(df["base_count"], errors="coerce").fillna(0.0)
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
    merged["count"] = merged[metric]
    merged["width"] = 1.0
    merged["height"] = 1.0
    data = merged.pivot(index="x", columns="y", values=metric)
    data = data.reindex(
        index=np.arange(data.index.min(), data.index.max() + 1.0, 1.0),
        columns=np.arange(data.columns.min(), data.columns.max() + 1.0, 1.0))
    data = data.rename_axis(None)
    data = data.rename_axis(None, axis=1)

    subtitle = (
        f"{int(merged['read_count'].sum())} reads across "
        f"{int((merged['read_count'] > 0).sum())} occupied channel/mux cells")

    return bokeh_matrix.heatmap(
        data,
        annot=False,
        vmax_quantile=0.995,
        cbar_kws={"title": metric.replace("_", " ")},
        subtitle=subtitle,
        hide_axes=True,
        match_aspect=True,
        hide_toolbar=True,
        hover_data=merged[["x", "y", "width", "height", "channel", "mux", "count"]],
        hover_items=[
            ("channel", "@channel{0}"),
            ("mux", "@mux{0}"),
            ("count", "@count{0}"),
        ])


def time_lineplot(df, y, ylabel, color="#2f7f6f"):
    """Create a line plot against elapsed run time in hours."""
    plot = ezc.lineplot(
        data=df,
        x="time_mid_h",
        y=y,
        estimator=None,
        errorbar=None,
        color=color)
    plot._fig.xaxis.axis_label = "elapsed time (hours)"
    plot._fig.yaxis.axis_label = ylabel
    return plot


def quantile_ribbonplot(df, prefix, ylabel):
    """Create a ribbon plot from KLL-estimated Q05-Q95 values."""
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
    if not plot_df.empty and plot_df.iloc[0][x] > 0.0:
        first = plot_df.iloc[[0]].copy()
        first[x] = 0.0
        plot_df = pd.concat([first, plot_df], ignore_index=True)
    return bokeh_ribbon.ribbonplot(
        plot_df,
        x=x,
        y=y,
        lower="q05",
        q1="q25",
        q2=y,
        q3="q75",
        upper="q95",
        linecolor="#222222",
        color="#6aa6f8",
        hover_items=[
            ("time (h)", f"@{{{x}}}{{0.00}}"),
            ("reads", "@read_count{0}"),
            ("Q05", "@q05{0.00}"),
            ("Q25", "@q25{0.00}"),
            ("Q50", f"@{{{y}}}{{0.00}}"),
            ("Q75", "@q75{0.00}"),
            ("Q95", "@q95{0.00}"),
        ])


def add_time_bin_plots(path):
    """Add time-bin plots for one sample to the current report section."""
    df = load_time_bin_summary(path)
    if df.empty or df["read_count"].sum() == 0:
        p("No timestamped reads were available for time-based plots.")
        return

    tabs = Tabs()
    with tabs.add_tab("Yield"):
        p("Cumulative base output through the run.")
        EZChart(
            time_lineplot(
                df, "cumulative_base_yield_gbp",
                "cumulative base yield (Gbp)"),
            height="360px")
    with tabs.add_tab("Reads"):
        p("Cumulative accepted read output through the run.")
        EZChart(
            time_lineplot(
                df, "cumulative_read_count",
                "cumulative read count",
                color="#1f6f9f"),
            height="360px")
    with tabs.add_tab("Read length"):
        p("KLL-estimated read length quantile ribbons by time bin.")
        length_plot = quantile_ribbonplot(df, "length", "read length (kb)")
        if length_plot is None:
            p("No read length quartiles were available.")
        else:
            EZChart(length_plot, height="420px")
    with tabs.add_tab("Qscore"):
        p("KLL-estimated qscore quantile ribbons by time bin.")
        qscore_plot = quantile_ribbonplot(df, "qscore", "qscore")
        if qscore_plot is None:
            p("No qscore quartiles were available.")
        else:
            EZChart(qscore_plot, height="420px")
    with tabs.add_tab("Translocation speed"):
        p("KLL-estimated translocation speed quantile ribbons by time bin.")
        speed_plot = quantile_ribbonplot(
            df, "speed", "translocation speed (nt/s)")
        if speed_plot is None:
            p("No translocation speed quartiles were available.")
        else:
            EZChart(speed_plot, height="420px")


def add_channel_summary_plots(path):
    """Add channel occupancy plots for one sample."""
    reads_plot = load_channel_summary_heatmap(path, "read_count")
    bases_plot = load_channel_summary_heatmap(path, "base_count")
    reads_aspect = f"{reads_plot._fig.width} / {reads_plot._fig.height}"
    bases_aspect = f"{bases_plot._fig.width} / {bases_plot._fig.height}"
    reads_max_width = f"{int(round(500 * reads_plot._fig.width / reads_plot._fig.height))}px"
    bases_max_width = f"{int(round(500 * bases_plot._fig.width / bases_plot._fig.height))}px"
    tabs = Tabs()
    with tabs.add_tab("Read counts"):
        EZChart(
            reads_plot,
            height="auto",
            additional_styles={
                "aspect-ratio": reads_aspect,
                "max-width": reads_max_width,
                "max-height": "500px",
                "margin-left": "auto",
                "margin-right": "auto",
            })
    with tabs.add_tab("Base counts"):
        EZChart(
            bases_plot,
            height="auto",
            additional_styles={
                "aspect-ratio": bases_aspect,
                "max-width": bases_max_width,
                "max-height": "500px",
                "margin-left": "auto",
                "margin-right": "auto",
            })


def add_quality_diagnostics_plots(length_qscore_path, speed_qscore_path):
    """Add the two quality diagnostic heatmaps to a shared tab set."""
    tabs = Tabs()
    with tabs.add_tab("Read length"):
        EZChart(load_length_qscore_heatmap(length_qscore_path), height="650px")
    with tabs.add_tab("Translocation speed"):
        EZChart(load_speed_qscore_heatmap(speed_qscore_path), height="650px")


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
        args.params, args.versions, args.wf_version)

    with open(args.metadata, "r") as f:
        metadata = json.load(f)
        sample_details = sorted([{
            'sample': d['alias'],
            'type': d['type'],
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

    if args.stats:
        with report.add_section("Read summary", "Read summary"):
            stats = tuple(args.stats)
            if len(stats) == 1:
                stats = stats[0]
                samples_with_stats = samples_with_stats[0]
            fastcat.SeqSummary(stats, sample_names=samples_with_stats)

    if args.length_qscore_2d and args.speed_qscore_2d:
        heatmap_samples = samples_with_stats
        if isinstance(heatmap_samples, str):
            heatmap_samples = (heatmap_samples,)

        with report.add_section("Quality diagnostics", "Quality diagnostics"):
            if len(args.length_qscore_2d) == 1 and len(args.speed_qscore_2d) == 1:
                add_quality_diagnostics_plots(
                    args.length_qscore_2d[0],
                    args.speed_qscore_2d[0])
            else:
                tabs = Tabs()
                for sample, length_path, speed_path in zip(
                        heatmap_samples,
                        args.length_qscore_2d,
                        args.speed_qscore_2d):
                    with tabs.add_tab(sample):
                        add_quality_diagnostics_plots(length_path, speed_path)

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

    with report.add_section("Sample Metadata", "Sample Metadata"):
        tabs = Tabs()
        for d in sample_details:
            with tabs.add_tab(d["sample"]):
                df = pd.DataFrame.from_dict(d, orient="index", columns=["Value"])
                df.index.name = "Key"
                DataTable.from_pandas(df)

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
        "--length-qscore-2d", nargs="+",
        help="Fastcat length_qscore_2d.tsv files, ordered as per entries in --metadata.")
    parser.add_argument(
        "--speed-qscore-2d", nargs="+",
        help="Fastcat speed_qscore_2d.tsv files, ordered as per entries in --metadata.")
    parser.add_argument(
        "--time-bin-summary", nargs="+",
        help="Fastcat time_bin_summary.tsv files, ordered as per entries in --metadata.")
    parser.add_argument(
        "--channel-summary", nargs="+",
        help="Fastcat channel_summary.tsv files, ordered as per entries in --metadata.")
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
