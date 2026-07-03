"""Create workflow report."""
import json

from dominate.tags import b, p
import ezcharts as ezc
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Tabs
from ezcharts.layout.snippets.table import DataTable
import pandas as pd

from . import bokeh_boxplot
from . import bokeh_matrix
from .util import get_main_logger, get_named_logger, wf_parser  # noqa: ABS101


def load_length_qscore_heatmap(path):
    """Load fastcat length/qscore TSV into a Bokeh heatmap matrix."""
    df = pd.read_csv(path, sep="\t")
    df["qscore"] = (df["qscore_bin_start"] + df["qscore_bin_end"]) * 0.5
    df["read length (kb)"] = (df["length_bin_start"] + df["length_bin_end"]) * 0.0005

    data = df.pivot(
        index="qscore",
        columns="read length (kb)",
        values="count")
    data = data.rename_axis("qscore")
    data = data.rename_axis("read length (kb)", axis=1)

    return bokeh_matrix.heatmap(
        data,
        annot=False,
        vmax_quantile=0.995,
        cbar_kws={"title": "count"},
        subtitle="Estimated counts")


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
        "qscore_q75", "qscore_q95")
    for col in numeric_cols:
        if col in df:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.sort_values("time_bin_start_s").reset_index(drop=True)
    df["time_mid_h"] = (
        df["time_bin_start_s"] + df["time_bin_end_s"]) / 7200.0
    labels = [
        f"{format_hours(row.time_bin_start_s / 3600.0)}-"
        f"{format_hours(row.time_bin_end_s / 3600.0)}"
        for row in df.itertuples()
    ]
    df["time_bin"] = pd.Categorical(labels, categories=labels, ordered=True)
    df["base_yield_mbp"] = df["base_count"] / 1e6
    df["cumulative_base_yield_gbp"] = df["cumulative_base_count"] / 1e9
    df["length_mean_kb"] = df["length_mean"] / 1000.0
    for col in ("q05", "q25", "q50", "q75", "q95"):
        df[f"length_plot_{col}"] = df[f"length_{col}"] / 1000.0
        df[f"qscore_plot_{col}"] = df[f"qscore_{col}"]
    return df


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


def quantile_boxplot(df, prefix, ylabel):
    """Create a box-style plot from KLL-estimated Q05-Q95 values."""
    x = "time bin (hours)"
    y = ylabel
    cols = {
        "time_bin": x,
        "read_count": "read_count",
        f"{prefix}_plot_q05": "q05",
        f"{prefix}_plot_q25": "q25",
        f"{prefix}_plot_q50": y,
        f"{prefix}_plot_q75": "q75",
        f"{prefix}_plot_q95": "q95",
    }
    plot_df = df[list(cols)].rename(columns=cols)
    return bokeh_boxplot.boxplot(
        plot_df,
        x=x,
        y=y,
        lower="q05",
        q1="q25",
        q2=y,
        q3="q75",
        upper="q95",
        linecolor="#222222",
        color="#2f7f6f",
        hover_items=[
            ("time bin", f"@{{{x}}}"),
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
        p("KLL-estimated read length quartile boxes by time bin.")
        length_plot = quantile_boxplot(df, "length", "read length (kb)")
        if length_plot is None:
            p("No read length quartiles were available.")
        else:
            EZChart(length_plot, height="420px")
    with tabs.add_tab("Qscore"):
        p("KLL-estimated qscore quartile boxes by time bin.")
        qscore_plot = quantile_boxplot(df, "qscore", "qscore")
        if qscore_plot is None:
            p("No qscore quartiles were available.")
        else:
            EZChart(qscore_plot, height="420px")


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

    if args.length_qscore_2d:
        heatmap_samples = samples_with_stats
        if isinstance(heatmap_samples, str):
            heatmap_samples = (heatmap_samples,)

        with report.add_section("Read length vs qscore", "Read length vs qscore"):
            if len(args.length_qscore_2d) == 1:
                EZChart(load_length_qscore_heatmap(args.length_qscore_2d[0]), height="650px")
            else:
                tabs = Tabs()
                for sample, path in zip(heatmap_samples, args.length_qscore_2d):
                    with tabs.add_tab(sample):
                        EZChart(load_length_qscore_heatmap(path), height="650px")

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
        "--time-bin-summary", nargs="+",
        help="Fastcat time_bin_summary.tsv files, ordered as per entries in --metadata.")
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


if __name__ == "__main__":
    raise SystemExit(main(argparser().parse_args()))
