"""Create workflow report."""
import json

from dominate.tags import b, p
from ezcharts.components import fastcat
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Tabs
from ezcharts.layout.snippets.table import DataTable
import pandas as pd

from .util import get_main_logger, get_named_logger, wf_parser  # noqa: ABS101


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
