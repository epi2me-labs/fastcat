"""Write wf-template-style report inputs for the fastcat report scaffold."""

import argparse
import json
from pathlib import Path


TEMPLATE_ALIAS = "__FASTCAT_SAMPLE_ALIAS_HTML__"
TEMPLATE_INPUT_TYPE = "__FASTCAT_INPUT_TYPE_HTML__"
TEMPLATE_BARCODE = "__FASTCAT_SAMPLE_BARCODE_HTML__"
TEMPLATE_BARCODE_ALIAS = "__FASTCAT_BARCODE_ALIAS_HTML__"
TEMPLATE_RUN_ID = "__FASTCAT_RUN_ID_HTML__"
TEMPLATE_BASECALLER = "__FASTCAT_BASECALLER_HTML__"
TEMPLATE_MODBASE = "__FASTCAT_MODBASE_HTML__"
TEMPLATE_FLOW_CELL_ID = "__FASTCAT_FLOW_CELL_ID_HTML__"
TEMPLATE_RUN_DATE = "__FASTCAT_RUN_DATE_HTML__"
TEMPLATE_RUN_DURATION = "__FASTCAT_RUN_DURATION_HTML__"
TEMPLATE_YIELD = "__FASTCAT_YIELD_HTML__"
TEMPLATE_READ_COUNT = "__FASTCAT_READ_COUNT_HTML__"
TEMPLATE_N50 = "__FASTCAT_N50_HTML__"
TEMPLATE_N10 = "__FASTCAT_N10_HTML__"
TEMPLATE_N90 = "__FASTCAT_N90_HTML__"
TEMPLATE_L50 = "__FASTCAT_L50_HTML__"
TEMPLATE_L10 = "__FASTCAT_L10_HTML__"
TEMPLATE_L90 = "__FASTCAT_L90_HTML__"
TEMPLATE_VERSION = "__FASTCAT_VERSION_HTML__"


def argparser():
    parser = argparse.ArgumentParser(
        description="write report scaffold inputs",
        add_help=False)
    parser.add_argument("--out-dir", required=True, help="output input directory")
    parser.add_argument(
        "--template", action="store_true",
        help="write placeholder values used when vendoring the report template")
    parser.add_argument("--alias", default="Synthetic reads", help="sample alias")
    parser.add_argument("--barcode", default="unclassified", help="sample barcode")
    parser.add_argument("--version", default="dev", help="fastcat/report version")
    return parser


def main(args):
    out_dir = Path(args.out_dir)
    versions_dir = out_dir / "versions"
    versions_dir.mkdir(parents=True, exist_ok=True)

    alias = TEMPLATE_ALIAS if args.template else args.alias
    barcode = TEMPLATE_BARCODE if args.template else args.barcode
    version = TEMPLATE_VERSION if args.template else args.version

    metadata = [{
        "alias": alias,
        "barcode": barcode,
        "has_stats": True,
    }]
    params = {
        "input": "synthetic",
        "report_template": args.template,
    }

    (out_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    if args.template:
        (out_dir / "run_summary.tsv").write_text(
            "sample_id\tinput_type\tbarcode_id\tbarcode_alias\trun_id\tbasecaller_model\tmodbase_models\tflow_cell_id\trun_date\t"
            "run_duration_s\tread_count\tbase_count\t"
            "n5_bp\tn10_bp\tn25_bp\tn50_bp\tn75_bp\tn90_bp\tn95_bp\t"
            "l5\tl10\tl25\tl50\tl75\tl90\tl95\n"
            f"{alias}\t{TEMPLATE_INPUT_TYPE}\t{barcode}\t{TEMPLATE_BARCODE_ALIAS}\t{TEMPLATE_RUN_ID}\t"
            f"{TEMPLATE_BASECALLER}\t{TEMPLATE_MODBASE}\t{TEMPLATE_FLOW_CELL_ID}\t{TEMPLATE_RUN_DATE}\t"
            f"{TEMPLATE_RUN_DURATION}\t{TEMPLATE_READ_COUNT}\t{TEMPLATE_YIELD}\t"
            f"\t{TEMPLATE_N10}\t\t{TEMPLATE_N50}\t\t{TEMPLATE_N90}\t\t"
            f"\t{TEMPLATE_L10}\t\t{TEMPLATE_L50}\t\t{TEMPLATE_L90}\t\n",
            encoding="utf-8")
    (out_dir / "params.json").write_text(
        json.dumps(params, indent=2) + "\n", encoding="utf-8")
    (versions_dir / "versions.txt").write_text(
        f"fastcat,{version}\n"
        f"workflow,{version}\n",
        encoding="utf-8",
    )
    return 0
