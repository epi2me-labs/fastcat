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
        overrides = {
            "sample_id": alias,
            "input_type": TEMPLATE_INPUT_TYPE,
            "barcode_id": barcode,
            "barcode_alias": TEMPLATE_BARCODE_ALIAS,
            "run_ids": TEMPLATE_RUN_ID,
            "basecaller_models": TEMPLATE_BASECALLER,
            "modbase_models": TEMPLATE_MODBASE,
            "flow_cell_ids": TEMPLATE_FLOW_CELL_ID,
            "run_date": TEMPLATE_RUN_DATE,
            "run_duration_s": TEMPLATE_RUN_DURATION,
            "read_count": TEMPLATE_READ_COUNT,
            "base_count": TEMPLATE_YIELD,
            "records_seen": "__FASTCAT_RECORDS_SEEN_HTML__",
            "records_too_long": "__FASTCAT_RECORDS_TOO_LONG_HTML__",
            "records_too_short": "__FASTCAT_RECORDS_TOO_SHORT_HTML__",
            "records_low_quality": "__FASTCAT_RECORDS_LOW_QUALITY_HTML__",
            "records_dust_masked": "__FASTCAT_RECORDS_DUST_MASKED_HTML__",
            "n10_bp": TEMPLATE_N10,
            "n50_bp": TEMPLATE_N50,
            "n90_bp": TEMPLATE_N90,
            "l10": TEMPLATE_L10,
            "l50": TEMPLATE_L50,
            "l90": TEMPLATE_L90,
        }
        (out_dir / "run_summary_overrides.json").write_text(
            json.dumps(overrides, indent=2) + "\n", encoding="utf-8")
    (out_dir / "params.json").write_text(
        json.dumps(params, indent=2) + "\n", encoding="utf-8")
    (versions_dir / "versions.txt").write_text(
        f"fastcat,{version}\n"
        f"workflow,{version}\n",
        encoding="utf-8",
    )
    return 0
