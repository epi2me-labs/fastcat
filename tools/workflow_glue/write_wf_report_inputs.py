"""Write wf-template-style report inputs for the fastcat report scaffold."""

import argparse
import json
from pathlib import Path


TEMPLATE_ALIAS = "__FASTCAT_SAMPLE_ALIAS_HTML__"
TEMPLATE_BARCODE = "__FASTCAT_SAMPLE_BARCODE_HTML__"
TEMPLATE_BARCODE_ALIAS = "__FASTCAT_BARCODE_ALIAS_HTML__"
TEMPLATE_RUN_ID = "__FASTCAT_RUN_ID_HTML__"
TEMPLATE_BASECALLER = "__FASTCAT_BASECALLER_HTML__"
TEMPLATE_FLOW_CELL_ID = "__FASTCAT_FLOW_CELL_ID_HTML__"
TEMPLATE_RUN_DATE = "__FASTCAT_RUN_DATE_HTML__"
TEMPLATE_RUN_DURATION = "__FASTCAT_RUN_DURATION_HTML__"
TEMPLATE_YIELD = "__FASTCAT_YIELD_HTML__"
TEMPLATE_READ_COUNT = "__FASTCAT_READ_COUNT_HTML__"
TEMPLATE_N50 = "__FASTCAT_N50_HTML__"
TEMPLATE_L50 = "__FASTCAT_L50_HTML__"
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
    (out_dir / "run_summary.tsv").write_text(
        "sample_id\tbarcode_id\tbarcode_alias\trun_id\tbasecaller_model\tflow_cell_id\trun_date\t"
        "run_duration_s\tread_count\tbase_count\tn50_bp\tl50\n"
        f"{alias}\t"
        f"{barcode}\t"
        f"{TEMPLATE_BARCODE_ALIAS if args.template else 'Synthetic reads'}\t"
        f"{TEMPLATE_RUN_ID if args.template else 'synthetic-run'}\t"
        f"{TEMPLATE_BASECALLER if args.template else 'dna_r10.4.1_e8.2_400bps_sup@v5.2.0'}\t"
        f"{TEMPLATE_FLOW_CELL_ID if args.template else 'synthetic-flowcell'}\t"
        f"{TEMPLATE_RUN_DATE if args.template else '2026-01-01 00:00:00'}\t"
        f"{TEMPLATE_RUN_DURATION if args.template else '172800'}\t"
        f"{TEMPLATE_READ_COUNT if args.template else '100000'}\t"
        f"{TEMPLATE_YIELD if args.template else '7620000000'}\t"
        f"{TEMPLATE_N50 if args.template else '3230'}\t"
        f"{TEMPLATE_L50 if args.template else '2063112'}\n",
        encoding="utf-8")
    (out_dir / "params.json").write_text(
        json.dumps(params, indent=2) + "\n", encoding="utf-8")
    (versions_dir / "versions.txt").write_text(
        f"fastcat,{version}\n"
        f"workflow,{version}\n",
        encoding="utf-8",
    )
    return 0
