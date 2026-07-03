"""Write wf-template-style report inputs for the fastcat report scaffold."""

import argparse
import json
from pathlib import Path


TEMPLATE_ALIAS = "__FASTCAT_SAMPLE_ALIAS_HTML__"
TEMPLATE_TYPE = "__FASTCAT_SAMPLE_TYPE_HTML__"
TEMPLATE_BARCODE = "__FASTCAT_SAMPLE_BARCODE_HTML__"
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
    parser.add_argument("--sample-type", default="synthetic", help="sample type")
    parser.add_argument("--barcode", default="unclassified", help="sample barcode")
    parser.add_argument("--version", default="dev", help="fastcat/report version")
    return parser


def main(args):
    out_dir = Path(args.out_dir)
    versions_dir = out_dir / "versions"
    versions_dir.mkdir(parents=True, exist_ok=True)

    alias = TEMPLATE_ALIAS if args.template else args.alias
    sample_type = TEMPLATE_TYPE if args.template else args.sample_type
    barcode = TEMPLATE_BARCODE if args.template else args.barcode
    version = TEMPLATE_VERSION if args.template else args.version

    metadata = [{
        "alias": alias,
        "type": sample_type,
        "barcode": barcode,
        "has_stats": True,
    }]
    params = {
        "input": "synthetic",
        "report_template": args.template,
    }

    (out_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (out_dir / "params.json").write_text(
        json.dumps(params, indent=2) + "\n", encoding="utf-8")
    (versions_dir / "versions.txt").write_text(
        f"fastcat,{version}\n"
        f"workflow,{version}\n",
        encoding="utf-8",
    )
    return 0
