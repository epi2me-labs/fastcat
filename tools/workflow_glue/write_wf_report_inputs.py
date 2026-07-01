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
    return parser


def main(args):
    out_dir = Path(args.out_dir)
    versions_dir = out_dir / "versions"
    versions_dir.mkdir(parents=True, exist_ok=True)

    metadata = [{
        "alias": TEMPLATE_ALIAS,
        "type": TEMPLATE_TYPE,
        "barcode": TEMPLATE_BARCODE,
        "has_stats": True,
    }]
    params = {
        "input": "synthetic",
        "report_template": True,
    }

    (out_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (out_dir / "params.json").write_text(
        json.dumps(params, indent=2) + "\n", encoding="utf-8")
    (versions_dir / "versions.txt").write_text(
        f"fastcat,{TEMPLATE_VERSION}\n"
        f"workflow,{TEMPLATE_VERSION}\n",
        encoding="utf-8",
    )
    return 0