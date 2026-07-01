"""Reads SAM from stdin and writes to STDOUT, sorting auxiliary tags by key."""

import argparse
import sys


def argparser():
    return argparse.ArgumentParser(
        description="Read SAM from stdin and write records with auxiliary tags sorted by key.",
        add_help=False,
    )


def main(args):
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 11:
            continue
        core_fields = fields[:11]
        aux_fields = fields[11:]
        aux_fields.sort()
        print("\t".join(core_fields + aux_fields))
