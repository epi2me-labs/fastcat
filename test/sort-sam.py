#!/usr/bin/env python3
"""Reads SAM from stdin and writes to STDOUT, sorting auxiliary tags by key."""
import sys
import re

if __name__ == "__main__":
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        else:
            fields = line.strip().split("\t")
            if len(fields) < 11:
                continue
            core_fields = fields[:11]
            aux_fields = fields[11:]
            aux_fields.sort()
            print("\t".join(core_fields + aux_fields))
