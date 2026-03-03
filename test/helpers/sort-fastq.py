#!/usr/bin/env python3
"""Read FASTQ from stdin and write records sorted by read id."""
import sys


def iter_records(stream):
    while True:
        header = stream.readline()
        if not header:
            return
        seq = stream.readline()
        plus = stream.readline()
        qual = stream.readline()
        if not qual:
            return
        yield (header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n"))


def record_key(rec):
    header = rec[0]
    read_id = header[1:].split()[0] if header.startswith("@") else header.split()[0]
    return (read_id, header, rec[1], rec[3])


def main():
    records = list(iter_records(sys.stdin))
    records.sort(key=record_key)
    for header, seq, plus, qual in records:
        sys.stdout.write(f"{header}\n{seq}\n{plus}\n{qual}\n")


if __name__ == "__main__":
    main()
