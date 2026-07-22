"""Generate the embedded summary-JSON schema C source from its JSON document."""

import argparse
import gzip
from pathlib import Path


def deterministic_gzip(data):
    """Compress data without Python-version-specific gzip header metadata."""
    compressed = bytearray(gzip.compress(data, compresslevel=9, mtime=0))
    compressed[9] = 255  # gzip OS field: unknown, rather than the host platform
    return bytes(compressed)


def c_byte_array(data):
    """Return a compact, readable C byte-array initializer."""
    lines = []
    for start in range(0, len(data), 12):
        chunk = ", ".join(f"0x{byte:02x}" for byte in data[start:start + 12])
        lines.append(f"    {chunk},")
    return "\n".join(lines)


def argparser():
    parser = argparse.ArgumentParser(
        add_help=False,
        description="Generate embedded C sources for the summary JSON Schema.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--header", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    return parser


def main(args):
    data = args.input.read_bytes()
    if not data.endswith(b"\n"):
        data += b"\n"
    compressed = deterministic_gzip(data)

    header = """#ifndef FASTCAT_SUMMARY_JSON_SCHEMA_DATA_H
#define FASTCAT_SUMMARY_JSON_SCHEMA_DATA_H

#include <stdio.h>

int fc_summary_json_schema_write(FILE *stream);

#endif
"""
    source = """/* generated from test/fixtures/misc/summary_json.schema.json; do not edit */
#include "cli/summary_json_schema_data.h"

#include <zlib.h>

static const unsigned char summary_json_schema_gz[] = {
"""
    source += c_byte_array(compressed)
    source += "\n};\n\n"
    source += """int fc_summary_json_schema_write(FILE *stream) {
  z_stream zs = {0};
  zs.next_in = (Bytef *)summary_json_schema_gz;
  zs.avail_in = (uInt)sizeof(summary_json_schema_gz);
  if (inflateInit2(&zs, 16 + MAX_WBITS) != Z_OK) {
    return -1;
  }

  int result = -1;
  int status = Z_OK;
  unsigned char buffer[8192];
  while (status == Z_OK) {
    zs.next_out = buffer;
    zs.avail_out = sizeof(buffer);
    status = inflate(&zs, Z_NO_FLUSH);
    if (status != Z_OK && status != Z_STREAM_END) {
      break;
    }
    size_t produced = sizeof(buffer) - zs.avail_out;
    if (produced > 0 && fwrite(buffer, 1, produced, stream) != produced) {
      break;
    }
  }
  if (status == Z_STREAM_END) {
    result = 0;
  }
  inflateEnd(&zs);
  return result;
}
"""

    args.header.write_text(header, encoding="utf-8")
    args.source.write_text(source, encoding="utf-8")
