"""Validate fastcat structured JSON summaries against their JSON Schema."""

import argparse
import base64
import copy
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from jsonschema import Draft202012Validator


DTYPE_WIDTHS = {
    "uint8": 1,
    "uint16": 2,
    "uint32": 4,
    "uint64": 8,
    "float64": 8,
}

DATAFRAME_PATHS = (
    "$.channels.__data__",
    "$.distributions.read_length.__data__",
    "$.distributions.quality.__data__",
    "$.alignment.distributions.accuracy.__data__",
    "$.alignment.distributions.query_coverage.__data__",
    "$.alignment.distributions.unmapped_read_length.__data__",
    "$.alignment.distributions.unmapped_quality.__data__",
    "$.alignment.distributions.poly_a_length.__data__",
    "$.time.__data__",
)


def argparser():
    parser = argparse.ArgumentParser(
        description="Validate fastcat summary.json files.",
        add_help=False)
    parser.add_argument(
        "--schema",
        type=Path,
        required=True,
        help="JSON Schema document to validate against")
    parser.add_argument(
        "summary_json",
        type=Path,
        nargs="+",
        help="summary.json or summary.json.gz files to validate")
    parser.add_argument(
        "--expected",
        type=Path,
        help="expected summary.json to validate and compare with the single input")
    parser.add_argument(
        "--mask-barcode",
        action="store_true",
        help="ignore sample and run-summary barcode IDs during comparison")
    parser.add_argument(
        "--path-root",
        type=Path,
        default=Path.cwd(),
        help="make absolute file paths relative to this directory during comparison")
    return parser


def load_json(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return json.load(handle)


def decode_payload(data, path, element_count, dtype):
    """Decode an ndarray payload and verify its expected byte length."""
    try:
        payload = base64.b64decode(data["__ndarray__"], validate=True)
    except (ValueError, TypeError) as error:
        return None, f"{path}.__ndarray__: invalid base64 payload ({error})"

    expected_bytes = element_count * DTYPE_WIDTHS[dtype]
    if len(payload) != expected_bytes:
        return None, (f"{path}.__ndarray__: expected {expected_bytes} bytes for "
                      f"{element_count} elements and {dtype}, got {len(payload)}")
    return payload, None


def validate_ndarray(data, path):
    """Return integrity errors for one schema-validated ndarray object."""
    element_count = 1
    for dimension in data["shape"]:
        element_count *= dimension
    _, error = decode_payload(data, path, element_count, data["dtype"])
    return [] if error is None else [error]


def decode_dataframe(data, path):
    """Decode a dataframe block exactly as a pandas consumer would."""
    row_count = data["shape"][0]
    arrays = {}
    errors = []
    for index, column in enumerate(data["columns"]):
        column_path = f"{path}.columns[{index}]"
        payload, error = decode_payload(column, column_path, row_count, column["dtype"])
        if error is not None:
            errors.append(error)
            continue
        arrays[column["name"]] = np.frombuffer(payload, dtype=column["dtype"])

    if errors:
        return None, errors

    dataframe = pd.DataFrame(arrays)
    columns = data["columns"]
    expected_shape = (row_count, len(columns))
    errors = []
    if dataframe.shape != expected_shape:
        errors.append(f"{path}: pandas dataframe has shape {dataframe.shape}, expected {expected_shape}")
    decoded_dtypes = tuple(str(dtype) for dtype in dataframe.dtypes)
    expected_dtypes = tuple(column["dtype"] for column in columns)
    if decoded_dtypes != expected_dtypes:
        errors.append(
            f"{path}: pandas dataframe dtypes expected {list(expected_dtypes)}, "
            f"got {list(decoded_dtypes)}")
    return dataframe, errors


def validate_dataframe(data, path):
    """Validate a dataframe block can be decoded exactly as a pandas consumer would."""
    _, errors = decode_dataframe(data, path)
    return errors


def validate_document(document, schema):
    """Return all schema and binary-integrity errors for one summary document."""
    validator = Draft202012Validator(schema)
    errors = [f"{error.json_path}: {error.message}" for error in validator.iter_errors(document)]
    if errors:
        return sorted(errors)

    for path in DATAFRAME_PATHS:
        value = document
        for key in path.removeprefix("$.").split("."):
            if value is None or key not in value:
                value = None
                break
            value = value[key]
        if value is None:
            continue
        errors.extend(validate_dataframe(value, path))

    for name, heatmap in document["heatmaps"].items():
        if heatmap is not None:
            errors.extend(validate_ndarray(heatmap["__data__"], f"$.heatmaps.{name}.__data__"))
    return sorted(errors)


def normalize_paths(document, root):
    """Return a copy with repository-local file paths made relative."""
    normalized = copy.deepcopy(document)
    root = root.resolve()
    for entry in normalized.get("files", []):
        path = Path(entry.get("path", ""))
        if path.is_absolute():
            try:
                entry["path"] = str(path.resolve().relative_to(root))
            except ValueError:
                pass
    return normalized


def compare_documents(actual, expected, mask_barcode, path_root):
    """Return an equality error after normalizing intentional differences."""
    actual = normalize_paths(actual, path_root)
    expected = normalize_paths(expected, path_root)
    if mask_barcode:
        for document in (actual, expected):
            for key in ("sample", "run_summary"):
                if key in document:
                    document[key]["barcode_id"] = "__barcode__"
    if actual != expected:
        return "document differs from expected result"
    return None


def validate_path(path, schema, schema_path):
    """Load and validate one document, returning it on success."""
    document = load_json(path)
    errors = validate_document(document, schema)
    if errors:
        messages = "\n".join(f"  {error}" for error in errors)
        raise SystemExit(f"{path} does not match {schema_path}:\n{messages}")
    print(f"validated {path}")
    return document


def main(args):
    schema = load_json(args.schema)
    Draft202012Validator.check_schema(schema)

    if args.expected is not None and len(args.summary_json) != 1:
        raise SystemExit("--expected requires exactly one summary.json input")

    for path in args.summary_json:
        actual = validate_path(path, schema, args.schema)

    if args.expected is not None:
        expected = validate_path(args.expected, schema, args.schema)
        error = compare_documents(actual, expected, args.mask_barcode, args.path_root)
        if error:
            raise SystemExit(f"{path}: {error} ({args.expected})")
