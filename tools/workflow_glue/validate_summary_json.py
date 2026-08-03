"""Validate fastcat structured JSON summaries against their JSON Schema."""

import argparse
import base64
import copy
from pathlib import Path

import numpy as np
import pandas as pd
from jsonschema import Draft202012Validator

from . import summary_json


DATAFRAME_PATHS = (
    "$.channels.data",
    "$.distributions.read_length.data",
    "$.distributions.quality.data",
    "$.alignment.distributions.accuracy.data",
    "$.alignment.distributions.query_coverage.data",
    "$.alignment.distributions.poly_a_length.data",
    "$.time.data",
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


def decode_payload(data, path, element_count, dtype):
    """Decode an ndarray payload and verify its expected byte length."""
    try:
        payload = base64.b64decode(data["data"], validate=True)
    except (ValueError, TypeError) as error:
        return None, f"{path}.data: invalid base64 payload ({error})"

    expected_bytes = element_count * summary_json.DTYPE_WIDTHS[dtype]
    if len(payload) != expected_bytes:
        return None, (f"{path}.data: expected {expected_bytes} bytes for "
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
    row_count, column_count = data["shape"]
    arrays = {}
    errors = []
    for index, column in enumerate(data["columns"]):
        column_path = f"{path}.columns[{index}]"
        try:
            arrays[column["name"]] = summary_json.decode_dataframe({
                "kind": "dataframe",
                "shape": [row_count, 1],
                "columns": [column],
            }).iloc[:, 0]
        except (KeyError, TypeError, ValueError) as error:
            errors.append(f"{column_path}: {error}")

    if errors:
        return None, errors

    dataframe = pd.DataFrame(arrays)
    columns = data["columns"]
    expected_shape = (row_count, column_count)
    errors = []
    if dataframe.shape != expected_shape:
        errors.append(f"{path}: pandas dataframe has shape {dataframe.shape}, expected {expected_shape}")
    decoded_dtypes = tuple(str(dtype) for dtype in dataframe.dtypes)
    expected_dtypes = tuple(
        "string" if column["dtype"] == "string" else str(np.dtype(column["dtype"]))
        for column in columns)
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

    for index, target in enumerate(document.get("coverage", {}).get("targets", [])):
        prefix = f"$.coverage.targets[{index}]"
        errors.extend(validate_dataframe(target["summary"]["data"], f"{prefix}.summary.data"))
        errors.extend(validate_dataframe(
            target["distribution"]["data"], f"{prefix}.distribution.data"))

    for name, heatmap in document["heatmaps"].items():
        if heatmap is not None:
            errors.extend(validate_ndarray(heatmap["data"], f"$.heatmaps.{name}.data"))
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
    document = summary_json.load(path)
    errors = validate_document(document, schema)
    if errors:
        messages = "\n".join(f"  {error}" for error in errors)
        raise SystemExit(f"{path} does not match {schema_path}:\n{messages}")
    print(f"validated {path}")
    return document


def main(args):
    schema = summary_json.load(args.schema)
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
