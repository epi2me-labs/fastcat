"""Decode typed payloads in fastcat summary.json documents."""

import base64
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


DTYPE_WIDTHS = {
    "<u1": 1,
    "<u2": 2,
    "<u4": 4,
    "<u8": 8,
    "<f4": 4,
    "<f8": 8,
}


def load(path):
    """Load a plain or gzip-compressed summary JSON document."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def decode(data):
    """Decode a self-describing array or dataframe container."""
    kind = data.get("kind")
    if kind == "array":
        return decode_ndarray(data)
    if kind == "dataframe":
        return decode_dataframe(data)
    raise ValueError(f"Unsupported summary.json container kind: {kind!r}")


def decode_ndarray(data):
    """Decode a self-describing array into a NumPy array."""
    if data.get("kind") != "array":
        raise ValueError("Expected an array container")
    dtype = data["dtype"]
    shape = tuple(data["shape"])
    if dtype == "string":
        values = data.get("data")
        if not isinstance(values, list):
            raise ValueError("String-array data must be a JSON array")
        if len(shape) != 1 or len(values) != shape[0]:
            raise ValueError("String-array data does not match its shape")
        return np.asarray(values, dtype=str)
    if dtype not in DTYPE_WIDTHS:
        raise ValueError(f"Unsupported summary.json dtype: {dtype}")
    payload = base64.b64decode(data["data"], validate=True)
    size = int(np.prod(shape, dtype=np.int64))
    expected = size * DTYPE_WIDTHS[dtype]
    if len(payload) != expected:
        raise ValueError(
            f"ndarray payload has {len(payload)} bytes, expected {expected} for {shape} {dtype}")
    return np.frombuffer(payload, dtype=dtype).reshape(shape)


def decode_dataframe(data):
    """Decode a self-describing dataframe container into a pandas DataFrame."""
    if data.get("kind") != "dataframe":
        raise ValueError("Expected a dataframe container")
    shape = data["shape"]
    if len(shape) != 2:
        raise ValueError("Dataframe shape must contain row and column counts")
    row_count, column_count = shape
    columns = data["columns"]
    if len(columns) != column_count:
        raise ValueError("Dataframe column count does not match shape")
    arrays = {}
    for column in columns:
        name = column["name"]
        if column.get("kind") != "array":
            raise ValueError(f"Column {name!r} is not an array container")
        column_shape = column.get("shape")
        if column_shape != [row_count]:
            raise ValueError(f"Column {name!r} shape does not match dataframe row count")
        values = decode_ndarray(column)
        arrays[name] = pd.Series(values, dtype="string") if column["dtype"] == "string" else values
    return pd.DataFrame(arrays)
