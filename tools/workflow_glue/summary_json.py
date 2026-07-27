"""Decode typed payloads in fastcat summary.json documents."""

import base64
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


DTYPE_WIDTHS = {
    "uint8": 1,
    "uint16": 2,
    "uint32": 4,
    "uint64": 8,
    "float32": 4,
    "float64": 8,
}


def load(path):
    """Load a plain or gzip-compressed summary JSON document."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def decode_ndarray(data):
    """Decode a schema ndarray into a NumPy array."""
    dtype = data["dtype"]
    if dtype not in DTYPE_WIDTHS:
        raise ValueError(f"Unsupported summary.json dtype: {dtype}")
    payload = base64.b64decode(data["__ndarray__"], validate=True)
    shape = tuple(data["shape"])
    size = int(np.prod(shape, dtype=np.int64))
    expected = size * DTYPE_WIDTHS[dtype]
    if len(payload) != expected:
        raise ValueError(
            f"ndarray payload has {len(payload)} bytes, expected {expected} for {shape} {dtype}")
    return np.frombuffer(payload, dtype=dtype).reshape(shape)


def decode_dataframe(data):
    """Decode a schema dataframe block into a pandas DataFrame."""
    row_count = data["shape"][0]
    arrays = {}
    for column in data["columns"]:
        name = column["name"]
        if column["dtype"] == "string":
            values = column.get("data")
            if not isinstance(values, list) or len(values) != row_count:
                raise ValueError(f"String column {name!r} has an invalid row count")
            arrays[name] = pd.Series(values, dtype="string")
        else:
            values = decode_ndarray({
                "__ndarray__": column["__ndarray__"],
                "dtype": column["dtype"],
                "shape": [row_count],
            })
            arrays[name] = values
    return pd.DataFrame(arrays)
