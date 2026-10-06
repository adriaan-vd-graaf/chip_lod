"""Hand-written, byte-deterministic output writers and small numeric helpers.

Formats (shared with the R port):
  float -> "%.12e", int -> "%d", missing (None) -> "NA" in CSV / null in JSON.
  CSV: header row, comma separator, no quoting, "\n" line endings, UTF-8 without BOM.
  JSON: fixed key order (insertion order), two-space indent, trailing newline.
Non-finite floats raise ValueError.
"""

import json
import math


def format_value(x):
    if x is None:
        return "NA"
    if isinstance(x, bool):
        raise TypeError("booleans are not allowed in CSV output")
    if isinstance(x, int):
        return "%d" % x
    if isinstance(x, float):
        if not math.isfinite(x):
            raise ValueError(f"refusing to write non-finite value {x!r}")
        return "%.12e" % x
    raise TypeError(f"unsupported value type {type(x).__name__}")


def write_csv(path, header, rows):
    lines = [",".join(header)]
    for row in rows:
        if len(row) != len(header):
            raise ValueError("row length does not match header")
        lines.append(",".join(format_value(v) for v in row))
    data = ("\n".join(lines) + "\n").encode("utf-8")
    with open(path, "wb") as fh:
        fh.write(data)


def _json_scalar(x):
    if x is None:
        return "null"
    if isinstance(x, bool):
        return "true" if x else "false"
    if isinstance(x, str):
        return '"' + x.replace("\\", "\\\\").replace('"', '\\"') + '"'
    return format_value(x)


def _json_lines(obj, indent):
    pad = "  " * indent
    pad_in = "  " * (indent + 1)
    if isinstance(obj, dict):
        if not obj:
            return "{}"
        items = [pad_in + _json_scalar(str(key)) + ": " + _json_lines(val, indent + 1)
                 for key, val in obj.items()]
        return "{\n" + ",\n".join(items) + "\n" + pad + "}"
    if isinstance(obj, list):
        if not obj:
            return "[]"
        items = [pad_in + _json_lines(val, indent + 1) for val in obj]
        return "[\n" + ",\n".join(items) + "\n" + pad + "]"
    return _json_scalar(obj)


def to_json(obj):
    return _json_lines(obj, 0) + "\n"


def write_json(path, obj):
    with open(path, "wb") as fh:
        fh.write(to_json(obj).encode("utf-8"))


def read_params(path):
    """Reading input is allowed to use the json module (only output is hand-written)."""
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def sum_loop(xs):
    s = 0.0
    for x in xs:
        s += x
    return s


def mean_loop(xs):
    if not xs:
        raise ValueError("mean of empty sequence")
    return sum_loop(xs) / len(xs)


def quantile_ecdf(xs, q):
    """Inverse empirical CDF: element at 1-based index max(1, ceil(q N)) of the sorted values."""
    n = len(xs)
    if n == 0:
        raise ValueError("quantile of empty sequence")
    s = sorted(xs)
    idx = math.ceil(q * n)
    if idx < 1:
        idx = 1
    if idx > n:
        idx = n
    return s[idx - 1]
