#!/usr/bin/env python3
"""Compare exact-model JSON fixtures across Python versions without byte fragility.

Keys, strings, booleans, sequences and integers must match exactly. Only finite
floating-point results tolerate tiny numerical roundoff (not model changes).
"""

import argparse
import json
import math
from pathlib import Path


def assert_same(a, b, path="$"):
    if isinstance(a, dict) and isinstance(b, dict):
        if a.keys() != b.keys():
            raise AssertionError(f"{path}: object keys changed")
        for key in a:
            assert_same(a[key], b[key], f"{path}.{key}")
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            raise AssertionError(f"{path}: list length changed")
        for i, (x, y) in enumerate(zip(a, b)):
            assert_same(x, y, f"{path}[{i}]")
    elif type(a) is float and type(b) is float:
        if not (math.isfinite(a) and math.isfinite(b)
                and math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-13)):
            raise AssertionError(f"{path}: float mismatch {a} vs {b}")
    elif type(a) is type(b) and a == b:
        return
    else:
        raise AssertionError(f"{path}: value/type mismatch {a!r} vs {b!r}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    args = parser.parse_args()
    expected = json.loads(args.reference.read_text(encoding="utf-8"))
    found = json.loads(args.candidate.read_text(encoding="utf-8"))
    assert_same(expected, found)
    print("SDEA fixture verified: exact structure and numerical agreement ≤ 1e-12 relative")


if __name__ == "__main__":
    main()
