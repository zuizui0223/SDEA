#!/usr/bin/env python3
"""Write the exact, deterministic SDEA benchmark fixture."""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sdea import run_benchmarks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/v1_exact_benchmarks.json"))
    args = parser.parse_args()
    result = run_benchmarks()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8")
    print(f"Written {args.output} ({len(result['assurance_tradeoff'])} prespecified contrasts)")
    h = result["history"]["brier"]
    a = result["patch_alignment"]["comparison"]
    print(f"Brier visible={h['visible']:.8f}; history={h['history']:.8f}; full={h['full']:.8f}")
    print(f"Matched-coarse-state joint endpoint: aligned={a['aligned']:.6f}, misaligned={a['misaligned']:.6f}")


if __name__ == "__main__":
    main()
