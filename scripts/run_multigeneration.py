"""Recompute v2 exactly, without sampling or external dependencies."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
from sdea.multigeneration import run_multigeneration


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output",default="results/v2_multigeneration_exact.json")
    args = parser.parse_args()
    result = run_multigeneration()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print(f"Generated {path}")
    b = result["baseline"]["complete"]
    print(f"At horizon 8: expression ever {b['first_expression_ever']:.6f}, three surviving expressions {b['three_consecutive_survived_ever']:.6f}, alive and three {b['alive_and_ever_three']:.6f}")
    print("Alignment:",result["alignment"])


if __name__ == "__main__":
    main()
