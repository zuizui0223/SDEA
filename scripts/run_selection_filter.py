"""Recompute the frozen SDEA v8 model-generated numerical fixture."""
import argparse
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.selection_filter import benchmark


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--replicates',type=int,default=2000)
    parser.add_argument('--horizon',type=int,default=7)
    parser.add_argument('--seed',type=int,default=1182026)
    parser.add_argument('--output',type=Path,default=ROOT/'results/v8_selection_filter_mc.json')
    args=parser.parse_args()
    results=benchmark(replicates=args.replicates,horizon=args.horizon,seed=args.seed)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results,sort_keys=True,indent=2,ensure_ascii=False)+'\n')
    print(f'Wrote {args.output}; {len(results["scenarios"])} scenarios × {args.replicates} histories')
    for name,row in results['scenarios'].items():
        print(name,
              f'gen1_produced={row["probabilities"]["gen1_joint_seed_produced"]:.4f}',
              f'gen1_recruited={row["probabilities"]["gen1_joint_adult_recruited"]:.4f}',
              f'final_joint={row["probabilities"]["final_joint_adult"]:.4f}')

if __name__=='__main__':
    main()
