"""Recompute deterministic SDEA v6 finite-source-budgets synthetic ensemble."""
import argparse
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.finite_budget import benchmark


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--replicates',type=int,default=3000)
    parser.add_argument('--output',type=Path,default=ROOT/'results/v6_finite_budget_mc.json')
    args=parser.parse_args()
    x=benchmark(args.replicates)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(x,sort_keys=True,indent=2,ensure_ascii=False)+'\n')
    print('SDEA v6:',args.output,'cases',len(x['scenarios']),'replicates',args.replicates)
    for name,case in x['scenarios'].items():
        p=case['probabilities']
        print(name,'patch0',round(p['patch0_occupied'],4),'patch1 first colonized',round(p['patch1_first_colonization'],4),
              'patch1 joint ever',round(p['patch1_ever_joint'],4))


if __name__=='__main__':
    main()
