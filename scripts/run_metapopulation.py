"""Reproduce independently seeded SDEA v5 two-patch Monte Carlo contrasts."""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.metapopulation import benchmark


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--replicates',type=int,default=9000)
    ap.add_argument('--output',type=Path,default=ROOT/'results/v5_metapopulation_mc.json')
    args=ap.parse_args()
    result=benchmark(args.replicates)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+'\n')
    print('Wrote',args.output)
    for name,cell in result['scenarios'].items():
        p=cell['probabilities']
        print(name, 'initial_vacancy_colonized',f'{p["patch0_colonized_initial_vacancy"]:.5f}',
              'joint_receiver',f'{p["patch0_joint_in_individual"]:.5f}',
              'both_occupied',f'{p["both_patches_occupied"]:.5f}')

if __name__=='__main__':
    main()
