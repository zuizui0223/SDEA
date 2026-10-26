"""Run SDEA v3 exact finite-diploid selfing/outcrossing comparison."""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.mendelian import benchmark


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,default=ROOT/'results/v3_mendelian_exact.json')
    args=ap.parse_args()
    result=benchmark()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+'\n')
    print('Wrote',args.output)
    for a,case in result['matched_founder_structure'].items():
        a0=case['by_founder']['two_heterozygotes']['first_step']['offspring_target_absence']
        a1=case['by_founder']['one_homozygote']['first_step']['offspring_target_absence']
        print('assurance',a,'allele-loss risk two het/homozygote',f'{a0:.9f}',f'{a1:.9f}')


if __name__=='__main__':main()
