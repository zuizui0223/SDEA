"""Run SDEA v4 exact genotype flow and source-provenance experiments."""
import argparse
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.landscape import benchmark

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,default=ROOT/'results/v4_landscape_exact.json')
    args=ap.parse_args()
    result=benchmark()
    args.output.parent.mkdir(exist_ok=True,parents=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+'\n')
    print('Wrote',args.output)
    for route in ('closed','pollen','seed_matched_carriers','seed_equal_nominal_rate'):
        c=result['ancestral_reopening_pulse'][route]
        print(route,'A expression gen1',c['immediate']['first_generation_expression']['ancestral'],
              'A support gen4',c['horizon_4_pulse']['allele_supported_ancestral'])

if __name__ == '__main__':main()
