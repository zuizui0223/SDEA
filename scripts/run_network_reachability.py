"""Run SDEA v7b directed-network contrasts with frozen independent RNG streams."""
import argparse
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.network_reachability import benchmark

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--replicates',type=int,default=4000)
    p.add_argument('--output',type=Path,default=ROOT/'results/v7b_network_reachability_mc.json')
    args=p.parse_args()
    result=benchmark(args.replicates)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+'\n')
    print('Wrote',args.output)
    e=result['exact']
    print('Exact t1 patch1=',e['cycle']['t1_patch1_colonized'],
          't2 patch2 cycle/return=',e['cycle']['t2_patch2_colonized'],
          e['return']['t2_patch2_colonized'])
    for k,v in result['scenarios'].items():
        p=v['probabilities']
        print(k,'patch2_by_t2=',f'{p["t2_patch2_occupied"]:.6f}',
              'patch2_ever=',f'{p["patch2_ever_occupied"]:.6f}',
              'both_sources_alive_at_h=',f'{p["any_patch_occupied_at_h"]:.6f}')
if __name__=='__main__':main()
