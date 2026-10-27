"""Regenerate SDEA v7 exact identity and paired 6-patch synthetic runs."""
import argparse, sys, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.network_topology import estimate

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--replicates',type=int,default=2500)
    p.add_argument('--horizon',type=int,default=7)
    p.add_argument('--output',type=Path,default=ROOT/'results/v7_network_topology_mc.json')
    args=p.parse_args()
    data=estimate(replicates=args.replicates,horizon=args.horizon)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,ensure_ascii=False,indent=2,sort_keys=True)+'\n')
    for name in ('seed_alt_minus_cluster','pollen_alt_minus_cluster','both_alt_minus_cluster'):
        c=data['paired_contrasts'][name]['ever_joint_individual_recruited']
        print(name,round(c['difference_variant_minus_baseline'],6),c['paired_mc_95'])
    print('wrote',args.output)
if __name__=='__main__':main()
