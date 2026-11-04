"""Original Zenodo 7335086 source-verified shelter-heldout floral-process-seeds audit."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.sinapis_process import ORIGINAL_MD5,get_sources,analyze

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fetch-pinned',action='store_true',help='retrieve original four Zenodo MD5-pinned source files')
    p.add_argument('--source-root',type=Path,help='directory containing four original author-MD5 files')
    p.add_argument('--output',type=Path,default=ROOT/'results/v17_sinapis_process.json')
    p.add_argument('--reps',type=int,default=200)
    p.add_argument('--seed',type=int,default=20261105)
    args=p.parse_args()
    if args.fetch_pinned==(args.source_root is not None):p.error('choose exactly --fetch-pinned OR --source-root')
    src=get_sources() if args.fetch_pinned else {name:(args.source_root/name).read_bytes() for name in ORIGINAL_MD5}
    r=analyze(src,placebo_reps=args.reps,seed=args.seed)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    print('Verified',len(r['source']['original_files_md5_verified']),'author-MD5 files; joined',r['audit']['individuals_with_all_model_fields'],'plants')
    print('Whole-shelter heldout MSE improvements:',r['heldout_results']['gains'])
    print('Process-identity placebo:',r['within_stratum_visitation_placebo'])
    print('Descriptive gate:',r['descriptive_process_readiness_gate'])
    if args.fetch_pinned:
        print('SDEA_V17_RESULT_JSON='+json.dumps(r,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False))
if __name__=='__main__':main()
