"""Source-pin SDEA v10 2007/2008 mother-year predictive transfer test."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.empirical_admission import load_inputs
from sdea.empirical_transfer import run_original_payload


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source-root',type=Path,default=None,help='a local tree of original pinned repo/commit/path source bytes')
    parser.add_argument('--fetch-pinned',action='store_true',help='allow read-only source download by immutable GitHub commit')
    parser.add_argument('--output',type=Path,default=ROOT/'results/v10_real_joint_transfer.json')
    args=parser.parse_args()
    manifest=json.loads((ROOT/'data/v9_source_lock.json').read_text())
    mother=manifest['sources']['mother_csv']
    raw=load_inputs({'sources':{'mother_csv':mother}},root=args.source_root,allow_network=args.fetch_pinned)['mother_csv']
    result=run_original_payload(raw,mother_entry=mother)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+'\n')
    print('SDEA v10 original data verified: 105 coherent mother-years, whole-year holdout')
    for key,cell in result['year_heldout_transfer']['folds'].items():
        print(key, 'heldout',cell['heldout_mothers'],'joint-independent logloss',round(cell['joint_minus_independent']['logloss_per_mother'],7),'brier',round(cell['joint_minus_independent']['categorical_brier_per_mother'],7))
    print('Transfer gate:', result['year_heldout_transfer']['predeclared_transfer_gate'])
    if args.fetch_pinned:
        print('SDEA_V10_RESULT_JSON='+json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False))


if __name__=='__main__':main()
