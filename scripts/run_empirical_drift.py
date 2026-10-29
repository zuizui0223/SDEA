"""SDEA v11 original-USDA source-pinned cross-year drift and shrinkage diagnosis."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.empirical_admission import load_inputs
from sdea.empirical_drift import run_original_payload


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--source-root',type=Path,default=None,
                    help='local original source bytes at repo/SHA/path, not redistributed')
    ap.add_argument('--fetch-pinned',action='store_true',
                    help='enable read-only downloading original Git blob by immutable SHA')
    ap.add_argument('--output',type=Path,default=ROOT/'results/v11_real_transfer_drift.json')
    args=ap.parse_args()
    lock=json.loads((ROOT/'data/v9_source_lock.json').read_text())
    mother=lock['sources']['mother_csv']
    raw=load_inputs({'sources':{'mother_csv':mother}},root=args.source_root,
                    allow_network=args.fetch_pinned)['mother_csv']
    result=run_original_payload(raw,mother_entry=mother)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+'\n')
    print('Source original',result['source']['source_mothers'],'admitted',result['source']['admitted_coherent_mothers'])
    for key,fold in result['train_year_only_association_shrinkage']['folds'].items():
        print(key,'chosen shrinkage',fold['model_selection']['selected_lambda'],
              'joint-minus-independent heldout:',fold['selected_minus_independent'])
    print('Decision',result['train_year_only_association_shrinkage']['transport_gate'])
    if args.fetch_pinned:
        print('SDEA_V11_RESULT_JSON='+json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False))

if __name__=='__main__':main()
