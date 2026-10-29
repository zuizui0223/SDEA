"""Recompute SDEA v12 2022 original-head source-pinned, four-site held-out test."""
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.site_transfer import run_original, verified_input


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--fetch-pinned',action='store_true',help='read exact immutable original author Git SHA and previous source receipt')
    p.add_argument('--original-file',type=Path,default=None,help='optional original byte-copy at pinned source SHA (not supplied by SDEA)')
    p.add_argument('--receipt-file',type=Path,default=None,help='optional prior aza3 receipt at pinned SHA')
    p.add_argument('--output',type=Path,default=ROOT/'results/v12_site_transfer_real.json')
    args=p.parse_args()
    if not args.fetch_pinned and (args.original_file is None or args.receipt_file is None):
        p.error('supply BOTH verified local source files or use --fetch-pinned')
    lock=json.loads((ROOT/'data/v12_source_lock.json').read_text())
    src=lock['sources']
    original=verified_input(src['original_author_seeds'],raw=args.original_file.read_bytes() if args.original_file else None,fetch_pinned=args.fetch_pinned)
    receipt=verified_input(src['aza3_independent_receipt'],raw=args.receipt_file.read_bytes() if args.receipt_file else None,fetch_pinned=args.fetch_pinned)
    result=run_original(original,receipt,lock)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+'\n')
    print('SDEA V12 source verified, 615 original heads, 36 complete site blocks')
    print('Decision:',result['decision'])
    for k,metrics in result['comparisons'].items():
        for metric in ('viable_fraction','emerged_insects_per_head'):
            m=metrics[metric]
            print(k,metric,'MSE zero',round(m['site_mean_mse_zero_reference'],7),'site heldout',round(m['site_mean_mse_train_other_sites'],7),'signs',m['site_sign_match_count'])
    if args.fetch_pinned:
        print('SDEA_V12_RESULT_JSON='+json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False))


if __name__=='__main__':
    main()
