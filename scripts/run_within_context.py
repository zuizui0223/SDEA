"""Source-verified SDEA v14 strict year x site x type transport challenge."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.external_stage_transfer import (
    ORIGINAL_URL,ORIGINAL_MD5,RECEIPT_URL,RECEIPT_BLOB_SHA,
    fetch_pinned,run_original_payload as run_v13,
)
from sdea.within_context import run_original_payload


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fetch-pinned',action='store_true')
    p.add_argument('--original',type=Path)
    p.add_argument('--aza3-receipt',type=Path)
    p.add_argument('--v13-fixture',type=Path,default=ROOT/'results/v13_ipomopsis_stage_transfer.json')
    p.add_argument('--output',type=Path,default=ROOT/'results/v14_within_context.json')
    a=p.parse_args()
    if a.fetch_pinned and (a.original is not None or a.aza3_receipt is not None):
        p.error('select pinned online original sources OR two supplied original files')
    if a.fetch_pinned:
        original=fetch_pinned(ORIGINAL_URL,expected_md5=ORIGINAL_MD5)
        prior=fetch_pinned(RECEIPT_URL,expected_blob=RECEIPT_BLOB_SHA)
    else:
        if a.original is None or a.aza3_receipt is None:
            p.error('source: provide --fetch-pinned or BOTH local source paths')
        original=a.original.read_bytes()
        prior=a.aza3_receipt.read_bytes()

    fixed=json.loads(a.v13_fixture.read_text())
    independently_recomputed_v13=run_v13(original,prior)
    if independently_recomputed_v13 != fixed:
        raise RuntimeError('PRIOR_V13_ORIGINAL_SOURCE_FIXTURE_RECOMPUTATION_FAILED')
    result=run_original_payload(original,prior,prior_v13=fixed)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    print('SDEA v14: source MD5 and prior same-study receipt verified; 367 paired-stage records')
    print('v13 original source result independently matches frozen v13 fixture')
    for mode,stages in result['within_context_forecasts'].items():
        for stage,analysis in stages.items():
            print(mode,stage,'gain row/equal',analysis['heldout_row_weighted_mse_gain'],
                  analysis['heldout_equal_year_mse_gain'],'years+',analysis['n_years_width_improves'],
                  'decision',analysis['readiness_gate'])
    if a.fetch_pinned:
        print('SDEA_V14_RESULT_JSON='+json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False))

if __name__=='__main__':main()
