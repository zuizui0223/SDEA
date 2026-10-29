"""Source-authenticated independent Ipomopsis original-study year-holdout gate."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.external_stage_transfer import ORIGINAL_URL,ORIGINAL_MD5,RECEIPT_URL,RECEIPT_BLOB_SHA,fetch_pinned,run_original_payload,ExternalStageError

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--fetch-pinned',action='store_true')
    p.add_argument('--original',type=Path)
    p.add_argument('--aza3-receipt',type=Path)
    p.add_argument('--output',type=Path,default=ROOT/'results/v13_ipomopsis_stage_transfer.json')
    a=p.parse_args()
    if a.fetch_pinned and (a.original or a.aza3_receipt):
        p.error('--fetch-pinned must not be combined with local source paths')
    if a.fetch_pinned:
        raw=fetch_pinned(ORIGINAL_URL,expected_md5=ORIGINAL_MD5)
        previous=fetch_pinned(RECEIPT_URL,expected_blob=RECEIPT_BLOB_SHA)
    else:
        if not a.original or not a.aza3_receipt:
            p.error('provide --fetch-pinned or BOTH explicitly supplied original input files')
        raw=a.original.read_bytes();previous=a.aza3_receipt.read_bytes()
    report=run_original_payload(raw,previous)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    print('Source MD5 and prior aza3 Git blob verified; records',report['audit']['original_records'],
          'paired stages',report['audit']['complete_paired_stage_records'])
    for endpoint,r in report['outcomes'].items():
        print(endpoint,'pooled MSE gain',r['heldout_weighted_mse_gain'],
              'equal-year gain',r['heldout_equal_year_mse_gain'],
              'improved years',r['n_years_width_improves'],'decision',r['readiness_gate'])
    print('Cirsium v12 transfer:',report['external_v12_transport']['decision'])
    if a.fetch_pinned:
        print('SDEA_V13_RESULT_JSON='+json.dumps(report,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False))
if __name__=='__main__':main()
