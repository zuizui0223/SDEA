"""Pinned v9 empirical admission and observed same-mother route replay."""
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.empirical_admission import load_inputs, run


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source-root',type=Path,default=None,help='locally supplied pinned original files by repo/sha/path; no copy is provided by SDEA')
    parser.add_argument('--fetch-pinned',action='store_true',help='download read-only original data at SHA-locked public GitHub raw URLs')
    parser.add_argument('--output',type=Path,default=ROOT/'results/v9_empirical_admission.json')
    x=parser.parse_args()
    manifest=json.loads((ROOT/'data/v9_source_lock.json').read_text())
    payloads=load_inputs(manifest,root=x.source_root,allow_network=x.fetch_pinned)
    report=run(manifest,payloads)
    x.output.parent.mkdir(parents=True,exist_ok=True)
    x.output.write_text(json.dumps(report,sort_keys=True,indent=2,ensure_ascii=False)+'\n')
    print('Verified',len(report['source_pins']),'pinned source blobs; wrote',x.output)
    print('USDA original mothers',report['original_reanalysis']['usda_maternal_good_seed_routes']['source_mothers'],
          'coherent',report['original_reanalysis']['usda_maternal_good_seed_routes']['admitted_coherent_mothers'])
    print('Full causal trait-regain contract:', report['full_contract_decision'])
    if x.fetch_pinned:
        print('SDEA_V9_RESULT_JSON=' + json.dumps(report,sort_keys=True,separators=(',',':'),ensure_ascii=False))


if __name__=='__main__':main()
