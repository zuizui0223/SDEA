"""SDEA v16 original Ulex identity-placebo and site-average-K comparison."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.kinship_identity import REPS,SEED,run
from sdea.ulex_partial import EXPECTED_MD5,retrieve_pinned


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--fetch-pinned',action='store_true')
    ap.add_argument('--source-root',type=Path)
    ap.add_argument('--reference-v15',type=Path,default=ROOT/'results/v15_ulex_partial.json')
    ap.add_argument('--reps',type=int,default=REPS)
    ap.add_argument('--seed',type=int,default=SEED)
    ap.add_argument('--output',type=Path,default=ROOT/'results/v16_ulex_genomic_placebo.json')
    x=ap.parse_args()
    if x.fetch_pinned == (x.source_root is not None):
        ap.error('Choose exactly one of --fetch-pinned or --source-root')
    payload=retrieve_pinned() if x.fetch_pinned else {
        name:(x.source_root/name).read_bytes() for name in EXPECTED_MD5}
    previous=json.loads(x.reference_v15.read_text())
    report=run(payload,expected_v15=previous,reps=x.reps,seed=x.seed)
    x.output.parent.mkdir(parents=True,exist_ok=True)
    x.output.write_text(json.dumps(report,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    print('Original source MD5s verified:',len(report['source']['source_MD5s_verified']))
    print('Original v15 frozen receipt recomputed:',report['original_v15_frozen_receipt_recomputed_exactly'])
    print('Real genomic increment:',report['real_genomic_gain'])
    print('Site-average genomic increment:',report['site_average_K_gain'])
    print('Genomic identity specificity:',report['identity_specific_increment_vs_site_average_K'])
    print('Site-label-preserving genotype reassignment placebo:',report['within_site_genotype_identity_placebo']['equal_site_genomic_gain'])
    print('Decision:',report['decision'])
    if x.fetch_pinned:
        print('SDEA_V16_RESULT_JSON='+json.dumps(report,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False))


if __name__=='__main__':main()
