"""Run pinned independent Ulex source partial genotype/trait/fruit test."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.ulex_partial import EXPECTED_MD5,UlexGateError,retrieve_pinned,analyze

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--fetch-pinned',action='store_true',help='allow source downloads after published MD5 audit')
    p.add_argument('--source-root',type=Path,help='directory with original same-MD5 registered Zenodo files')
    p.add_argument('--output',type=Path,default=ROOT/'results/v15_ulex_partial.json')
    a=p.parse_args()
    if a.fetch_pinned == (a.source_root is not None):p.error('choose exactly one of --fetch-pinned / --source-root')
    payloads=retrieve_pinned() if a.fetch_pinned else {name:(a.source_root/name).read_bytes() for name in EXPECTED_MD5}
    report=analyze(payloads)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    print('Ulex original file MD5 checks',len(report['source']['source_file_md5_verified']))
    print('Joined genotypes/tratis/fruits:',report['audit']['matched_genomic_flower_fruit_individuals'])
    print('Trait incremental:',report['trait_increment'])
    print('Genetic residual incremental:',report['genetic_increment_conditional_on_traits'])
    if a.fetch_pinned:
        print('SDEA_V15_RESULT_JSON='+json.dumps(report,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False))
if __name__=='__main__':main()
