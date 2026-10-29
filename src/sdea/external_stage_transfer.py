"""SDEA v13: original-author Ipomopsis stage-specific held-out-year prediction.

This is an independent ORIGINAL STUDY from the Cirsium v9-v12 studies, but
NOT a test of the Cirsium bagging intervention or the SDEA genetic-regain theory.
Only original author-registered data (MD5), and a pinned prior source-audit
receipt (Git blob SHA) are fetched read-only. No raw observation is redistributed.
All choices are retrospective and exploratory in light of earlier aza3 audits.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from collections import Counter
from urllib.request import Request, urlopen

ORIGINAL_URL = ('https://zenodo.org/records/5083399/files/'
                'Selection_by_pollinators_and_predators.csv?download=1')
ORIGINAL_MD5 = '12ed3dbe6d269c5f808f9c6b230d468c'
RECEIPT_URL = ('https://raw.githubusercontent.com/zuizui0223/aza3/'
               'd54aba2dcac3ecbfd4c4fb6bd5b0b88562f3595d/'
               'data/evidence/ipomopsis_two_stage_trait_selection_independent_original_js_v1.json')
RECEIPT_BLOB_SHA = 'df57a0f13ac30ab2253969917eb5617b1cf64552'
HEADER = ('Site','Year','Type','width','color','nectar','sepal',
          'totpinene','indoleN','seeds','seedsinit','flyescape')
YEARS=tuple(str(x) for x in range(2010,2019))
SITES=('agg','hyb','ten')
TYPES=('A','H','T')
SOURCE_YEAR_N=(6,35,37,27,49,56,40,96,91)
SOURCE_SITE_N=(221,144,72)
SOURCE_TYPE_N=(175,246,16)
SOURCE_MISSING=(0,0,0,18,90,160,329,234,265,32,28,51)
COMPLETE_PAIRED_ROWS=367
RIDGE_PENALTY=1.0
ENDPOINTS=('seed_initiation_log1p','fly_escape_fraction')

class ExternalStageError(ValueError):
    """Data, provenance, or identifiability problem: never record as a zero effect."""


def git_blob_sha(raw:bytes)->str:
    return hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()


def fetch_pinned(url:str, *, expected_md5:str|None=None, expected_blob:str|None=None)->bytes:
    try:
        req=Request(url,headers={'User-Agent':'SDEA-v13-auditable-public-data/1.0',
                                 'Accept':'text/csv,application/json,*/*'})
        with urlopen(req,timeout=75) as response:raw=response.read()
    except Exception as exc:
        raise ExternalStageError('ORIGINAL_SOURCE_UNAVAILABLE: '+str(exc)) from exc
    if expected_md5 and hashlib.md5(raw).hexdigest()!=expected_md5:
        raise ExternalStageError('ORIGINAL_SOURCE_MD5_MISMATCH')
    if expected_blob and git_blob_sha(raw)!=expected_blob:
        raise ExternalStageError('PRIOR_SOURCE_RECEIPT_GIT_BLOB_MISMATCH')
    return raw


def ismissing(x:str|None)->bool:
    return x is None or x.strip() in ('.','','NA')


def read_original(raw:bytes, *, strict:bool=True)->tuple[list[dict],dict]:
    if strict and hashlib.md5(raw).hexdigest()!=ORIGINAL_MD5:
        raise ExternalStageError('ORIGINAL_SOURCE_MD5_MISMATCH')
    try:
        table=csv.DictReader(io.StringIO(raw.decode('utf-8-sig')))
        if tuple(table.fieldnames or ())!=HEADER:
            raise ExternalStageError('UNEXPECTED_SOURCE_HEADER')
        records=list(table)
    except UnicodeDecodeError as exc:
        raise ExternalStageError('NON_UTF8_SOURCE') from exc
    for row in records:
        if any(v is None for v in row.values()):
            raise ExternalStageError('TRUNCATED_OR_MALFORMED_ORIGINAL_SOURCE_ROW')
    if strict and len(records)!=437:
        raise ExternalStageError('ORIGINAL_437_ROWS_CHANGED')
    allmissing={k:sum(ismissing(r[k]) for r in records) for k in HEADER}
    counts_year=Counter(r['Year'].strip() for r in records)
    counts_site=Counter(r['Site'].strip() for r in records)
    counts_type=Counter(r['Type'].strip() for r in records)
    if any(r['Year'].strip() not in YEARS or r['Site'].strip() not in SITES or
           r['Type'].strip() not in TYPES for r in records):
        raise ExternalStageError('UNEXPECTED_ORIGINAL_SITE_YEAR_TYPE')
    audit={'original_records':len(records),'original_columns':len(HEADER),
           'missing_per_source_column':allmissing,
           'year_counts':dict(sorted(counts_year.items())),
           'site_counts':dict(sorted(counts_site.items())),
           'type_counts':dict(sorted(counts_type.items())),
           'source_has_stable_individual_id':False,
           'source_site_labels_not_independent_physical_sites':True}
    if strict:
        if (tuple(counts_year.get(y,0) for y in YEARS)!=SOURCE_YEAR_N or
            tuple(counts_site.get(s,0) for s in SITES)!=SOURCE_SITE_N or
            tuple(counts_type.get(t,0) for t in TYPES)!=SOURCE_TYPE_N or
            tuple(allmissing[k] for k in HEADER)!=SOURCE_MISSING):
            raise ExternalStageError('ORIGINAL_SOURCE_PROFILE_DIFFERS_FROM_PINNED_AZA3_AUDIT')
    source_seed_gt=0
    eligible=[]
    for row_no,r in enumerate(records, start=2):
        vals={}
        for k in ('width','seeds','seedsinit','flyescape'):
            if ismissing(r[k]):vals[k]=None
            else:
                try: f=float(r[k].strip())
                except ValueError as exc:
                    raise ExternalStageError('NONNUMERIC_'+k+'_'+str(row_no)) from exc
                if not math.isfinite(f):raise ExternalStageError('NONFINITE_'+k+'_'+str(row_no))
                if k in ('seeds','seedsinit') and f<0:raise ExternalStageError('NEGATIVE_SEED_COUNT')
                if k=='width' and f<=0:raise ExternalStageError('NONPOSITIVE_COROLLA_WIDTH')
                if k=='flyescape' and not 0<=f<=1:raise ExternalStageError('INVALID_FLYESCAPE_FRACTION')
                vals[k]=f
        if vals['seeds'] is not None and vals['seedsinit'] is not None and vals['seeds']>vals['seedsinit']:
            source_seed_gt+=1
        if all(vals[k] is not None for k in ('width','seedsinit','flyescape')):
            eligible.append({'year':r['Year'].strip(),'site_code':r['Site'].strip(),
                             'physical_site':'agg' if r['Site'].strip()=='agg' else 'upper',
                             'type':r['Type'].strip(),'width':vals['width'],
                             'initiated':math.log1p(vals['seedsinit']),
                             'escape':vals['flyescape']})
    audit['seeds_greater_than_seedsinit_rows']=source_seed_gt
    audit['complete_paired_stage_records']=len(eligible)
    audit['eligible_counts_by_year']=dict(sorted(Counter(r['year'] for r in eligible).items()))
    if strict and (source_seed_gt!=41 or len(eligible)!=COMPLETE_PAIRED_ROWS):
        raise ExternalStageError('ORIGINAL_ELIGIBLE_PROFILE_CHANGED')
    if len(eligible)<18 or any(sum(r['year']==y for r in eligible)<2 for y in YEARS):
        raise ExternalStageError('INSUFFICIENT_YEAR_HOLDOUT_SUPPORT')
    return eligible,audit


def validate_prior_receipt(raw:bytes)->dict:
    if git_blob_sha(raw)!=RECEIPT_BLOB_SHA:
        raise ExternalStageError('PRIOR_SOURCE_RECEIPT_GIT_BLOB_MISMATCH')
    try: receipt=json.loads(raw.decode('utf-8-sig'))
    except (UnicodeDecodeError,ValueError) as exc:
        raise ExternalStageError('INVALID_PRIOR_RECEIPT_JSON') from exc
    a=receipt['audit'];s=receipt['source'];model=receipt['analyses']['width']
    if (s['metadata_author_md5']!=ORIGINAL_MD5 or
        a['n_original_data_records']!=437 or
        a['n_original_fields']!=12 or
        a['seeds_greater_than_seedsinit']!=41 or
        model['summary']['n_source_complete_records']!=COMPLETE_PAIRED_ROWS):
        raise ExternalStageError('PRIOR_RECEIPT_SOURCE_PROFILE_DRIFT')
    return {'repo':'zuizui0223/aza3','commit':'d54aba2dcac3ecbfd4c4fb6bd5b0b88562f3595d',
            'blob_sha':RECEIPT_BLOB_SHA,'original_n':a['n_original_data_records'],
            'complete_two_stage_n':model['summary']['n_source_complete_records'],
            'same_original_study_not_an_independent_replication':True}


def _solve(a:list[list[float]],b:list[float])->list[float]:
    """Pivoted Gauss-Jordan for the small penalized normal equations."""
    n=len(b)
    if n==0 or any(len(row)!=n for row in a):raise ExternalStageError('INVALID_LINEAR_DESIGN')
    z=[list(row)+[float(b[i])] for i,row in enumerate(a)]
    for i in range(n):
        p=max(range(i,n),key=lambda j:abs(z[j][i]))
        if abs(z[p][i])<1e-12:raise ExternalStageError('SINGULAR_TRAINING_DESIGN')
        z[i],z[p]=z[p],z[i]
        d=z[i][i]
        z[i]=[v/d for v in z[i]]
        for j in range(n):
            if j==i:continue
            v=z[j][i]
            z[j]=[q-v*t for q,t in zip(z[j],z[i])]
    return [z[i][-1] for i in range(n)]


def _stats(rows:list[dict],feature:str)->tuple[float,float]:
    if not rows:raise ExternalStageError('EMPTY_TRAINING_FOLD')
    vs=[float(row[feature]) for row in rows]
    mu=sum(vs)/len(vs)
    sd=math.sqrt(sum((v-mu)**2 for v in vs)/len(vs))
    return mu,sd if sd>1e-12 else 1.


def _features(row:dict, *, width:bool, endpoint:str, stats:dict)->list[float]:
    x=[1.,float(row['physical_site']=='upper'),float(row['type']=='H'),float(row['type']=='T')]
    if endpoint=='fly_escape_fraction':
        mu,sd=stats['initiated'];x.append((row['initiated']-mu)/sd)
    if width:
        mu,sd=stats['width'];x.append((row['width']-mu)/sd)
    return x


def fit(train:list[dict],*,width:bool,endpoint:str,penalty:float=RIDGE_PENALTY)->dict:
    """Penalty fixed in advance; statistics estimated solely from training rows."""
    if endpoint not in ENDPOINTS:raise ExternalStageError('UNKNOWN_STAGE_ENDPOINT')
    if not isinstance(penalty,(float,int)) or not math.isfinite(penalty) or penalty<=0:
        raise ExternalStageError('INVALID_POSITIVE_RIDGE_PENALTY')
    if not train:raise ExternalStageError('EMPTY_TRAINING_FOLD')
    stats={'initiated':_stats(train,'initiated'),'width':_stats(train,'width')}
    xx=[_features(row,width=width,endpoint=endpoint,stats=stats) for row in train]
    ys=[row['initiated'] if endpoint=='seed_initiation_log1p' else row['escape'] for row in train]
    p=len(xx[0]); a=[[0.]*p for _ in range(p)];b=[0.]*p
    for x,y in zip(xx,ys):
        for j,v in enumerate(x):
            b[j]+=v*y
            for k,z in enumerate(x):a[j][k]+=v*z
    for j in range(1,p):a[j][j]+=float(penalty)
    beta=_solve(a,b)
    return {'beta':beta,'stats':stats,'endpoint':endpoint,'has_width':width,'training_rows':len(train)}


def predict(model:dict,row:dict)->float:
    x=_features(row,width=model['has_width'],endpoint=model['endpoint'],stats=model['stats'])
    y=sum(v*b for v,b in zip(x,model['beta']))
    if model['endpoint']=='fly_escape_fraction':y=min(1.,max(0.,y))
    return y


def heldout_year_score(rows:list[dict],*,year:str,endpoint:str)->dict:
    if year not in YEARS:raise ExternalStageError('UNKNOWN_HELDOUT_YEAR')
    tr=[r for r in rows if r['year']!=year]
    te=[r for r in rows if r['year']==year]
    if not te:raise ExternalStageError('EMPTY_TEST_YEAR')
    if {r['year'] for r in tr}&{r['year'] for r in te}:raise AssertionError('train/test year leakage')
    base=fit(tr,width=False,endpoint=endpoint)
    full=fit(tr,width=True,endpoint=endpoint)
    ykey='initiated' if endpoint=='seed_initiation_log1p' else 'escape'
    a=sum((predict(base,r)-r[ykey])**2 for r in te)/len(te)
    b=sum((predict(full,r)-r[ykey])**2 for r in te)/len(te)
    return {'year':year,'heldout_rows':len(te),'train_rows':len(tr),
            'baseline_mse':a,'plus_width_mse':b,
            'width_information_gain_mse':a-b,
            'width_training_standardized_coefficient':full['beta'][-1],
            'training_only_standardization':True}


def analyze(rows:list[dict],audit:dict,receipt:dict)->dict:
    outcomes={}
    for endpoint in ENDPOINTS:
        folds=[heldout_year_score(rows,year=year,endpoint=endpoint) for year in YEARS]
        n=sum(x['heldout_rows'] for x in folds)
        if n!=len(rows) or any(x['train_rows']!=len(rows)-x['heldout_rows'] for x in folds):
            raise AssertionError('heldout year accounting')
        equal_year_gain=sum(x['width_information_gain_mse'] for x in folds)/len(folds)
        weighted_gain=sum(x['width_information_gain_mse']*x['heldout_rows'] for x in folds)/n
        pos=sum(x['width_information_gain_mse']>1e-12 for x in folds)
        gate=(weighted_gain>0 and equal_year_gain>0 and pos>=7)
        outcomes[endpoint]={
            'heldout_years':len(YEARS),'heldout_original_records_once_each':n,
            'folds':folds,'heldout_weighted_mse_gain':weighted_gain,
            'heldout_equal_year_mse_gain':equal_year_gain,
            'n_years_width_improves':pos,
            'readiness_gate':'WIDTH_ADDS_ROBUST_CROSS_YEAR_INFORMATION_IN_THIS_SOURCE'
                             if gate else 'NO_ROBUST_YEAR_HELDOUT_WIDTH_GAIN'
        }
    return {
      'protocol':'SDEA-v13-2026-10-29-independently-sourced-stage-specific-year-holdout',
      'scientific_status':'POSTHOC_EXPLORATORY_original-study prediction; not prospective validation or source Cirsium treatment replication',
      'source':{'original_author_study':'Campbell et al. 2022 Ipomopsis',
                'paper_doi':'10.1086/716740','original_data_doi':'10.7280/D1KM49',
                'zenodo_record':5083399,'original_file':'Selection_by_pollinators_and_predators.csv',
                'author_registered_file_md5':ORIGINAL_MD5,'source_replay_verified':True,
                'earlier_aza3_same_study_receipt':receipt},
      'audit':audit,
      'design':{'unit':'original published row; stable individual ID absent',
                'physical_site_groups':{'agg':'agg','hyb':'upper','ten':'upper'},
                'holdout':'entire original Year 2010-2018, 9 folds',
                'paired_stage_eligibility':'width, seedsinit and flyescape all observed',
                'stage1_endpoint':'log1p(seedsinit)',
                'stage2_endpoint':'flyescape (fruit fraction, NOT seed survival fraction)',
                'stage1_baseline':'type + physical site',
                'stage2_baseline':'type + physical site + log1p(seedsinit)',
                'candidate_added_variable':'width, standardized on training years only',
                'ridge_penalty_non_intercept':RIDGE_PENALTY,
                'mse_aggregation':'one prediction per original eligible row then equal-year and row-weighted aggregates',
                'cutoff':'strict positive gain for pooled & equal year, at least seven of nine years for stage-specific robust gate',
                'not_an_original_pre_registration':True},
      'outcomes':outcomes,
      'external_v12_transport':{
          'decision':'NOT_IDENTIFIABLE_NO_COMMON_INTERVENTION_OR_ENDPOINT',
          'v12':'Cirsium arvense 2022 four-arm physical bagging vs Control, viable achenes/flowers and insect emergence/head',
          'v13':'Ipomopsis 2010-2018 measured width, initiated seeds, fruit-fly escape fraction; no original bagging indicator or same denominator',
          'no_pooled_cirsium_ipomopsis_treatment_effect':True},
      'scientific_ceiling':'Independent botanical study but single original dataset; does NOT identify causal width selection, sex-specific pollination, predator identity, gene routes, trait regain or intergenerational establishment'
    }


def run_original_payload(original:bytes,prior:bytes)->dict:
    r=validate_prior_receipt(prior)
    rows,audit=read_original(original,strict=True)
    return analyze(rows,audit,r)
