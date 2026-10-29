"""SDEA v15: source-verified Ulex individual-level genotype/trait/fruit partial gate.

Independent original Castellanos et al. 2023 study. Do not infer same-individual
pollinator state (censuses have NO individual key), evolutionary regain,
causal flower selection, genotype heritability or temporally future fruit set.

The predeclared 6-site LOSO models use a fixed linear ridge model and a fixed
nonnegative genomic-similarity residual smoother; neither is tuned by test-site
outcomes. Raw records never appear in outputs.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from collections import Counter, defaultdict
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

RECORD_ID = 7761289
API_URL = 'https://zenodo.org/api/records/7761289'
DOI = '10.5061/dryad.9zw3r22jp'
EXPECTED_MD5 = {
    'fruits.and.mean.floral.traits.csv': '8ca72b7d9721f6d6279e793a1526c43e',
    'floral.traits.csv': '14216bb46af39bfd0b959c1aa092485b',
    'relatedness.matrix.csv': 'ecaa79680ffa25cb45ebdb73c11fe407',
    'pollinator.census.csv': '78494af547de11c1800af50a74c5642f',
    'README.md': '9a6d62c907a93fe7b67268922f756995',
}
EXPECTED_FRUIT_HEADER = ('site','ind','elev','weight','area','scars','fruits','nofruit')
EXPECTED_TRAIT_HEADER = ('site','ind','flower','weight','area','elev')
EXPECTED_POLLINATOR_HEADER = ('year','locality','flowers','tot.visits','Apis','Bombus','other','visits.flower','observations')
RIDGE = 5.0
KINSHIP_SHRINKAGE = 5.0
MIN_MATCHED = 150

class UlexGateError(ValueError):
    """Missing/changed source, key, or estimand is never a natural null."""


def md5(raw:bytes)->str:
    return hashlib.md5(raw).hexdigest()


def valid_url(url: str)->bool:
    u=urlsplit(url)
    return u.scheme=='https' and u.hostname=='zenodo.org' and not u.username and not u.password


def retrieve_pinned(*, timeout: int=90)->dict[str,bytes]:
    """Fetch canonical Zenodo record and only MD5-registered original files."""
    try:
        with urlopen(Request(API_URL,headers={'User-Agent':'SDEA-Ulex-partial-source-lock/1.0'}),timeout=timeout) as resp:
            meta=json.load(resp)
    except Exception as exc:
        raise UlexGateError('ORIGINAL_SOURCE_METADATA_UNAVAILABLE') from exc
    if str(meta.get('id'))!=str(RECORD_ID):
        raise UlexGateError('WRONG_ZENODO_RECORD')
    files=meta.get('files',[])
    if isinstance(files,dict):files=files.get('entries',files)
    if isinstance(files,dict):files=[dict(value, key=name) for name,value in files.items()]
    if not isinstance(files,list):raise UlexGateError('UNEXPECTED_RECORD_FILES_SCHEMA')
    lookup={str(f.get('key',f.get('filename'))):f for f in files}
    if any(name not in lookup for name in EXPECTED_MD5):raise UlexGateError('REQUIRED_ORIGINAL_FILE_MISSING')
    out={}
    for name,want in EXPECTED_MD5.items():
        f=lookup[name]
        checksum=str(f.get('checksum','')).lower().split(':')[-1]
        if checksum and checksum!=want:raise UlexGateError('PUBLISHED_CHECKSUM_DRIFT:'+name)
        urls=f.get('links') or {}
        url=urls.get('content') or urls.get('download') or urls.get('self')
        if not isinstance(url,str) or not valid_url(url):raise UlexGateError('UNTRUSTED_SOURCE_URL:'+name)
        cap=1_000_000 if name=='relatedness.matrix.csv' else 150_000
        try:
            with urlopen(Request(url,headers={'User-Agent':'SDEA-Ulex-partial-source-lock/1.0','Accept':'text/csv,text/plain,*/*'}),timeout=timeout) as resp:
                raw=resp.read(cap+1)
        except Exception as exc:
            raise UlexGateError('ORIGINAL_SOURCE_FILE_UNAVAILABLE:'+name) from exc
        if len(raw)>cap:raise UlexGateError('SOURCE_FILE_EXCEEDS_SIZE_CAP:'+name)
        if md5(raw)!=want:raise UlexGateError('ORIGINAL_MD5_MISMATCH:'+name)
        out[name]=raw
    return out


def verify_payloads(payloads:dict[str,bytes])->None:
    if set(payloads)!=set(EXPECTED_MD5):raise UlexGateError('MISSING_OR_EXTRA_SOURCES')
    for name,want in EXPECTED_MD5.items():
        if md5(payloads[name])!=want:raise UlexGateError('ORIGINAL_MD5_MISMATCH:'+name)


def _csv(raw:bytes, *, header:tuple[str,...]|None=None)->tuple[tuple[str,...],list[dict[str,str]]]:
    try:
        d=csv.DictReader(io.StringIO(raw.decode('utf-8-sig')))
        cols=tuple(d.fieldnames or ())
        if header and cols!=header:raise UlexGateError('ORIGINAL_CSV_HEADER_CHANGED')
        rows=list(d)
    except (UnicodeDecodeError,csv.Error) as exc:raise UlexGateError('ORIGINAL_CSV_ENCODING_OR_PARSE_FAILURE') from exc
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise UlexGateError('MALFORMED_CSV_ROW')
    return cols,rows


def _number(value: str, *, positive:bool=False, nonnegative:bool=False)->float|None:
    if value is None or value.strip() in ('','NA','.'):return None
    try: v=float(value)
    except ValueError as exc:raise UlexGateError('NONNUMERIC_SOURCE_VALUE') from exc
    if not math.isfinite(v) or positive and v<=0 or nonnegative and v<0:
        raise UlexGateError('INVALID_SOURCE_NUMBER')
    return v


def _key(row:dict)->tuple[str,str]:
    return (row['site'].strip(),row['ind'].strip())


def parse_fruit_and_traits(payloads:dict[str,bytes])->tuple[list[dict],dict]:
    fc,fruit=_csv(payloads['fruits.and.mean.floral.traits.csv'],header=EXPECTED_FRUIT_HEADER)
    tc,traits=_csv(payloads['floral.traits.csv'],header=EXPECTED_TRAIT_HEADER)
    pc,poll=_csv(payloads['pollinator.census.csv'],header=EXPECTED_POLLINATOR_HEADER)
    if 'ind' in pc or 'plant' in pc or 'plant_id' in pc:
        raise UlexGateError('POLLINATOR_INDIVIDUAL_ID_SCHEMA_CHANGED')
    if not fruit or not traits or not poll:raise UlexGateError('EMPTY_ORIGINAL_SOURCE_TABLE')
    if len({_key(r) for r in fruit})!=len(fruit):raise UlexGateError('DUPLICATE_SITE_INDIVIDUAL_FRUIT')
    seen=defaultdict(list)
    for r in traits:
        site,ind=_key(r)
        if not site or not ind:raise UlexGateError('EMPTY_TRAIT_KEY')
        seen[(site,ind)].append(r)
    nall=len(fruit); usable=[]; bad=[]; reasons=Counter()
    for r in fruit:
        site,ind=_key(r)
        if not site or not ind:raise UlexGateError('EMPTY_FRUIT_KEY')
        nums={c:_number(r[c],nonnegative=c in ('scars','fruits','nofruit'),positive=c in ('weight','area'))
              for c in ('elev','weight','area','scars','fruits','nofruit')}
        if any(v is None for v in nums.values()):
            reasons['missing_fruit_or_trait_field']+=1
            continue
        scars,fruits,nofruit=(nums[k] for k in ('scars','fruits','nofruit'))
        if not all(v.is_integer() for v in (scars,fruits,nofruit)) or fruits>scars or not math.isclose(nofruit,scars-fruits,abs_tol=1e-9):
            bad.append((site,ind))
            continue
        if scars==0:
            reasons['no_flower_scars_denominator']+=1;continue
        # No overlapping flower-row data is silently assumed to verify a fruit-only trait.
        if (site,ind) not in seen:
            reasons['no_flower_replicate_key']+=1
        usable.append({'site':site,'ind':ind,'elev':nums['elev'],
                       'area':nums['area'],'weight':nums['weight'],
                       'scars':int(scars),'fruits':int(fruits),'fruit_set':fruits/scars})
    if bad:raise UlexGateError('SOURCE_FRUIT_SCARS_ACCOUNTING_INCONSISTENCY')
    sites=Counter(r['site'] for r in usable)
    if len(sites)!=6 or min(sites.values())<10 or len(usable)<180:
        raise UlexGateError('INSUFFICIENT_OR_UNEXPECTED_POPULATION_SUPPORT')
    return usable,{
        'fruit_source_rows':nall,'individual_fruit_rows_usable':len(usable),
        'eligible_by_site':dict(sorted(sites.items())),
        'trait_flower_replicate_rows':len(traits),
        'fruit_records_with_any_repeated_flower_key':sum(_key(r) in seen for r in fruit),
        'pollinator_census_rows_quarantined_unjoined':len(poll),
        'pollinator_plant_id_in_original_header':False,
        'ineligible_reasons':dict(sorted(reasons.items())),
        'fruit_set_definition':'original fruits / original flower scars, source nofruit=scars-fruits',
        'breeding_and_descendant_fitness_not_measured':True,
    }


def parse_relatedness(payloads:dict[str,bytes])->tuple[dict[str,dict[str,float]],dict]:
    cols,rows=_csv(payloads['relatedness.matrix.csv'])
    if len(cols)!=226 or cols[0]!='ind' or len(rows)!=225 or len(set(cols))!=len(cols):
        raise UlexGateError('EXPECTED_225_INDIVIDUAL_RELATEDNESS_MATRIX_CHANGED')
    labels=list(cols[1:]); matrix={}
    for r in rows:
        ident=r['ind'].strip()
        if not ident or ident in matrix:raise UlexGateError('DUPLICATE_GENOTYPE_ID')
        values={}
        for col in labels:
            val=_number(r[col])
            if val is not None:values[col]=val
        matrix[ident]=values
    if set(matrix)!=set(labels):raise UlexGateError('RELATEDNESS_LABELS_DISAGREE_WITH_ROW_IDENTIFIERS')
    pair_count=0;max_asym=0.
    for i in labels:
        for j in labels:
            a=matrix[i].get(j);b=matrix[j].get(i)
            if a is None or b is None:continue
            max_asym=max(max_asym,abs(a-b));pair_count+=1
    if pair_count<.9*225*225 or max_asym>1e-5:
        raise UlexGateError('INSUFFICIENT_OR_ASYMMETRIC_RELATEDNESS_MATRIX')
    return matrix,{'genotyped_individuals':225,'relatedness_observed_cells':pair_count,
                   'relatedness_max_asymmetry':max_asym,
                   'negative_relatedness_values_are_valid':True}


def _solve(a:list[list[float]],b:list[float])->list[float]:
    n=len(b)
    aug=[list(a[i])+[b[i]] for i in range(n)]
    for col in range(n):
        pivot=max(range(col,n),key=lambda k:abs(aug[k][col]))
        if abs(aug[pivot][col])<1e-12:raise UlexGateError('SINGULAR_FIT')
        aug[col],aug[pivot]=aug[pivot],aug[col]
        div=aug[col][col]
        aug[col]=[v/div for v in aug[col]]
        for row in range(n):
            if row==col:continue
            factor=aug[row][col]
            aug[row]=[v-factor*w for v,w in zip(aug[row],aug[col])]
    return [a[-1] for a in aug]


def _stats(train:list[dict],field:str)->tuple[float,float]:
    a=[r[field] for r in train]
    mu=sum(a)/len(a)
    sd=math.sqrt(sum((x-mu)**2 for x in a)/len(a))
    return mu,max(sd,1e-10)


def _vector(r:dict, *, stats:dict, traits:bool)->list[float]:
    elev=(r['elev']-stats['elev'][0])/stats['elev'][1]
    x=[1.,elev]
    if traits:
        a=(r['area']-stats['area'][0])/stats['area'][1]
        w=(r['weight']-stats['weight'][0])/stats['weight'][1]
        x.extend([a,w,a*a,w*w])
    return x


def fit(train:list[dict],*,traits:bool,penalty:float=RIDGE)->dict:
    if len(train)<30 or not math.isfinite(penalty) or penalty<=0:
        raise UlexGateError('INVALID_TRAINING_SUPPORT_OR_PENALTY')
    stats={k:_stats(train,k) for k in ('elev','area','weight')}
    X=[_vector(r,stats=stats,traits=traits) for r in train]
    y=[r['fruit_set'] for r in train]
    p=len(X[0]);A=[[0.]*p for _ in range(p)];b=[0.]*p
    for row,yi in zip(X,y):
        for i in range(p):
            b[i]+=row[i]*yi
            for j in range(p):A[i][j]+=row[i]*row[j]
    for i in range(1,p):A[i][i]+=penalty
    return {'beta':_solve(A,b),'stats':stats,'traits':traits}


def predict(fitobj:dict,r:dict)->float:
    return max(0.,min(1.,sum(a*b for a,b in zip(fitobj['beta'],_vector(r,stats=fitobj['stats'],traits=fitobj['traits'])))))


def related_residual_prediction(test:dict,train:list[dict],fitobj:dict,K:dict[str,dict[str,float]],*,shrinkage:float=KINSHIP_SHRINKAGE)->float:
    if not math.isfinite(shrinkage) or shrinkage<=0:raise UlexGateError('INVALID_GENETIC_SHRINKAGE')
    base=predict(fitobj,test)
    row=K[test['ind']]
    weights=[max(0.,row.get(r['ind'],0.)) for r in train]
    # This is a fixed cross-site K-neighbour residual smoother, not heritability.
    adjustment=sum(w*(r['fruit_set']-predict(fitobj,r)) for w,r in zip(weights,train))/(shrinkage+sum(weights))
    return max(0.,min(1.,base+adjustment))


def analyze(payloads:dict[str,bytes])->dict:
    verify_payloads(payloads)
    fruit,source=parse_fruit_and_traits(payloads)
    K,genetics=parse_relatedness(payloads)
    joined=[r for r in fruit if r['ind'] in K]
    joined_keys=Counter(r['site'] for r in joined)
    if len(joined)<MIN_MATCHED or len(joined_keys)!=6 or min(joined_keys.values())<8:
        raise UlexGateError('INSUFFICIENT_SAME_PLANT_GENOMIC_TRAIT_FRUIT_ALIGNMENT')
    # An individual ID reused in separate sites would be a non-unique genotype join.
    idsites=defaultdict(set)
    for r in fruit:idsites[r['ind']].add(r['site'])
    if any(len(v)>1 for v in idsites.values()):
        raise UlexGateError('REUSED_GENOTYPE_KEY_ACROSS_POPULATIONS')
    sites=sorted(joined_keys)
    per_site=[]
    for site in sites:
        tr=[r for r in joined if r['site']!=site]
        te=[r for r in joined if r['site']==site]
        f0=fit(tr,traits=False)
        f1=fit(tr,traits=True)
        residual=[[0.,0.,0.] for _ in te]
        for i,r in enumerate(te):
            y=r['fruit_set']
            residual[i]=[(predict(f0,r)-y)**2,
                         (predict(f1,r)-y)**2,
                         (related_residual_prediction(r,tr,f1,K)-y)**2]
        mse=[sum(z[i] for z in residual)/len(residual) for i in range(3)]
        per_site.append({'site':site,'heldout_individuals':len(te),'train_individuals':len(tr),
                         'baseline_elevation_mse':mse[0],
                         'plus_floral_area_weight_quadratic_mse':mse[1],
                         'plus_relatedness_residual_mse':mse[2],
                         'trait_mse_gain':mse[0]-mse[1],
                         'relatedness_mse_gain_conditional_on_traits':mse[1]-mse[2]})
    n=sum(d['heldout_individuals'] for d in per_site)
    def aggregate(before:str,after:str)->dict:
        gain=[x[before]-x[after] for x in per_site]
        return {'equal_site_mean_mse_gain':sum(gain)/len(gain),
                'individual_weighted_mse_gain':sum(g*x['heldout_individuals'] for g,x in zip(gain,per_site))/n,
                'sites_with_positive_gain':sum(g>1e-12 for g in gain),
                'site_count':len(per_site),
                'descriptive_readiness_gate':'SIX_SITE_ALL_POSITIVE' if all(g>0 for g in gain) else 'NO_UNIVERSAL_SITE_TRANSFER'}
    return {'protocol':'SDEA-v15-2026-10-29-first-source-admitted-Ulex-partial-G-T-F',
      'status':'RETROSPECTIVE_INDEPENDENT_PUBLIC_STUDY_ORIGINAL_DATA_REANALYSIS_NOT_CAUSAL_OR_LONGITUDINAL',
      'source':{'study':'Castellanos et al. 2023, Ulex parviflorus',
                'paper_doi':'10.1098/rspb.2023.0141','original_data_doi':DOI,
                'zenodo_id':RECORD_ID,'source_file_md5_verified':dict(EXPECTED_MD5),
                'original_raw_observations_in_SDEA':False},
      'audit':{**source,**genetics,'matched_genomic_flower_fruit_individuals':n,
               'matched_by_site':dict(sorted(joined_keys.items())),
               'pollinator_census_individual_link':False,
               'future_genetic_trait_reacquisition_endpoint':False},
      'fixed_predictive_contract':{'validation':'leave_entire_one_original_site_out, six sites, one individual scored once',
             'outcome':'fruit_set=fruits/scars (binomial aggregate proportion, only contemporaneous twig season)',
             'M0':'train-only elevation ridge, penalty 5',
             'M1':'M0 + area + weight + area^2 + weight^2 from train-standardized floral traits',
             'M2':'M1 + fixed nonnegative genomically-related train-residual smoother, shrink 5',
             'prior_individual_genotype_needed_at_test':True,
             'relatedness_does_not_estimate_heritability':True,
             'source_census_pollinator_visits_not_individual_linked':True,
             'mse_endpoint':'unweighted individual fruit-set fraction squared error; report equal-site and individual-weighted summary',
             'all_model_hyperparameters_predeclared':True},
      'heldout_sites':per_site,
      'trait_increment':aggregate('baseline_elevation_mse','plus_floral_area_weight_quadratic_mse'),
      'genetic_increment_conditional_on_traits':aggregate('plus_floral_area_weight_quadratic_mse','plus_relatedness_residual_mse'),
      'full_ANTO_EGWE_izu_model_admitted':False,
      'scientific_ceiling':'Only conditional same-plant G–T–F association, no pollinator plant-ID, causal trait fitness, field multigeneration fate, pedigree direction, ancestral pathway regain or independent source reciprocity.'}
