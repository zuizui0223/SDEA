"""SDEA v17: source-verified original Sinapis same-plant floral/process/reproduction test.

Original 2022 rainfall experiment. Descriptive out-of-shelter predictions only;
no causal mediation, genotype, ancestral trait, or descendant fitness inference.
No raw source rows exported. Standard library only.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import random
from collections import Counter, defaultdict
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

RECORD_ID=7335086
API_URL='https://zenodo.org/api/records/7335086'
ORIGINAL_MD5={
    'Floral_traits.csv':'24f1f01435fa22faa832c415bbd0b77e',
    'Int_data.csv':'d67c32cc5aa357fe6e6d952a801ea38f',
    'Plantlist.csv':'630eb7d3d7c871f65c2bd62890e1597c',
    'README.md':'68bf11faaa9786fbb39f9e32006b336d',
}
TREATMENTS=('Control','Drought','PrecipRed')
SHELTERS=('1','2','3','4','5')
N_ORIGINAL_PLANTS=53
N_ORIGINAL_VISIT_GROUP_ROWS=427
RIDGE=4.0
PLACEBO_REPS=200
PLACEBO_SEED=20261105
FRUIT_FIELDS=('plantID','plotID','shelterID','motherID','treatment','SWP','meanflower','sumflower','meanheight','meaninflor','meandisplayflower','meanstamen','meanstyle','meancalyx','meanpetalwidth','meanpetallength','meaninflorsize','timetoflower','floweringtime','meannectar','pollen_flower','pods','avgratiopods','avgseed','seedsperplant','seed_mass','AGB_mass','root_mass','ratioAGB_root','pollenstigma','meanleaves')
VISIT_FIELDS=('plantID','visitorgroup','motherID','treatment','shelter','plotID','meanflower','observetime','sumday','NrInteraction','NrInt_day','NrInt_day_flower')
PLANT_FIELDS=('plantID','motherID','treatment','shelter','plotID')
class SinapisGateError(ValueError):
    """Original source or biological grain invalid: do not impute a null."""

def _md5(b:bytes)->str:
    return hashlib.md5(b).hexdigest()

def _safe_url(u:str)->bool:
    x=urlsplit(u)
    return x.scheme=='https' and x.hostname=='zenodo.org' and not x.username and not x.password

def get_sources(*, timeout=65)->dict[str,bytes]:
    """Original registered Zenodo files, MD5 verified before biological parse."""
    try:
        with urlopen(Request(API_URL,headers={'User-Agent':'SDEA-v17-original-source/1.0'}),timeout=timeout) as resp:
            meta=json.load(resp)
    except Exception as e:
        raise SinapisGateError('ORIGINAL_SOURCE_METADATA_UNAVAILABLE') from e
    if str(meta.get('id'))!=str(RECORD_ID):
        raise SinapisGateError('WRONG_ORIGINAL_RECORD')
    fs=meta.get('files',[])
    if isinstance(fs,dict):fs=fs.get('entries',fs)
    if isinstance(fs,dict):fs=[dict(item,key=name) for name,item in fs.items()]
    if not isinstance(fs,list):raise SinapisGateError('ORIGINAL_FILE_INDEX_SCHEMA')
    idx={str(x.get('key',x.get('filename'))):x for x in fs}
    if not set(ORIGINAL_MD5)<=set(idx):raise SinapisGateError('ORIGINAL_FILE_MISSING')
    out={}
    for name,want in ORIGINAL_MD5.items():
        ent=idx[name]
        checksum=str(ent.get('checksum','')).lower().split(':')[-1]
        if checksum and checksum!=want:
            raise SinapisGateError('ORIGINAL_METADATA_MD5_CHANGED:'+name)
        links=ent.get('links') or {}
        url=links.get('content') or links.get('download') or links.get('self')
        if not isinstance(url,str) or not _safe_url(url):
            raise SinapisGateError('ORIGINAL_DOWNLOAD_URL_REJECTED:'+name)
        try:
            with urlopen(Request(url,headers={'User-Agent':'SDEA-v17-original-source/1.0','Accept':'*/*'}),timeout=timeout) as resp:
                raw=resp.read(90_001)
        except Exception as e:
            raise SinapisGateError('ORIGINAL_FILE_UNAVAILABLE:'+name) from e
        if len(raw)>90_000 or _md5(raw)!=want:
            raise SinapisGateError('ORIGINAL_FILE_MD5_OR_SIZE_DRIFT:'+name)
        out[name]=raw
    return out

def validate_bytes(src:dict[str,bytes])->None:
    if set(src)!=set(ORIGINAL_MD5):raise SinapisGateError('ORIGINAL_FILE_SET_CHANGED')
    for f,want in ORIGINAL_MD5.items():
        if _md5(src[f])!=want:raise SinapisGateError('ORIGINAL_FILE_MD5_DRIFT:'+f)

def _table(raw:bytes,cols:tuple[str,...])->list[dict[str,str]]:
    try:
        table=csv.DictReader(io.StringIO(raw.decode('utf-8-sig')),delimiter=';')
        actual=tuple(table.fieldnames or ())
        if actual!=cols:raise SinapisGateError('ORIGINAL_COLUMN_SCHEMA_CHANGED:'+str(actual[:6]))
        rows=list(table)
    except (UnicodeDecodeError,csv.Error) as e:
        raise SinapisGateError('ORIGINAL_CSV_PARSE_FAILURE') from e
    if any(None in x or any(v is None for v in x.values()) for x in rows):
        raise SinapisGateError('ORIGINAL_MALFORMED_ROW')
    return rows

def _num(x:str,*,positive=False,nonnegative=False,integer=False)->float|None:
    if x is None or x.strip() in ('','NA','.','#N/A'):
        return None
    try:v=float(x.strip())
    except ValueError as e:raise SinapisGateError('NONNUMERIC_ORIGINAL_VALUE') from e
    if not math.isfinite(v) or (positive and v<=0) or (nonnegative and v<0) or (integer and not v.is_integer()):
        raise SinapisGateError('INVALID_ORIGINAL_NUMBER')
    return v

def _ident(r:dict)->tuple[str,str,str,str,str]:
    return (r['plantID'].strip(),r['motherID'].strip(),r['treatment'].strip(),(r.get('shelter',r.get('shelterID'))).strip(),r['plotID'].strip())

def source_join(src:dict[str,bytes],*,strict=True)->tuple[list[dict],dict]:
    """ID-exact 3-table join, groups are NOT independent experimental plants."""
    if strict:validate_bytes(src)
    plants=_table(src['Plantlist.csv'],PLANT_FIELDS)
    floral=_table(src['Floral_traits.csv'],FRUIT_FIELDS)
    visits=_table(src['Int_data.csv'],VISIT_FIELDS)
    if strict and (len(plants)!=N_ORIGINAL_PLANTS or len(floral)!=N_ORIGINAL_PLANTS or len(visits)!=N_ORIGINAL_VISIT_GROUP_ROWS):
        raise SinapisGateError('ORIGINAL_SOURCE_ROW_COUNTS_CHANGED')
    if not plants or not floral or not visits:
        raise SinapisGateError('EMPTY_ORIGINAL_COHORT')
    plant_map={}
    for r in plants:
        idinfo=_ident(r);id=idinfo[0]
        if not all(idinfo) or id in plant_map or idinfo[2] not in TREATMENTS or idinfo[3] not in SHELTERS:
            raise SinapisGateError('INVALID_OR_DUPLICATE_ORIGINAL_PLANT_ID')
        plant_map[id]=idinfo
    flower_map={};missing=Counter()
    for r in floral:
        key=_ident(r);id=key[0]
        if id not in plant_map or key!=plant_map[id] or id in flower_map:
            raise SinapisGateError('FLORAL_ID_OR_TREATMENT_MISMATCH')
        flower_map[id]=r
    if set(flower_map)!=set(plant_map):
        raise SinapisGateError('ORIGINAL_FLORAL_PLANT_COVERAGE')
    vi=defaultdict(list)
    for r in visits:
        key=_ident(r);id=key[0]
        if id not in plant_map or key!=plant_map[id]:
            raise SinapisGateError('VISITOR_ID_OR_TREATMENT_MISMATCH')
        if not r['visitorgroup'].strip():raise SinapisGateError('MISSING_VISITOR_GROUP')
        vi[id].append(r)
    if strict and set(vi)!=set(plant_map):
        raise SinapisGateError('ORIGINAL_VISITOR_PLANT_COVERAGE_CHANGED')
    result=[];counts=Counter();raw_groups=Counter();visit_time_conflicts=0
    for id,info in plant_map.items():
        if id not in vi:
            counts['no_direct_visitation_record']+=1;continue
        r=flower_map[id];recs=vi[id]
        numeric={}
        for field,kwargs in {
            'sumflower':{'positive':True},
            'seedsperplant':{'nonnegative':True},
            'meandisplayflower':{'positive':True},
            'meanpetalwidth':{'positive':True},
        }.items():
            val=_num(r[field],**kwargs)
            numeric[field]=val
            if val is None:counts['missing_'+field]+=1
        time_vals=[];total_visits=0;groups=set()
        for v in recs:
            t=_num(v['observetime'],positive=True)
            n=_num(v['NrInteraction'],nonnegative=True,integer=True)
            if t is None or n is None:
                raise SinapisGateError('MISSING_VISITOR_TIME_OR_COUNT')
            time_vals.append(t);total_visits+=int(n)
            if n>0:groups.add(v['visitorgroup'].strip())
            raw_groups[v['visitorgroup'].strip()]+=1
        if max(time_vals)-min(time_vals)>1e-6:
            # The many visitor-group rows share the SAME plant observation effort.
            raise SinapisGateError('INCONSISTENT_WITHIN_PLANT_VISITOR_EXPOSURE')
        if any(q is None for q in numeric.values()):
            counts['ineligible_complete_source_outcomes_and_trait']+=1
            continue
        y=math.log1p(numeric['seedsperplant'])
        logflowers=math.log1p(numeric['sumflower'])
        visits_per_hour=total_visits/(time_vals[0]/60.)
        result.append({'shelter':info[3],'treatment':info[2],
                       'log_seed':y,'log_flower':logflowers,
                       'display':numeric['meandisplayflower'],
                       'petal_width':numeric['meanpetalwidth'],
                       'log_visit_rate':math.log1p(visits_per_hour),
                       'n_visitor_groups_positive':len(groups)})
    n_by_shelter=Counter(r['shelter'] for r in result)
    if (len(result)<35 or set(n_by_shelter)!=set(SHELTERS) or
            min(n_by_shelter.values())<5):
        raise SinapisGateError('INSUFFICIENT_ORIGINAL_SAME_PLANT_STAGE_COVERAGE')
    return result,{'plantlist_source_rows':len(plants),
                   'floral_source_rows':len(floral),
                   'visitor_group_source_rows':len(visits),
                   'plants_with_visit_groups':len(vi),
                   'individuals_with_all_model_fields':len(result),
                   'observed_shelter_counts':dict(sorted(n_by_shelter.items())),
                   'ineligibility_counts':dict(sorted(counts.items())),
                   'original_visitor_groups':dict(sorted(raw_groups.items())),
                   'visitor_observation_time_reused_once_per_plant':True,
                   'shelter_is_holdout_unit_not_individual_visitor_row':True,
                   'experimental_original_maternal_family_key_exposed':True,
                   'plant_genetic_alleles_or_ancestors_measured':False}

def _solve(a:list[list[float]],b:list[float])->list[float]:
    n=len(b);m=[list(row)+[b[i]] for i,row in enumerate(a)]
    if not n or any(len(row)!=n+1 for row in m):raise SinapisGateError('INVALID_RIDGE_MODEL')
    for i in range(n):
        k=max(range(i,n),key=lambda j:abs(m[j][i]))
        if abs(m[k][i])<1e-11:raise SinapisGateError('SINGULAR_TRAINING_DESIGN')
        m[i],m[k]=m[k],m[i]
        piv=m[i][i];m[i]=[z/piv for z in m[i]]
        for j in range(n):
            if j==i:continue
            f=m[j][i];m[j]=[u-f*v for u,v in zip(m[j],m[i])]
    return [m[i][-1] for i in range(n)]

MODEL_FEATURES={
 'M0':('log_flower',),
 'M1':('log_flower','display','petal_width'),
 'M2':('log_flower','display','petal_width','log_visit_rate'),
 'M3':('log_flower','log_visit_rate'),
}

def train_model(train:list[dict],*,model='M0',penalty=RIDGE)->dict:
    if model not in MODEL_FEATURES or len(train)<20 or not math.isfinite(penalty) or penalty<=0:
        raise SinapisGateError('INVALID_TRAINING_DESIGN_OR_HYPERPARAMETER')
    feats=MODEL_FEATURES[model]
    stats={}
    for k in feats:
        vals=[float(r[k]) for r in train]
        mu=sum(vals)/len(vals)
        sd=math.sqrt(sum((a-mu)**2 for a in vals)/len(vals))
        stats[k]=(mu,max(1e-12,sd))
    def x(r):return [1.,float(r['treatment']=='Drought'),float(r['treatment']=='PrecipRed')]+[(r[k]-stats[k][0])/stats[k][1] for k in feats]
    X=[x(r) for r in train]
    y=[r['log_seed'] for r in train]
    p=len(X[0]);A=[[0.]*p for _ in range(p)];B=[0.]*p
    for row,yi in zip(X,y):
        for i in range(p):
            B[i]+=row[i]*yi
            for j in range(p):A[i][j]+=row[i]*row[j]
    for i in range(1,p):A[i][i]+=penalty
    return {'model':model,'coefs':_solve(A,B),'stats':stats,'fit_n':len(train)}

def predict(m:dict,r:dict)->float:
    x=[1.,float(r['treatment']=='Drought'),float(r['treatment']=='PrecipRed')]+[(r[k]-m['stats'][k][0])/m['stats'][k][1] for k in MODEL_FEATURES[m['model']]]
    return sum(a*b for a,b in zip(x,m['coefs']))

def score_shelters(rows:list[dict])->dict:
    if {r['shelter'] for r in rows}!=set(SHELTERS):raise SinapisGateError('MISSING_ORIGINAL_SHELTERS')
    folds=[]
    for hold in SHELTERS:
        train=[r for r in rows if r['shelter']!=hold]
        test=[r for r in rows if r['shelter']==hold]
        if len(test)<5 or len(train)<20 or any(r['shelter']==hold for r in train):raise SinapisGateError('HELDOUT_SHELTER_GRAIN_FAILURE')
        fits={k:train_model(train,model=k) for k in MODEL_FEATURES}
        mse={k:sum((predict(fits[k],r)-r['log_seed'])**2 for r in test)/len(test) for k in MODEL_FEATURES}
        folds.append({'shelter':hold,'heldout_plants':len(test),'train_plants':len(train),
                      'mse_models':mse,
                      'trait_gain_m0_to_m1':mse['M0']-mse['M1'],
                      'process_gain_m1_to_m2':mse['M1']-mse['M2'],
                      'process_gain_m0_to_m3':mse['M0']-mse['M3']})
    return {'folds':folds,'total_heldout_individuals_once_each':sum(x['heldout_plants'] for x in folds),
            'gains':{k:{'equal_shelter':sum(z[k] for z in folds)/len(folds),
                         'individual_weighted':sum(z[k]*z['heldout_plants'] for z in folds)/sum(z['heldout_plants'] for z in folds),
                         'improved_shelters':sum(z[k]>0 for z in folds)}
                     for k in ('trait_gain_m0_to_m1','process_gain_m1_to_m2','process_gain_m0_to_m3')}}

def _quantile(xs:list[float],p:float)->float:
    a=sorted(xs);z=(len(a)-1)*p;i=int(z)
    return a[i]*(1-(z-i))+a[min(i+1,len(a)-1)]*(z-i)

def shuffle_placebo(rows:list[dict],*,reps=PLACEBO_REPS,seed=PLACEBO_SEED)->dict:
    if type(reps)is not int or not 10<=reps<=1000 or type(seed)is not int or seed<=0:raise SinapisGateError('INVALID_PLACEBO_PARAMS')
    orig=score_shelters(rows)
    real=orig['gains']['process_gain_m1_to_m2']['equal_shelter']
    rng=random.Random(seed)
    site_trt=defaultdict(list)
    for i,r in enumerate(rows):site_trt[(r['shelter'],r['treatment'])].append(i)
    null=[]
    for _ in range(reps):
        ss=[dict(r) for r in rows]
        for idxs in site_trt.values():
            vv=[rows[i]['log_visit_rate'] for i in idxs]
            rng.shuffle(vv)
            for i,v in zip(idxs,vv):ss[i]['log_visit_rate']=v
        null.append(score_shelters(ss)['gains']['process_gain_m1_to_m2']['equal_shelter'])
    return {'replicates':reps,'seed':seed,'unit':'shuffle visitation predictor across original individuals within exactly each original shelter × rainfall treatment stratum, never response',
            'observed_equal_shelter_increment':real,
            'placebo_distribution':{'q025':_quantile(null,.025),'median':_quantile(null,.5),'q975':_quantile(null,.975),
                                    'min':min(null),'max':max(null)},
            'fraction_placebo_greater_or_equal_observed_descriptive':sum(v>=real for v in null)/len(null),
            'is_randomized_causal_p_value':False}

def analyze(src:dict[str,bytes],*,placebo_reps=PLACEBO_REPS,seed=PLACEBO_SEED)->dict:
    validate_bytes(src)
    rows,audit=source_join(src)
    outcome=score_shelters(rows)
    placebo=shuffle_placebo(rows,reps=placebo_reps,seed=seed)
    g=outcome['gains']['process_gain_m1_to_m2']
    enough=(g['equal_shelter']>0 and g['individual_weighted']>0 and g['improved_shelters']>=4 and
            g['equal_shelter']>placebo['placebo_distribution']['q975'])
    return {
      'version':'SDEA-v17-2026-11-05-original-Sinapis-plant-visitation-floral-seed-stage',
      'status':'RETROSPECTIVE_SOURCE_AUTHENTICATED_ORIGINAL_STUDY_NOT_CAUSAL_INHERITANCE',
      'source':{'original_study':'Hoefer et al. 2022 Sinapis arvensis field rainfall manipulation',
                'zenodo_record':RECORD_ID,'dataset_doi':'10.5061/dryad.66t1g1k5f',
                'original_files_md5_verified':dict(ORIGINAL_MD5),
                'raw_source_observations_committed_to_SDEA':False},
      'audit':audit,
      'estimand':{'same_plant_join':'Plantlist plantID to Floral_traits plantID and Int_data plantID; strict motherID/treatment/shelter/plot ID equality',
                  'holdout':'whole original rain shelter; five folds, one plant scored once',
                  'primary_endpoint':'log1p(original source estimated seedsperplant) – contemporaneous seasonal output, not future offspring',
                  'M0':'rainfall treatment indicators and log1p(total source seasonal flowers)',
                  'M1':'M0 + floral display size + mean petal width',
                  'M2':'M1 + log1p(sum visitor-group NrInteraction / single plant observation effort in hours)',
                  'M3':'M0 + the same observed visit intensity (secondary)',
                  'model':'fixed ridge least squares, penalty 4, train-only standardization',
                  'model_selection_or_tuning_on_test_outcome':False,
                  'site_level_evaluation_not_visitation_row_pseudoreplication':True,
                  'placebo':'within shelter×treatment original-plant visit-rate relabeling, no outcome relabeling'},
      'heldout_results':outcome,
      'within_stratum_visitation_placebo':placebo,
      'descriptive_process_readiness_gate':'PROCESS_ADDS_PORTABLE_INFORMATION_IN_FIVE_SHELTERS' if enough else 'NO_ROBUST_PROCESS_INFORMATION_BEYOND_FLORAL_AND_FLOWER_SUPPLY',
      'cross_study':'CANNOT_NUMERICALLY_TRANSPORT_CIRSIUM_MESH_BAG_OR_IPOMOPSIS_STAGE_EFFECTS',
      'full_genetic_eco_evolutionary_ancestral_regain_admitted':False,
      'scientific_ceiling':'One original 5-shelter experimental rainfall study, contemporaneous 53-ish plant responses; no molecular genotype, actual pollen origin, causal mediation, experimentally separated access, independent new plants, descendants or ancestral trait reexpression.'
    }
