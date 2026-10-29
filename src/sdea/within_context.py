"""SDEA v14. Strict context-conditioned width forecast challenge.

Source provenance and original row admissibility are owned by SDEA v13.  This
module DOES NOT implement ancestral regain, selection, causal pollination or
plant pedigree inference.  It tests a new retrospective, post-v13 conditional
forecast estimand within year x sampling location x type.

For an entirely unseen held-out Year, the TEST COHORT'S UNLABELED covariates
(width, and in stage two observed initiated seeds) may determine within-context
feature means.  Neither test response may enter those means or any fit.  This is
cohort-conditional / transductive prediction, not a single new-flower forecast.
"""
from __future__ import annotations

from collections import defaultdict
from math import isfinite

from .external_stage_transfer import (
    ENDPOINTS, YEARS, ExternalStageError, RIDGE_PENALTY,
    fit, predict, read_original, validate_prior_receipt,
)

MODES=('physical_location','original_site_code')
HOLDOUT_GATE_MIN_YEARS=7


def _group(row:dict, mode:str)->tuple:
    if mode not in MODES:
        raise ExternalStageError('UNKNOWN_CONTEXT_STRATIFICATION')
    site = row['physical_site'] if mode=='physical_location' else row['site_code']
    return (row['year'],site,row['type'])


def _features_by_group(rows:list[dict],mode:str,*,stage:str,stats:dict)->list[dict]:
    """Generate normalized feature-only deviations for an unseen cohort.

    Note: initiated is the TARGET in stage1 and MUST NOT be accessed there.
    In stage2 initiated is explicitly a prior-stage predictor, not flyescape.
    """
    if stage not in ENDPOINTS:
        raise ExternalStageError('UNKNOWN_STAGE_ENDPOINT')
    if not rows:
        raise ExternalStageError('EMPTY_CONDITIONAL_CONTEXT')
    groups=defaultdict(list)
    for i,r in enumerate(rows):
        groups[_group(r,mode)].append(i)
    wd=stats['width'][1]
    sd=stats['initiated'][1]
    result=[None]*len(rows)
    for key,indices in groups.items():
        mw=sum(rows[i]['width'] for i in indices)/len(indices)
        mi=(sum(rows[i]['initiated'] for i in indices)/len(indices)
            if stage=='fly_escape_fraction' else None)
        for i in indices:
            r=rows[i]
            w=(r['width']-mw)/wd
            a=(r['initiated']-mi)/sd if mi is not None else 0.
            result[i]={'group':key,'n_group':len(indices),'w_centered':w,
                       'prior_centered':a}
    if any(x is None for x in result):
        raise AssertionError('lost feature row')
    for key,indices in groups.items():
        if abs(sum(result[i]['w_centered'] for i in indices))>1e-8:
            raise AssertionError('test covariates not centered')
        if stage=='fly_escape_fraction' and abs(sum(result[i]['prior_centered'] for i in indices))>1e-8:
            raise AssertionError('test prior input not centered')
    return result


def _train_context_slope(train:list[dict],mode:str,stage:str,penalty:float=RIDGE_PENALTY)->dict:
    """Estimate only within year x location x type variation from TRAIN rows.

    For the escape stage, make width orthogonal to previous seed initiation
    *within each training cell*, with ridge-stabilized one-column partialling.
    Escape outcomes only enter the final trained slope, not test preprocessing.
    """
    if stage not in ENDPOINTS:
        raise ExternalStageError('UNKNOWN_STAGE_ENDPOINT')
    if not train or not isfinite(penalty) or penalty<=0:
        raise ExternalStageError('BAD_TRAINING_DATA_OR_PENALTY')
    from .external_stage_transfer import _stats
    stats={'width':_stats(train,'width'),'initiated':_stats(train,'initiated')}
    x=_features_by_group(train,mode,stage=stage,stats=stats)
    ys=defaultdict(list)
    key='initiated' if stage=='seed_initiation_log1p' else 'escape'
    for i,row in enumerate(train):ys[x[i]['group']].append(float(row[key]))
    means={g:sum(v)/len(v) for g,v in ys.items()}
    a=[d['prior_centered'] for d in x]
    w=[d['w_centered'] for d in x]
    y=[float(row[key])-means[d['group']] for row,d in zip(train,x)]
    if stage=='fly_escape_fraction':
        den=sum(z*z for z in a)+penalty
        gamma=sum(z*b for z,b in zip(a,w))/den
        by=sum(z*b for z,b in zip(a,y))/den
        xres=[wi-gamma*ai for wi,ai in zip(w,a)]
        yres=[yi-by*ai for yi,ai in zip(y,a)]
    else:
        gamma=0.;xres=w;yres=y
    beta=sum(xi*yi for xi,yi in zip(xres,yres))/(sum(xi*xi for xi in xres)+penalty)
    if not isfinite(beta) or not isfinite(gamma):
        raise ExternalStageError('NONFINITE_CONTEXT_SLOPE')
    return {'beta':beta,'gamma_width_on_prior':gamma,
            'training_group_count':len(ys),'n_training_rows':len(train),
            'n_train_non_singleton_groups':len([v for v in ys.values() if len(v)>1]),
            'n_training_rows_in_non_singletons':sum(len(v) for v in ys.values() if len(v)>1),
            'stats':stats}


def _feature_only_test_rows(rows:list[dict],stage:str)->list[dict]:
    """Never hand a test response column to feature-group centering."""
    if stage=='seed_initiation_log1p':
        return [{k:r[k] for k in ('year','physical_site','site_code','type','width')} for r in rows]
    if stage=='fly_escape_fraction':
        return [{k:r[k] for k in ('year','physical_site','site_code','type','width','initiated')} for r in rows]
    raise ExternalStageError('UNKNOWN_STAGE_ENDPOINT')


def _prediction(rows:list[dict], train:list[dict], mode:str,stage:str)->tuple[list[float],list[float],dict]:
    """Return baseline and augmented predictions using only train Y and test X."""
    if any(r['year'] in {x['year'] for x in train} for r in rows):
        raise ExternalStageError('YEAR_HOLDOUT_LEAKAGE')
    base=fit(train,width=False,endpoint=stage)
    slope=_train_context_slope(train,mode,stage)
    features=_features_by_group(_feature_only_test_rows(rows,stage),mode,stage=stage,stats=slope['stats'])
    pred0=[predict(base,r) for r in rows]
    pred1=[]
    for b,d in zip(pred0,features):
        correction=slope['beta']*(d['w_centered']-slope['gamma_width_on_prior']*d['prior_centered'])
        v=b+correction
        if stage=='fly_escape_fraction':v=max(0.,min(1.,v))
        pred1.append(v)
    info={k:v for k,v in slope.items() if k!='stats'}
    info.update({'heldout_original_rows':len(rows),'heldout_groups':len(set(d['group'] for d in features)),
                 'heldout_singleton_groups':len(set(d['group'] for d in features if d['n_group']==1)),
                 'heldout_rows_with_non_singleton_context':sum(d['n_group']>1 for d in features),
                 'test_centers_use_unlabeled_covariates_only':True,
                 'heldout_outcomes_never_used_in_prediction':True})
    return pred0,pred1,info


def _mse(rows:list[dict], pred:list[float], stage:str)->float:
    key='initiated' if stage=='seed_initiation_log1p' else 'escape'
    if len(rows)!=len(pred) or not rows:raise ExternalStageError('MSE_ROW_ALIGNMENT_ERROR')
    if not all(isfinite(v) for v in pred):raise ExternalStageError('NONFINITE_PREDICTION')
    return sum((float(r[key])-v)**2 for r,v in zip(rows,pred))/len(rows)


def heldout(rows:list[dict], *, mode:str,stage:str)->dict:
    if stage not in ENDPOINTS or mode not in MODES:
        raise ExternalStageError('INVALID_CONTEXT_FORECAST')
    folds=[]
    for year in YEARS:
        train=[r for r in rows if r['year']!=year]
        test=[r for r in rows if r['year']==year]
        if not train or not test:raise ExternalStageError('YEAR_COHORT_MISSING')
        pred0,pred1,info=_prediction(test,train,mode,stage)
        b=_mse(test,pred0,stage);s=_mse(test,pred1,stage)
        folds.append({'year':year,'n_heldout':len(test),'baseline_mse':b,
                      'plus_within_context_width_mse':s,'gain_mse':b-s,
                      **info})
    n=sum(f['n_heldout'] for f in folds)
    if n!=len(rows):raise AssertionError('whole-year accounting broke')
    weighted=sum(x['gain_mse']*x['n_heldout'] for x in folds)/n
    equal=sum(x['gain_mse'] for x in folds)/len(folds)
    positive=sum(x['gain_mse']>1e-12 for x in folds)
    success=weighted>0 and equal>0 and positive>=HOLDOUT_GATE_MIN_YEARS
    return {'mode':mode,'endpoint':stage,'folds':folds,
            'heldout_years':len(folds),'heldout_original_records_once_each':n,
            'n_years_width_improves':positive,
            'heldout_row_weighted_mse_gain':weighted,
            'heldout_equal_year_mse_gain':equal,
            'readiness_gate':('WITHIN_CONTEXT_WIDTH_INFO_RETAINS_YEAR_TRANSFER_IN_THIS_SOURCE'
                              if success else 'NO_ROBUST_WITHIN_CONTEXT_YEAR_TRANSFER'),
            'criterion':'both mean gains >0 and improved in >=7/9 years; POSTHOC source-internal only'}


def analyze(rows:list[dict],audit:dict,prior:dict,prior_v13:dict|None=None)->dict:
    if len(rows)!=367 or audit.get('original_records')!=437:
        raise ExternalStageError('IPOMOPSIS_DATASET_NOT_SOURCE_PINNED')
    results={mode:{stage:heldout(rows,mode=mode,stage=stage) for stage in ENDPOINTS} for mode in MODES}
    if prior_v13 is not None:
        if prior_v13['audit']['complete_paired_stage_records']!=len(rows):
            raise ExternalStageError('PRIOR_V13_ROWS_DRIFT')
        source_base=prior_v13['outcomes']
        earlier={stage:{'v13_heldout_years_improved':source_base[stage]['n_years_width_improves'],
                          'v13_row_weighted_gain':source_base[stage]['heldout_weighted_mse_gain'],
                          'v13_equal_year_gain':source_base[stage]['heldout_equal_year_mse_gain']}
                 for stage in ENDPOINTS}
    else:
        earlier=None
    return {
        'protocol':'SDEA-v14-2026-10-29-posthoc-within-year-site-type-width-transfer',
        'status':'RETROSPECTIVE_SAME_ORIGINAL_IPOMOPSIS_SOURCE; strict context stress-test after v13 observed',
        'source':{'original_study':'Campbell et al 2022 Ipomopsis',
                  'original_data_doi':'10.7280/D1KM49',
                  'author_csv_md5':'12ed3dbe6d269c5f808f9c6b230d468c',
                  'source_verified':True,'prior_receipt':prior},
        'audit':audit,
        'design':{'unit':'original published row (NOT independent verified plant)',
                  'physical_site_merges_original_hyb_and_ten':'same physical upper field; no false 3-site replication',
                  'modes':list(MODES),
                  'within_context':'one original source Year × physical-field or original site label × Type',
                  'holdout':'entire original Year, 9 folds, no heldout Y for model fitting',
                  'test_feature_centering':'unlabeled full heldout-year cohort width and stage2 prior initiation ONLY (transductive)',
                  'training_response_group_means':'training years only',
                  'stage2_prior_known':'observed initiated seeds used as a prior-stage feature; it is not an intervention',
                  'within_width_slope_penalty':RIDGE_PENALTY,
                  'model_selection':'none, all modes + stages declared; retrospective post-v13 test',
                  'endpoint_scoring':'MSE separately for log1p initiated seeds and fruit-fly escaped FRUIT proportion'},
        'v13_previous_same_source_record':earlier,
        'within_context_forecasts':results,
        'cross_original_study_numerical_treatment_effect':'NOT_IDENTIFIABLE_NO_MATCHED_TREATMENT_OR_DENOMINATOR',
        'claim_ceiling':'unlabeled-cohort conditional within-year response forecast ONLY: no causal pollen effect, evolution, stable individual IDs, gene pathways, ancestry, heritable restoration or independent cross-species replication'
    }


def run_original_payload(original:bytes,prior_receipt:bytes,prior_v13:dict|None=None)->dict:
    prior=validate_prior_receipt(prior_receipt)
    rows,audit=read_original(original,strict=True)
    return analyze(rows,audit,prior,prior_v13=prior_v13)
