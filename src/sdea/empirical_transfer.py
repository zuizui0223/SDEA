"""SDEA v10: externally held-out mother-route *forecast* comparison.

Original USDA mother records are read through v9's SHA-locked source auditor.
Two models share EXACTLY THE SAME fitted Bernoulli marginals; their sole
contrast is whether same-mother route association is retained. It is not a
heritability, pigment, ancestral-regain, or causal seed-removal test.
"""
from __future__ import annotations

import math
from collections import defaultdict

from .empirical_admission import SourceGateError, parse_maternal_joint

CATEGORIES = ('neither', 'apical_only', 'lateral_only', 'both')
TREATMENTS = ('di', 'dni', 'ndi', 'ndni')
YEARS = ('2007', '2008')
ALPHA = 0.5  # Jeffreys-like symmetric categorical smoothing per 4-state cell
EXPECTED_N = {'2007': 38, '2008': 67}


def _counts(cell: dict) -> tuple[int, int, int, int]:
    try:
        c = tuple(cell[k] for k in CATEGORIES)
        n = cell['n']
    except (KeyError, TypeError) as exc:
        raise SourceGateError('missing maternal four-state cell count') from exc
    if any(type(x) is not int or x < 0 for x in c) or type(n) is not int or n <= 0 or sum(c) != n:
        raise SourceGateError('invalid four-state maternal count; positive n and conservation required')
    return c


def _input_table(summary: dict, *, require_original: bool = True) -> dict[str, tuple[int,...]]:
    try:
        data = summary['year_treatment_strata']
    except (KeyError, TypeError) as exc:
        raise SourceGateError('missing original year-treatment table') from exc
    required = {year + '|' + trt for year in YEARS for trt in TREATMENTS}
    if not isinstance(data, dict) or set(data) != required:
        raise SourceGateError('unexpected years/treatments or dropped stratum')
    table = {k: _counts(v) for k,v in data.items()}
    if require_original:
        if summary.get('source_mothers') != 106 or summary.get('admitted_coherent_mothers') != 105:
            raise SourceGateError('original USDA cohort count changed')
        omitted = summary.get('excluded_incoherent')
        if type(omitted) is not list or len(omitted) != 1 or (omitted[0]['year'], omitted[0]['plant_id']) != ('2007', '327'):
            raise SourceGateError('source coherence exclusion changed')
        for year,n in EXPECTED_N.items():
            if sum(sum(table[year+'|'+t]) for t in TREATMENTS) != n:
                raise SourceGateError('held-out year mother total changed')
    return table


def forecast(counts: tuple[int, ...], *, alpha: float = ALPHA) -> dict[str, tuple[float,...]]:
    """Compare product marginals against joint four-cell posterior predictive.

    Both predictors have identical apical/lateral marginals, so the only
    hypothesised missing state is within-mother joint association. Additive
    alpha ensures every category has positive predictive density.
    """
    if type(counts) not in (tuple,list) or len(counts) != 4 or any(type(n) is not int or n<0 for n in counts):
        raise SourceGateError('counts must be four nonnegative integers')
    if type(alpha) not in (int,float) or not math.isfinite(alpha) or alpha <= 0:
        raise SourceGateError('symmetric smoothing must be finite positive')
    n = sum(counts)
    total = n + 4*alpha
    joint = tuple((x+alpha)/total for x in counts)
    apical = joint[1] + joint[3]
    lateral = joint[2] + joint[3]
    independent = ((1-apical)*(1-lateral), apical*(1-lateral),
                   (1-apical)*lateral, apical*lateral)
    if abs(sum(joint)-1)>1e-12 or abs(sum(independent)-1)>1e-12:
        raise AssertionError('forecast must be a valid categorical distribution')
    # These are exactly *shared* marginals, even after the smoothing prior.
    if abs(independent[1]+independent[3]-apical)>1e-12 or abs(independent[2]+independent[3]-lateral)>1e-12:
        raise AssertionError('marginal matching failed')
    return {'independent':independent, 'joint':joint, 'apical_marginal':apical,
            'lateral_marginal':lateral}


def _score(counts: tuple[int,...], probability: tuple[float,...]) -> dict[str,float]:
    """Summed proper four-category log-loss and unhalved categorical Brier."""
    if len(counts)!=4 or len(probability)!=4 or any(x<=0 or not math.isfinite(x) for x in probability):
        raise SourceGateError('invalid categorical prediction')
    norm_sq=sum(p*p for p in probability)
    return {
        'logloss_sum':sum(n*(-math.log(probability[i])) for i,n in enumerate(counts)),
        'brier_sum':sum(n*(1-2*probability[i]+norm_sq) for i,n in enumerate(counts)),
    }


def _score_single(category: int, probability: tuple[float,...]) -> dict[str,float]:
    if type(category) is not int or not 0<=category<4:
        raise SourceGateError('unrecognized four-category outcome')
    return _score(tuple(int(j==category) for j in range(4)),probability)


def heldout_forecast(table: dict[str,tuple[int,...]]) -> dict:
    """Train on all mothers in one year, test only mothers from the other.

    Same named treatment used for matching, no target-year outcomes enter
    training. Reversing direction uses the other, disjoint year cohort.
    """
    folds = {}
    for train, test in ((YEARS[0],YEARS[1]),(YEARS[1],YEARS[0])):
        cells=[]
        sums = {'independent':{'logloss_sum':0.,'brier_sum':0.},
                'joint':{'logloss_sum':0.,'brier_sum':0.}}
        n_total = 0
        for t in TREATMENTS:
            x = table[train+'|'+t]
            y = table[test+'|'+t]
            fits = forecast(x)
            row={'treatment':t, 'training_mothers':sum(x),'heldout_mothers':sum(y),
                 'train_category_counts':dict(zip(CATEGORIES,x)),
                 'heldout_category_counts':dict(zip(CATEGORIES,y)),
                 'same_marginals':{
                    'apical':fits['apical_marginal'],'lateral':fits['lateral_marginal']},
                 'prediction':{},'heldout_scores':{}}
            for model in ('independent','joint'):
                row['prediction'][model]=dict(zip(CATEGORIES,fits[model]))
                metrics=_score(y,fits[model]);row['heldout_scores'][model]=metrics
                for k,v in metrics.items():sums[model][k]+=v
            n_total += sum(y)
            cells.append(row)
        folds[train+'_to_'+test]={
            'train_year':train,'heldout_year':test, 'heldout_mothers':n_total,
            'by_treatment':cells,
            'model_scores':{name:{'logloss_per_mother':z['logloss_sum']/n_total,
                                  'categorical_brier_per_mother':z['brier_sum']/n_total}
                            for name,z in sums.items()},
            'joint_minus_independent':{
                'logloss_per_mother':(sums['joint']['logloss_sum']-sums['independent']['logloss_sum'])/n_total,
                'categorical_brier_per_mother':(sums['joint']['brier_sum']-sums['independent']['brier_sum'])/n_total}
        }
        if n_total != EXPECTED_N[test]:
            raise SourceGateError('held-out year denominator differs from locked mothers')
    nall=sum(f['heldout_mothers'] for f in folds.values())
    pooled={}
    for metric in ('logloss_per_mother','categorical_brier_per_mother'):
        pooled[metric]=sum(f['joint_minus_independent'][metric]*f['heldout_mothers'] for f in folds.values())/nall
    supported=all(f['joint_minus_independent'][key]<0 for f in folds.values()
                  for key in ('logloss_per_mother','categorical_brier_per_mother'))
    return {'folds':folds,'tested_mother_year_units':nall,
            'pooled_joint_minus_independent':pooled,
            'predeclared_transfer_gate': 'JOINT_TRANSFER_GAIN_IN_BOTH_YEARS_AND_BOTH_SCORES' if supported
                                         else 'NO_REPLICATED_JOINT_FORECAST_GAIN',
            'interpretation':'retrospective year-heldout 4-state prediction; no independent taxa/lineages or future descendant data'}


def leave_one_mother_out_within_year(table: dict[str,tuple[int,...]]) -> dict:
    """Exploratory secondary validation. Scores one left-out maternal record.

    Exactly subtract each observed mother-category from the training stratum.
    Never pass a mother's own event into its prediction. This is a within-year
    *same-treatment* generalization, not temporal transport.
    """
    report={}
    for year in YEARS:
        summed={model:{'logloss_sum':0.,'brier_sum':0.} for model in ('independent','joint')}
        n_year=0
        for t in TREATMENTS:
            whole=table[year+'|'+t]
            for cat,freq in enumerate(whole):
                if not freq:continue
                train=list(whole)
                train[cat]-=1
                preds=forecast(tuple(train))
                for model in ('independent','joint'):
                    score=_score_single(cat,preds[model])
                    for k,z in score.items():summed[model][k]+=freq*z
                n_year+=freq
        report[year]={
            'heldout_mothers':n_year,
            'model_scores':{model:{'logloss_per_mother':z['logloss_sum']/n_year,
                                   'categorical_brier_per_mother':z['brier_sum']/n_year}
                            for model,z in summed.items()},
        }
        report[year]['joint_minus_independent']={key:report[year]['model_scores']['joint'][key]-report[year]['model_scores']['independent'][key]
                                             for key in ('logloss_per_mother','categorical_brier_per_mother')}
    return report


def association_by_stratum(table: dict[str,tuple[int,...]]) -> dict:
    result={}
    for year in YEARS:
        for t in TREATMENTS:
            n,a_only,l_only,both=table[year+'|'+t]
            total=n+a_only+l_only+both
            pa=(a_only+both)/total; pl=(l_only+both)/total
            joint=both/total
            cov=joint-pa*pl
            success=1-n/total
            independent=pa+pl-pa*pl
            if abs(success-independent+cov)>1e-12:
                raise AssertionError('observed union vs independence must equal -covariance')
            result[year+'|'+t]={
                'n':total,'observed_any_good_seed':success,
                'apical_success':pa,'lateral_success':pl,
                'same_mother_covariance':cov,
                'independent_marginal_union':independent,
                'observed_minus_independent':success-independent,
            }
    return result


def run_from_summary(summary: dict, *, mother_entry: dict) -> dict:
    table=_input_table(summary)
    transfer=heldout_forecast(table)
    within=leave_one_mother_out_within_year(table)
    association=association_by_stratum(table)
    s2007=association['2007|ndni']['same_mother_covariance'];s2008=association['2008|ndni']['same_mother_covariance']
    return {
        'protocol':'SDEA v10 (2026-10-28) original-USDA held-out joint seed route transfer audit',
        'status':'RETROSPECTIVE EXPLORATORY ORIGINAL-DATA PREDICTION; not independent taxon or ancestral re-evolution validation',
        'source':{'repository':mother_entry['repo'], 'commit':mother_entry['commit'],
                  'path':mother_entry['path'],'verified_git_blob_sha':mother_entry['blob_sha'],
                  'original_doi':mother_entry['doi'],
                  'source_mothers':summary['source_mothers'],
                  'admitted_coherent_mothers':summary['admitted_coherent_mothers'],
                  'omitted_original_mother_year':summary['excluded_incoherent']},
        'design':{'four_categories':list(CATEGORIES),'years':list(YEARS),
                  'matched_treatments':list(TREATMENTS), 'symmetric_cell_pseudocount':ALPHA,
                  'same_marginals_in_both_predictive_models':True,
                  'primary_holdout':'whole source year, all mothers from that year in test; same-treatment match only',
                  'primary_metrics':['logloss_per_mother','categorical_brier_per_mother'],
                  'secondary':'leave-one-mother-out inside year+treatment',
                  'test_unit':'one maternal plant-year, not plant/head/offspring independent replicate',
                  'never_claim':'ancestral trait re-expression, gene inheritance, adaptation or permanent organ rescue'},
        'observed_strata':association,
        'observed_same_treatment_year_covariance_sign_change':{'treatment':'ndni',
           '2007':s2007,'2008':s2008,'opposite_sign':s2007*s2008<0},
        'year_heldout_transfer':transfer,
        'same_year_leave_one_mother_out_secondary':within,
        'biological_decision':'OBSERVED_SAME_MOTHER_ROUTE_DEPENDENCE_BUT_NO_REPLICATED_YEAR_TRANSFER' if transfer['predeclared_transfer_gate']=='NO_REPLICATED_JOINT_FORECAST_GAIN' else 'JOINT_PREDICTIVE_GAIN_WITHIN_THIS_TWO_YEAR_ARCHIVE_ONLY',
        'interpretation_ceiling':'two anatomical good-seed routes on one maternal plant, no measured pollinator co-action, gene toolkit, independently dated ancestral history, re-expression or offspring establishment',
    }


def run_original_payload(raw: bytes, *, mother_entry: dict) -> dict:
    """Raw integrity must be validated BEFORE this call by v9 verify source gate."""
    from .empirical_admission import validate_source_bytes
    validate_source_bytes(mother_entry,raw)
    return run_from_summary(parse_maternal_joint(raw,check_original=True),mother_entry=mother_entry)
