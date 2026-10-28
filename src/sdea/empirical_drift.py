"""SDEA v11: post-hoc source-pinned cross-year route-state drift audit.

This is a deliberately EXPLORATORY extension of an observed negative transfer.
The only outcome is same-mother apical/lateral *good-seed production*, not
ancestral trait reactivation, molecular function, descendant survival or fitness.
No held-out year data enter model selection. A covariance decomposition is an
exact descriptive probability identity, NOT a causal attribution of year drift.
"""
from __future__ import annotations

import math
from .empirical_admission import SourceGateError, parse_maternal_joint, validate_source_bytes
from .empirical_transfer import (_input_table, forecast, _score, CATEGORIES,
                                TREATMENTS, YEARS, EXPECTED_N, ALPHA)

LAMBDAS = (0.0, 0.1, 0.25, 0.5, 0.75, 1.0)
SCORES = ('logloss_per_mother', 'categorical_brier_per_mother')


def _check_lambda(lam: float) -> None:
    if type(lam) not in (int, float) or not math.isfinite(lam) or not 0 <= lam <= 1:
        raise SourceGateError('association shrinkage lambda must be a finite probability')


def shrink_association(counts: tuple[int, ...], lam: float,
                       *, alpha: float = ALPHA) -> tuple[float, ...]:
    """Convexly shrink only the *joint association*, preserving BOTH marginals.

    lambda=0: independence product; lambda=1: four-cell joint training law.
    Both endpoint laws are exactly those in frozen SDEA v10, with same alpha.
    """
    _check_lambda(lam)
    p = forecast(counts, alpha=alpha)
    ind, joint = p['independent'], p['joint']
    out = tuple((1-lam)*ind[i] + lam*joint[i] for i in range(4))
    if any(not (0 < z <= 1) for z in out) or abs(sum(out)-1) > 2e-12:
        raise AssertionError('shrunken prediction not a strictly valid four-class probability')
    for indices in ((1,3), (2,3)):
        if abs(sum(out[i] for i in indices)-sum(ind[i] for i in indices)) > 2e-12:
            raise AssertionError('shrinking association changed a marginal')
    return out


def _cell_distribution(counts: tuple[int, ...]) -> dict:
    if (type(counts) is not tuple or len(counts) != 4 or
        any(type(x) is not int or x < 0 for x in counts) or sum(counts) < 1):
        raise SourceGateError('real descriptive cell needs positive integer mother count')
    n=sum(counts)
    p=tuple(v/n for v in counts)
    ap=p[1]+p[3]; lat=p[2]+p[3]
    product=((1-ap)*(1-lat),ap*(1-lat),(1-ap)*lat,ap*lat)
    covariance=p[3]-ap*lat
    residual=(covariance,-covariance,-covariance,covariance)
    if any(abs(p[i]-product[i]-residual[i])>1e-12 for i in range(4)):
        raise AssertionError('exact marginal/association vector identity violated')
    union=1-p[0]
    lower=max(ap,lat); upper=min(1.,ap+lat)
    if not lower-1e-12 <= union <= upper+1e-12:
        raise AssertionError('Fréchet bounds violated by observed union')
    return {'n':n, 'joint_distribution':p,'product_marginals':product,
            'apical_probability':ap,'lateral_probability':lat,
            'covariance':covariance,'association_residual':residual,
            'observed_any_route':union,'frechet_any_route':[lower,upper]}


def _tv(x, y):
    if len(x)!=4 or len(y)!=4:
        raise SourceGateError('TV needs four equally indexed states')
    return sum(abs(a-b) for a,b in zip(x,y))/2


def observed_drift(table: dict[str, tuple[int,...]]) -> dict:
    """Exact *signed-vector* P2008-P2007 decomposition within each treatment.

    Delta(P)=Delta(product of observed marginals) + Delta(association residual).
    Scalar total-variation distances are NOT additive in general and this is
    not causal partition of environmental mechanisms.
    """
    if set(table) != {y+'|'+t for y in YEARS for t in TREATMENTS}:
        raise SourceGateError('all eight original cells required for drift decomposition')
    out={}
    for t in TREATMENTS:
        first=_cell_distribution(table['2007|'+t]);second=_cell_distribution(table['2008|'+t])
        direct=[b-a for a,b in zip(first['joint_distribution'],second['joint_distribution'])]
        marginal=[b-a for a,b in zip(first['product_marginals'],second['product_marginals'])]
        assoc=[b-a for a,b in zip(first['association_residual'],second['association_residual'])]
        if any(abs(a-b-c)>1e-12 for a,b,c in zip(direct,marginal,assoc)):
            raise AssertionError('year shift decomposition does not close')
        tv_j=_tv(first['joint_distribution'],second['joint_distribution'])
        tv_m=_tv(first['product_marginals'],second['product_marginals'])
        tv_a=sum(abs(x) for x in assoc)/2
        if tv_j > tv_m+tv_a+1e-12 or tv_j+1e-12 < abs(tv_m-tv_a):
            raise AssertionError('triangle bounds on total variation violated')
        out[t]={
          'year_2007':first,'year_2008':second,
          'year_difference_2008_minus_2007':{
            'apical_success':second['apical_probability']-first['apical_probability'],
            'lateral_success':second['lateral_probability']-first['lateral_probability'],
            'covariance':second['covariance']-first['covariance'],
            'full_four_state_vector':direct,
            'product_marginal_vector':marginal,
            'joint_association_vector':assoc,
            'total_variation_joint':tv_j,
            'total_variation_product_marginals':tv_m,
            'total_variation_association_component':tv_a,
            'association_sign_changed':first['covariance']*second['covariance']<0,
          },
          'claim_limit':'exact arithmetic partition of two cross-sectional distributions; TV components overlap, cannot be interpreted as causal shares'
        }
    return out


def _train_cells_from_one_year(table, year):
    if year not in YEARS:
        raise SourceGateError('unexpected training year')
    return {t:table[year+'|'+t] for t in TREATMENTS}


def choose_lambda_train_only(train: dict[str, tuple[int,...]],
                             grid: tuple[float,...] = LAMBDAS) -> dict:
    """Nested training-year leave-one-*mother*-out. Never reads test-year rows.

    A SINGLE lambda shared by all four treatments is chosen with LOO log loss,
    no tie-breaking on Brier or held-out year. Grid declared, not optimized to
    retrospective held-out outcome. A model-selection stress test, NOT a new
    prospectively registered natural phenotype prediction trial.
    """
    if type(train) is not dict or set(train) != set(TREATMENTS):
        raise SourceGateError('training must contain four named treatment cells only')
    for cell in train.values():
        _cell_distribution(cell)
    if type(grid) is not tuple or not grid or len(set(grid))!=len(grid):
        raise SourceGateError('lambda grid must be nonempty tuple with unique values')
    for lam in grid:_check_lambda(lam)
    n=sum(sum(train[t]) for t in TREATMENTS)
    rows=[]
    for lam in grid:
        loss=0.;brier=0.
        for t in TREATMENTS:
            counts=train[t]
            for i,number in enumerate(counts):
                if number == 0:continue
                reduced=list(counts);reduced[i]-=1
                p=shrink_association(tuple(reduced),lam)
                loss -= number*math.log(p[i])
                brier += number*(1-2*p[i]+sum(x*x for x in p))
        rows.append({'lambda':lam,'train_only_loo_logloss':loss/n,
                     'train_only_loo_brier':brier/n})
    # deterministic: log loss first, smallest lambda on exact tie
    chosen=min(rows,key=lambda d:(d['train_only_loo_logloss'],d['lambda']))['lambda']
    return {'training_mothers':n,'grid':rows,'selected_lambda':chosen,
            'selection_rule':'minimum same-year mother-LOO log loss; ties to smaller lambda; no target-year data'}


def _heldout_scores(train: dict[str, tuple[int,...]], test:dict[str,tuple[int,...]],
                    lam:float) -> dict:
    _check_lambda(lam)
    if set(train)!=set(TREATMENTS) or set(test)!=set(TREATMENTS):
        raise SourceGateError('different or missing treatment strata')
    total={'logloss_sum':0.,'brier_sum':0.}
    n=0
    for t in TREATMENTS:
        predictions=shrink_association(train[t],lam)
        s=_score(test[t],predictions)
        for k in total:total[k]+=s[k]
        n+=sum(test[t])
    return {'lambda':lam,'heldout_mothers':n,
            'logloss_per_mother':total['logloss_sum']/n,
            'categorical_brier_per_mother':total['brier_sum']/n}


def _score_delta(candidate, baseline):
    return {k:candidate[k]-baseline[k] for k in SCORES}


def train_only_shrinkage_transfer(table:dict[str,tuple[int,...]]) -> dict:
    """Fit lambda on one year, score held-out other year for two proper metrics."""
    if set(table)!={y+'|'+t for y in YEARS for t in TREATMENTS}:
        raise SourceGateError('all eight cells required')
    out={}
    for train_year,test_year in (('2007','2008'),('2008','2007')):
        tr=_train_cells_from_one_year(table,train_year)
        te=_train_cells_from_one_year(table,test_year)
        chosen=choose_lambda_train_only(tr)
        sensitivity=[_heldout_scores(tr,te,lam) for lam in LAMBDAS]
        selected=_heldout_scores(tr,te,chosen['selected_lambda'])
        neutral=_heldout_scores(tr,te,0.)
        joint=_heldout_scores(tr,te,1.)
        if selected['heldout_mothers']!=EXPECTED_N[test_year]:
            raise SourceGateError('heldout denominator changed')
        out[f'{train_year}_to_{test_year}']={
            'train_year':train_year,'heldout_year':test_year,
            'model_selection':chosen,
            'heldout_grid_posthoc_diagnostic_not_for_selection':sensitivity,
            'heldout_independent':neutral,
            'heldout_full_joint':joint,
            'heldout_selected':selected,
            'selected_minus_independent':_score_delta(selected,neutral),
            'full_joint_minus_independent':_score_delta(joint,neutral),
        }
    improved=all(fold['selected_minus_independent'][k]<0 for fold in out.values() for k in SCORES)
    return {
        'folds':out,
        'heldout_mother_years':sum(v['heldout_selected']['heldout_mothers'] for v in out.values()),
        'transport_gate':('TRAIN_ONLY_SHRINKAGE_GAIN_IN_BOTH_YEARS_BOTH_SCORES' if improved
                         else 'NO_REPLICATED_TRAIN_ONLY_SHRINKAGE_TRANSFER_GAIN'),
        'status':'retrospective follow-up after v10 negative result; posthoc model-class sensitivity, not preregistered discovery',
        'no_heldout_year_data_in_tuning':True,
    }


def run_from_summary(summary:dict, *, mother_entry:dict) -> dict:
    table=_input_table(summary)
    decomposition=observed_drift(table)
    transfer=train_only_shrinkage_transfer(table)
    if transfer['heldout_mother_years']!=105:
        raise SourceGateError('original data denominator mismatch')
    return {
        'protocol':'SDEA v11 source-pinned anatomical-route ecological transport decomposition (2026-10-28)',
        'scope':'POSTHOC exploratory two-year same-original-study comparison; no external taxa/genomic reexpression',
        'source':{'repository':mother_entry['repo'],'commit':mother_entry['commit'],
                  'path':mother_entry['path'],'git_blob_sha':mother_entry['blob_sha'],
                  'original_doi':mother_entry['doi'],
                  'source_mothers':summary['source_mothers'],
                  'admitted_coherent_mothers':summary['admitted_coherent_mothers'],
                  'excluded_incoherent':summary['excluded_incoherent']},
        'year_treatment_drift':decomposition,
        'train_year_only_association_shrinkage':transfer,
        'scientific_decision':('REAL_MARGINAL_AND_ASSOCIATION_DRIFT_NO_REPLICATED_SHRINKAGE_TRANSFER'
          if transfer['transport_gate']=='NO_REPLICATED_TRAIN_ONLY_SHRINKAGE_TRANSFER_GAIN'
          else 'EXPLORATORY_TWO_YEAR_SHRINKAGE_GAIN_NEEDS_EXTERNAL_REPLICATION'),
        'decision_ceiling':'proof of descriptive distribution drift and failed forecast transfer in this original study only; not causal influence, ancestor trait re-expression or genetic toolkit portability',
    }


def run_original_payload(raw:bytes, *, mother_entry:dict) -> dict:
    validate_source_bytes(mother_entry,raw)
    return run_from_summary(parse_maternal_joint(raw,check_original=True),mother_entry=mother_entry)
