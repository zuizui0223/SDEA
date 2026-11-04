"""SDEA v16: source-pinned individual-vs-site genomic identity placebo.

A *post hoc* negative-control challenge to the tiny Ulex v15 genomic
predictive increment. Permute the genotype-to-phenotype allocation *within*
each original source population, preserving each population's genomic-label
inventory and the entire observed relatedness matrix. Do not permute outcomes,
source site, plant genotyped status, flower traits or heldout-site fold splits.

The randomization ensemble is descriptive, NOT a causal genotype-effect p-value.
There are no new biological observations, no independently verified offspring,
no same-plant pollinator join and no molecular ancestral pathway measurements.
"""
from __future__ import annotations

import json
import math
import random
from collections import Counter, defaultdict

from . import ulex_partial as u

REPS = 1000
SEED = 20261029
EXPECTED_JOINED = 202
EXPECTED_SITES = ('AresM', 'Cheste', 'Chiva', 'Montserrat', 'Simat', 'Sot.Chera')


class IdentityGateError(u.UlexGateError):
    """Fail closed rather than implying no genetic influence from bad source."""


def _mse(pred: list[float], target: list[float]) -> float:
    if len(pred) != len(target) or not pred:
        raise IdentityGateError('INVALID_MSE_INPUT')
    return sum((x-y)**2 for x, y in zip(pred, target))/len(pred)


def _quantile(values: list[float], q: float) -> float:
    if not values or not 0 <= q <= 1 or not all(math.isfinite(x) for x in values):
        raise IdentityGateError('INVALID_QUANTILE')
    v = sorted(values)
    x = (len(v)-1)*q
    i = int(x)
    return v[i]*(1-(x-i))+v[min(i+1, len(v)-1)]*(x-i)


def _stats(values: list[float]) -> dict:
    return {'min':min(values), 'q025':_quantile(values,.025),
            'median':_quantile(values,.5), 'q975':_quantile(values,.975),
            'max':max(values)}


def _weighted_fold_mean(per_site: list[float], folds: list[dict]) -> float:
    total = sum(len(f['test']) for f in folds)
    return sum(v*len(f['test']) for v,f in zip(per_site,folds))/total


def _fold_data(rows: list[dict], K: dict) -> tuple[list[dict], list[list[float]], dict]:
    """One true genotype ID per source individual; nothing sorted by response."""
    original = sorted([dict(r) for r in rows if r['ind'] in K],
                      key=lambda r:(r['site'], r['ind']))
    if len(original)!=EXPECTED_JOINED:
        raise IdentityGateError('EXPECTED_202_GENOTYPED_INDIVIDUALS')
    if len({r['ind'] for r in original}) != len(original):
        raise IdentityGateError('REUSED_INDIVIDUAL_ID_IN_SOURCE_SITES')
    sites = sorted({r['site'] for r in original})
    if tuple(sites)!=EXPECTED_SITES:
        raise IdentityGateError('EXPECTED_SIX_SOURCE_POPULATIONS')
    grouped = {s:[i for i,r in enumerate(original) if r['site']==s] for s in sites}
    if min(len(idx) for idx in grouped.values())<8:
        raise IdentityGateError('SOURCE_POPULATION_UNDER_COHORT_GATE')
    # Precompute all exact positive-relatedness weights used by the original M2.
    W = [[max(0.,K[r['ind']].get(t['ind'],0.)) for t in original] for r in original]
    folds = []
    for site in sites:
        test = grouped[site]
        train = [i for s in sites if s != site for i in grouped[s]]
        if not test or len(train)+len(test)!=len(original):
            raise IdentityGateError('INVALID_SITE_HELDOUT_SPLIT')
        tr = [original[j] for j in train]
        fit = u.fit(tr,traits=True)
        baseline = [u.predict(fit, original[i]) for i in test]
        truth = [original[i]['fruit_set'] for i in test]
        residual = [original[i]['fruit_set']-u.predict(fit,original[i]) for i in train]
        folds.append({'site':site, 'test':test, 'train':train,
                      'baseline':baseline, 'truth':truth,
                      'residual':residual,
                      'trait_only_mse':_mse(baseline,truth)})
    return folds,W,grouped


def _predict_mapping_fold(fold:dict, W:list[list[float]],
                          mapping:list[int], *,shrinkage:float=u.KINSHIP_SHRINKAGE)->list[float]:
    """Test predictions can use genotypes, traits, and TRAINING outcomes only."""
    if not math.isfinite(shrinkage) or shrinkage<=0:
        raise IdentityGateError('INVALID_SHRINKAGE')
    n=len(mapping)
    if len(W)!=n or sorted(mapping)!=list(range(n)):
        raise IdentityGateError('GENOMIC_LABEL_BIJECTION_REQUIRED')
    mapped_tr=[mapping[j] for j in fold['train']]
    out=[]
    for i,pred_trait in zip(fold['test'],fold['baseline']):
        genome_row=W[mapping[i]]
        numerator=0.; mass=0.
        for j,training_residual in zip(mapped_tr,fold['residual']):
            weight=genome_row[j]
            numerator+=weight*training_residual
            mass+=weight
        out.append(max(0.,min(1.,pred_trait+numerator/(shrinkage+mass))))
    return out


def _score_mapping(folds: list[dict], W: list[list[float]],
                   mapping: list[int], *,shrinkage:float=u.KINSHIP_SHRINKAGE)->list[float]:
    """Calculate six site MSEs under the supplied genotype label assignment."""
    return [_mse(_predict_mapping_fold(f,W,mapping,shrinkage=shrinkage),f['truth'])
            for f in folds]


def _site_average_K(folds: list[dict],W: list[list[float]],*,shrinkage:float=u.KINSHIP_SHRINKAGE)->list[float]:
    """Transductive *site averaged* K: no within-test-site genomic identity."""
    if shrinkage<=0:raise IdentityGateError('INVALID_SHRINKAGE')
    scores=[]
    for fold in folds:
        test,train=fold['test'],fold['train']
        wavg=[sum(W[i][j] for i in test)/len(test) for j in train]
        mass=sum(wavg)
        correction=sum(w*r for w,r in zip(wavg,fold['residual']))/(shrinkage+mass)
        yhat=[max(0.,min(1.,p+correction)) for p in fold['baseline']]
        scores.append(_mse(yhat,fold['truth']))
    return scores


def _within_site_relabeling(grouped:dict[str,list[int]],rng:random.Random,n:int)->list[int]:
    """Permute SNP label assignments only within the six real populations."""
    mapping=list(range(n))
    for site in sorted(grouped):
        indices=grouped[site]
        assigned=indices.copy()
        rng.shuffle(assigned)
        for i,label in zip(indices,assigned):mapping[i]=label
    if sorted(mapping)!=list(range(n)):
        raise IdentityGateError('NONBIJECTIVE_GENOMIC_IDENTITY_PLACEBO')
    return mapping


def analyze_rows(rows:list[dict],K:dict,*,reps:int=REPS,seed:int=SEED,
                 original_v15:dict|None=None)->dict:
    if type(reps) is not int or not 10<=reps<=5000 or type(seed) is not int or seed<1:
        raise IdentityGateError('INVALID_FIXED_GENOTYPE_PLACEBO_BUDGET')
    folds,W,grouped=_fold_data(rows,K)
    baseline=[f['trait_only_mse'] for f in folds]
    real_mse=_score_mapping(folds,W,list(range(len(W))))
    site_coarse_mse=_site_average_K(folds,W)
    # Source v15 was implemented separately and MUST reproduce before a new
    # result is interpreted. Omit only in isolated unit tests with synthetic K.
    if original_v15 is not None:
        old={f['site']:f for f in original_v15['heldout_sites']}
        if set(old)!=set(grouped):raise IdentityGateError('V15_FOLD_SITE_IDENTITY_CHANGED')
        for f,orig_m2 in zip(folds,real_mse):
            v=old[f['site']]
            if abs(f['trait_only_mse']-v['plus_floral_area_weight_quadratic_mse'])>1e-11 or abs(orig_m2-v['plus_relatedness_residual_mse'])>1e-11:
                raise IdentityGateError('ORIGINAL_V15_M1_M2_REPRODUCTION_FAILED')
    per_site=[]
    for fold,trait,original,coarse in zip(folds,baseline,real_mse,site_coarse_mse):
        per_site.append({'site':fold['site'],'heldout_original_genotyped_individuals':len(fold['test']),
            'trait_only_mse':trait,'original_genomic_mse':original,
            'site_average_genomic_mse':coarse,
            'original_genomic_gain_over_traits':trait-original,
            'site_average_genomic_gain_over_traits':trait-coarse,
            'individual_genomic_identity_gain_over_site_averaged_K':coarse-original})
    real_gen_gain=[a-b for a,b in zip(baseline,real_mse)]
    real_specific_gain=[a-b for a,b in zip(site_coarse_mse,real_mse)]
    rng=random.Random(seed)
    null_equal=[];null_weight=[]
    per_site_null=[[] for _ in folds]
    for _ in range(reps):
        mapping=_within_site_relabeling(grouped,rng,len(W))
        perm_mse=_score_mapping(folds,W,mapping)
        gains=[a-b for a,b in zip(baseline,perm_mse)]
        null_equal.append(sum(gains)/len(gains))
        null_weight.append(_weighted_fold_mean(gains,folds))
        for arr,v in zip(per_site_null,gains):arr.append(v)
    observed_equal=sum(real_gen_gain)/len(real_gen_gain)
    observed_weight=_weighted_fold_mean(real_gen_gain,folds)
    avg_site_equal=sum(x['site_average_genomic_gain_over_traits'] for x in per_site)/len(per_site)
    avg_site_weight=_weighted_fold_mean([x['site_average_genomic_gain_over_traits'] for x in per_site],folds)
    specific_equal=sum(real_specific_gain)/len(real_specific_gain)
    specific_weight=_weighted_fold_mean(real_specific_gain,folds)
    n_real_better=sum(x>1e-12 for x in real_specific_gain)
    null_gte=sum(v>=observed_equal for v in null_equal)
    # This is a *descriptive* reference, not a valid exchangeability-based
    # causal p-value or the significance of heritable floral variation.
    readiness=(observed_equal>_quantile(null_equal,.975) and
               specific_equal>0 and specific_weight>0 and n_real_better>=5)
    return {
      'protocol':'SDEA-v16-2026-10-29-posthoc-Ulex-site-preserving-genomic-label-placebo',
      'status':'SOURCE_PINNED_POSTHOC_ORIGINAL_STUDY_IDENTITY_SPECIFICITY_STRESS_TEST',
      'source_n':len(rows),'matched_genotypes':sum(len(f['test']) for f in folds),
      'n_original_study_sites':len(folds),
      'null_design':{
        'permutation_unit':'assign original matched genomic ID to a different individual within the exact same original plant site; one-to-one; no response or trait permutations',
        'replicates':reps,'random_seed':seed,'preserves_site_label_counts':True,
        'preserves_genomic_relatedness_matrix':True,'preserves_holdout_site_identity':True,
        'test_fruit_outcomes_never_used_as_predictor':True,
        'reference_not_exact_genotype_experiment_or_p_value':True,
        'exchangeability_of_identities_within_sites_not_established':True,
        'not_heritable_effect_or_natural_selection_test':True},
      'site_average_comparator':{'mean_K_weights_across_unlabeled_test_site_individuals':True,
        'identity_specific_G_not_used_for_test_personalization':True,
        'requires_test_site_preobserved_genotypes':True},
      'site_folds':per_site,
      'real_genomic_gain':{'equal_site':observed_equal,'individual_weighted':observed_weight,
        'sites_positive_over_traits':sum(x>1e-12 for x in real_gen_gain)},
      'site_average_K_gain':{'equal_site':avg_site_equal,'individual_weighted':avg_site_weight,
        'sites_positive_over_traits':sum(x['site_average_genomic_gain_over_traits']>1e-12 for x in per_site)},
      'identity_specific_increment_vs_site_average_K':{
        'equal_site':specific_equal,'individual_weighted':specific_weight,
        'sites_real_genomic_prediction_better_than_site_K':n_real_better},
      'within_site_genotype_identity_placebo':{
        'equal_site_genomic_gain':_stats(null_equal),
        'individual_weighted_genomic_gain':_stats(null_weight),
        'per_site_genomic_gain':{f['site']:_stats(a) for f,a in zip(folds,per_site_null)},
        'null_replicates_ge_observed_equal_site_gain':null_gte,
        'null_ge_observed_fraction_descriptive_only':null_gte/reps},
      'decision':'IDENTITY_SPECIFIC_INCREMENT_CLEARS_EXPLORATORY_ADVERSE_CONTROLS' if readiness else
                 'NO_ROBUST_INDIVIDUAL_GENOMIC_IDENTITY_INFORMATION_BEYOND_SITE_STRUCTURE',
      'full_ANCESTRAL_GENOTYPE_INTERACTION_FUTURE_FITNESS_TEST':False,
      'scientific_ceiling':'One original Ulex field study; same-site permutations are not experimentally valid genotypic randomization. No source pollinator ind, no measured pigment ancestry, no progeny and no causal selection estimates.'}


def run(payloads:dict[str,bytes],*,expected_v15:dict,reps:int=REPS,seed:int=SEED)->dict:
    """Never run the new placebo before the frozen v15 original receipt passes."""
    u.verify_payloads(payloads)
    recomputed=u.analyze(payloads)
    if recomputed != expected_v15:
        raise IdentityGateError('ORIGINAL_V15_FROZEN_RECEIPT_MISMATCH')
    fruit,_=u.parse_fruit_and_traits(payloads)
    K,_=u.parse_relatedness(payloads)
    result=analyze_rows(fruit,K,reps=reps,seed=seed,original_v15=recomputed)
    if result['matched_genotypes'] != recomputed['audit']['matched_genomic_flower_fruit_individuals']:
        raise IdentityGateError('MATCHED_GENOMIC_COHORT_DRIFT')
    return {
       'original_v15_frozen_receipt_recomputed_exactly':True,
       'source':{'original_study':'Castellanos et al. 2023 Ulex parviflorus',
                 'original_zenodo_record':u.RECORD_ID,
                 'source_MD5s_verified':dict(u.EXPECTED_MD5),
                 'source_read_only':True,'original_observation_rows_redistributed':False,
                 'pollinator_plant_id_available':False},
       **result}
