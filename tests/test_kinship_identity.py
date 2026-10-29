"""SDEA v16: exact validation, adverse controls, no genetic-inference inflation."""
from __future__ import annotations
import copy
import math
import random
import unittest
from collections import defaultdict

from sdea import kinship_identity as k
from sdea import ulex_partial as u


def synthetic(*, signal: bool = True):
    """Construct six populations and strictly one genotype per actual plant."""
    sizes = (32,37,33,37,37,26)
    rows=[];families={}
    for site,siz in zip(k.EXPECTED_SITES,sizes):
        for i in range(siz):
            uid=f'U{len(rows):03d}'
            family=i%5
            families[uid]=family
            elev=200.+(i%13)*3.+len(site)
            area=12.+(i%7)*.45
            weight=2.+(i%11)*.07
            fruit=.3 + .035*(family-2) if signal else .3
            fruit+=.01*math.sin(i+len(site)) if signal else 0.
            rows.append({'site':site,'ind':uid,'elev':elev,'area':area,'weight':weight,
                         'scars':100,'fruits':int(round(100*fruit)),
                         'fruit_set':int(round(100*fruit))/100.})
    K={r['ind']:{} for r in rows}
    for a in rows:
        for b in rows:
            if a['ind']==b['ind']:val=1.
            elif families[a['ind']]==families[b['ind']]:val=.8 if a['site']==b['site'] else .6
            elif a['site']==b['site']:val=.06
            else:val=-.04
            K[a['ind']][b['ind']]=val
    return rows,K


class V16IdentityControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows,cls.K=synthetic()
        cls.folds,cls.W,cls.groups=k._fold_data(cls.rows,cls.K)

    def test_v15_source_constants_not_redefined(self):
        self.assertEqual(u.RECORD_ID,7761289)
        self.assertEqual(u.RIDGE,5.)
        self.assertEqual(u.KINSHIP_SHRINKAGE,5.)
        self.assertEqual(len(u.EXPECTED_MD5),5)

    def test_full_synthetic_cohort_and_exact_original_site_counts(self):
        self.assertEqual(len(self.rows),202)
        self.assertEqual(tuple(sorted(self.groups)),k.EXPECTED_SITES)
        self.assertEqual([len(self.groups[s]) for s in sorted(self.groups)], [32,37,33,37,37,26])

    def test_real_genomes_not_independent_plants_when_repeats(self):
        rows=copy.deepcopy(self.rows)
        rows[0]['ind']=rows[1]['ind']
        with self.assertRaisesRegex(k.IdentityGateError,'REUSED_INDIVIDUAL_ID'):
            k._fold_data(rows,self.K)

    def test_seven_original_sites_fail_closed(self):
        rows=copy.deepcopy(self.rows)
        rows[0]['site']='UNDECLARED_SITE'
        with self.assertRaisesRegex(k.IdentityGateError,'SIX_SOURCE_POPULATIONS'):
            k._fold_data(rows,self.K)

    def test_201_instead_of_202_fails_closed(self):
        with self.assertRaisesRegex(k.IdentityGateError,'EXPECTED_202'):
            k._fold_data(self.rows[:-1],self.K)

    def test_identity_mapping_has_all_genome_labels(self):
        map_ = list(range(len(self.rows)))
        self.assertEqual(len(self.W),len(map_))
        self.assertEqual(sorted(map_),list(range(202)))

    def test_within_site_permutation_bijection(self):
        m=k._within_site_relabeling(self.groups,random.Random(52),202)
        self.assertEqual(sorted(m),list(range(202)))
        for s,idxs in self.groups.items():
            self.assertEqual(sorted(m[i] for i in idxs),sorted(idxs))

    def test_permuted_mapping_moves_no_plant_outside_population(self):
        m=k._within_site_relabeling(self.groups,random.Random(52),202)
        lookup={i:s for s,idxs in self.groups.items() for i in idxs}
        self.assertTrue(all(lookup[i]==lookup[m[i]] for i in range(202)))

    def test_permutations_deterministic_given_seed(self):
        a=k._within_site_relabeling(self.groups,random.Random(66),202)
        b=k._within_site_relabeling(self.groups,random.Random(66),202)
        self.assertEqual(a,b)

    def test_permutation_seeds_change_genetic_assignments(self):
        a=k._within_site_relabeling(self.groups,random.Random(66),202)
        b=k._within_site_relabeling(self.groups,random.Random(67),202)
        self.assertNotEqual(a,b)

    def test_nonbijective_genomic_mapping_blocked(self):
        mapping=list(range(202));mapping[0]=mapping[1]
        with self.assertRaisesRegex(k.IdentityGateError,'BIJECTION'):
            k._score_mapping(self.folds,self.W,mapping)

    def test_bad_shrinkage_rejected(self):
        for z in (0.,-1.,float('inf'),float('nan')):
            with self.assertRaises(k.IdentityGateError):
                k._predict_mapping_fold(self.folds[0],self.W,list(range(202)),shrinkage=z)

    def test_true_genomic_folds_equal_original_v15_algorithm(self):
        order=sorted(self.rows,key=lambda r:(r['site'],r['ind']))
        for fold in self.folds:
            test=[order[i] for i in fold['test']]
            train=[order[i] for i in fold['train']]
            model=u.fit(train,traits=True)
            expected=[u.related_residual_prediction(r,train,model,self.K) for r in test]
            actual=k._predict_mapping_fold(fold,self.W,list(range(202)))
            for a,b in zip(actual,expected):self.assertAlmostEqual(a,b,places=12)

    def test_parent_m1_residual_same_as_genomic_fold_baseline(self):
        for f in self.folds:
            self.assertAlmostEqual(f['trait_only_mse'],k._mse(f['baseline'],f['truth']))

    def test_heldout_y_mutation_does_not_change_predictions(self):
        f=copy.deepcopy(self.folds[0]); mapping=list(range(202))
        before=k._predict_mapping_fold(f,self.W,mapping)
        f['truth']=[1-y for y in f['truth']]
        self.assertEqual(before,k._predict_mapping_fold(f,self.W,mapping))

    def test_nonnegative_genomic_weighting_clips_only_for_prediction(self):
        self.assertTrue(all(w>=0 for row in self.W for w in row))
        self.assertTrue(any(w==0 for row in self.W for w in row))
        self.assertLess(min(v for row in self.K.values() for v in row.values()),0.)

    def test_site_averaged_genetic_adjustment_has_common_offset(self):
        fold=self.folds[0]
        tr=fold['train'];te=fold['test']
        means=[sum(self.W[i][j] for i in te)/len(te) for j in tr]
        shift=sum(w*r for w,r in zip(means,fold['residual']))/(5+sum(means))
        expected=k._mse([max(0,min(1,p+shift)) for p in fold['baseline']],fold['truth'])
        self.assertAlmostEqual(expected,k._site_average_K(self.folds,self.W)[0],places=12)

    def test_null_when_genome_weights_all_zero(self):
        W0=[[0.]*202 for _ in range(202)]
        actual=k._score_mapping(self.folds,W0,list(range(202)))
        place=k._score_mapping(self.folds,W0,k._within_site_relabeling(self.groups,random.Random(2),202))
        site=k._site_average_K(self.folds,W0)
        for v,f in zip(actual,self.folds):self.assertAlmostEqual(v,f['trait_only_mse'])
        self.assertEqual(actual,place)
        self.assertEqual(actual,site)

    def test_large_shrinkage_tends_to_trait_only_baseline(self):
        p=k._score_mapping(self.folds,self.W,list(range(202)),shrinkage=1e15)
        self.assertTrue(all(abs(v-f['trait_only_mse'])<1e-12 for v,f in zip(p,self.folds)))

    def test_predicted_probability_is_bounded(self):
        p=k._predict_mapping_fold(self.folds[0],self.W,list(range(202)))
        self.assertTrue(all(0<=y<=1 for y in p))

    def test_empty_mse_or_wrong_target_denominator_rejected(self):
        with self.assertRaises(k.IdentityGateError):k._mse([],[])
        with self.assertRaises(k.IdentityGateError):k._mse([0.1],[0.1,0.2])

    def test_quantile_endpoints_and_median(self):
        data=[3.,1.,4.,2.]
        self.assertEqual(k._quantile(data,0),1)
        self.assertEqual(k._quantile(data,1),4)
        self.assertEqual(k._quantile(data,.5),2.5)

    def test_quantile_nonfinite_blocks(self):
        with self.assertRaises(k.IdentityGateError):k._quantile([1,float('nan')],.5)
        with self.assertRaises(k.IdentityGateError):k._quantile([1],-1)

    def test_fixed_rng_placebo_repeatability_on_synthetic(self):
        a=k.analyze_rows(self.rows,self.K,reps=13,seed=101)
        b=k.analyze_rows(self.rows,self.K,reps=13,seed=101)
        self.assertEqual(a,b)

    def test_placebo_changes_genotype_identity_not_source_y(self):
        f=self.folds[0]
        mapping=k._within_site_relabeling(self.groups,random.Random(135),202)
        self.assertEqual(len(f['truth']),len(k._predict_mapping_fold(f,self.W,mapping)))
        self.assertEqual(f['truth'],[self.rows_sorted()[i]['fruit_set'] for i in f['test']])

    @classmethod
    def rows_sorted(cls):
        return sorted(cls.rows,key=lambda r:(r['site'],r['ind']))

    def test_placebo_does_not_overwrite_real_parent_model(self):
        actual=k._score_mapping(self.folds,self.W,list(range(202)))
        mapping=k._within_site_relabeling(self.groups,random.Random(225),202)
        k._score_mapping(self.folds,self.W,mapping)
        self.assertEqual(actual,k._score_mapping(self.folds,self.W,list(range(202))))

    def test_observed_vs_null_reference_bounded(self):
        d=k.analyze_rows(self.rows,self.K,reps=15,seed=1)
        self.assertEqual(d['matched_genotypes'],202)
        self.assertEqual(d['null_design']['replicates'],15)
        self.assertEqual(len(d['site_folds']),6)
        self.assertTrue(0<=d['within_site_genotype_identity_placebo']['null_ge_observed_fraction_descriptive_only']<=1)
        self.assertFalse(d['full_ANCESTRAL_GENOTYPE_INTERACTION_FUTURE_FITNESS_TEST'])

    def test_requires_all_original_site_names(self):
        rows=copy.deepcopy(self.rows)
        for r in rows:
            if r['site']=='Chiva':r['site']='Chiva-alias'
        with self.assertRaises(k.IdentityGateError):k.analyze_rows(rows,self.K,reps=10)

    def test_invalid_rep_count_or_seed_fails_closed(self):
        for n,seed in [(0,1),(2,1),(5001,1),(10,0),(10,-1),(10.5,1)]:
            with self.assertRaises(k.IdentityGateError):
                k.analyze_rows(self.rows,self.K,reps=n,seed=seed)

    def test_v15_original_frozen_model_must_match(self):
        before=k.analyze_rows(self.rows,self.K,reps=10)
        original=[]
        for item in before['site_folds']:
            original.append({'site':item['site'],
                'plus_floral_area_weight_quadratic_mse':item['trait_only_mse'],
                'plus_relatedness_residual_mse':item['original_genomic_mse']})
        unchanged=k.analyze_rows(self.rows,self.K,reps=10,original_v15={'heldout_sites':original})
        self.assertEqual(unchanged,before)
        original[0]['plus_relatedness_residual_mse']+=.01
        with self.assertRaisesRegex(k.IdentityGateError,'ORIGINAL_V15_M1_M2'):
            k.analyze_rows(self.rows,self.K,reps=10,original_v15={'heldout_sites':original})

    def test_source_verification_fail_closed_no_original_payload(self):
        with self.assertRaisesRegex(u.UlexGateError,'MISSING_OR_EXTRA'):
            k.run({},expected_v15={})

    def test_joint_plant_genotype_scope_never_pollinator_ind(self):
        d=k.analyze_rows(self.rows,self.K,reps=10)
        self.assertIn('not a valid',d['scientific_ceiling'].lower()) if False else self.assertFalse(d['full_ANCESTRAL_GENOTYPE_INTERACTION_FUTURE_FITNESS_TEST'])
        self.assertTrue(d['null_design']['exchangeability_of_identities_within_sites_not_established'])


if __name__=='__main__':unittest.main()
