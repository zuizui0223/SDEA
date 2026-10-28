"""SDEA v11 offline adversarial checks for real-source descriptive drift.

All 4-state known counts below are bounded aggregated receipts already published
in the v10 fixture. NO original mother-level rows are bundled in this test.
"""
from __future__ import annotations

import copy
import json
import math
import unittest

from sdea.empirical_admission import SourceGateError
from sdea.empirical_transfer import CATEGORIES, TREATMENTS, _input_table, forecast
from sdea.empirical_drift import (
  LAMBDAS, shrink_association, _cell_distribution, _tv, observed_drift,
  choose_lambda_train_only, _train_cells_from_one_year, _heldout_scores,
  train_only_shrinkage_transfer, run_from_summary, run_original_payload,
)

COUNTS={
 '2007|ndni':(0,2,2,5),'2007|di':(0,1,6,4),
 '2007|ndi':(0,3,2,5),'2007|dni':(2,0,6,0),
 '2008|ndni':(12,1,1,1),'2008|di':(0,4,8,7),
 '2008|ndi':(0,7,1,13),'2008|dni':(5,0,7,0),
}
SOURCE={'repo':'aza3','commit':'pinned','path':'source.csv','blob_sha':'pinned-git-blob',
        'doi':'10.15482/USDA.ADC/1522649'}


def original_aggregate_summary():
    return {'source_mothers':106,'admitted_coherent_mothers':105,
            'excluded_incoherent':[{'year':'2007','plant_id':'327'}],
            'year_treatment_strata':{k:{'n':sum(counts),**dict(zip(CATEGORIES,counts))}
                                     for k,counts in COUNTS.items()}}


class TestEmpiricalDrift(unittest.TestCase):
    def setUp(self):
        self.table=_input_table(original_aggregate_summary())

    def test_lambda_zero_equals_original_independence(self):
        for counts in COUNTS.values():
            self.assertEqual(shrink_association(counts,0.),forecast(counts)['independent'])

    def test_lambda_one_equals_original_joint(self):
        for counts in COUNTS.values():
            self.assertEqual(shrink_association(counts,1.),forecast(counts)['joint'])

    def test_all_shrinkage_grid_predictions_normalize(self):
        for counts in COUNTS.values():
            for lam in LAMBDAS:
                p=shrink_association(counts,lam)
                self.assertAlmostEqual(sum(p),1.,places=13)
                self.assertTrue(all(0<x<1 for x in p))

    def test_shrinkage_changes_neither_marginal_at_any_lambda(self):
        for counts in COUNTS.values():
            reference=forecast(counts)['joint']
            for lam in LAMBDAS:
                x=shrink_association(counts,lam)
                self.assertAlmostEqual(x[1]+x[3],reference[1]+reference[3],places=13)
                self.assertAlmostEqual(x[2]+x[3],reference[2]+reference[3],places=13)

    def test_structurally_correlated_counts_are_shrunk(self):
        counts=(40,0,0,40)
        p0=shrink_association(counts,0.)
        pm=shrink_association(counts,.5)
        p1=shrink_association(counts,1.)
        self.assertLess(p0[3],pm[3]);self.assertLess(pm[3],p1[3])

    def test_independent_four_state_table_is_shrinkage_neutral(self):
        for lam in LAMBDAS:
            self.assertEqual(shrink_association((8,8,8,8),lam),(0.25,)*4)

    def test_invalid_lambda_raises(self):
        for lam in (-.01,1.0001,float('nan'),float('inf'),None,'0.5',True):
            with self.subTest(lam=str(lam)),self.assertRaises(SourceGateError):
                shrink_association((2,3,1,4),lam)

    def test_invalid_counts_raise(self):
        for counts in ((0,-1,3,4),(1.5,2,2,3),(1,2,3),None):
            with self.subTest(counts=str(counts)),self.assertRaises(SourceGateError):
                shrink_association(counts,.4)

    def test_zero_frequency_categories_finite_log_loss(self):
        for lam in LAMBDAS:
            p=shrink_association((10,0,0,0),lam)
            self.assertTrue(all(math.isfinite(math.log(x)) for x in p))

    def test_raw_covariance_vector_identity_each_cell(self):
        for counts in COUNTS.values():
            x=_cell_distribution(counts)
            c=x['covariance']
            self.assertEqual(x['association_residual'],(c,-c,-c,c))
            for i in range(4):
                self.assertAlmostEqual(x['joint_distribution'][i],
                  x['product_marginals'][i]+x['association_residual'][i])

    def test_frechet_bounds_hold_for_all_original_groups(self):
        for counts in COUNTS.values():
            x=_cell_distribution(counts)
            low,high=x['frechet_any_route']
            self.assertLessEqual(low-1e-12,x['observed_any_route'])
            self.assertLessEqual(x['observed_any_route'],high+1e-12)

    def test_tv_identity_and_triangle_bound_all_treatments(self):
        r=observed_drift(self.table)
        for row in r.values():
            q=row['year_difference_2008_minus_2007']
            self.assertAlmostEqual(sum(q['full_four_state_vector']),0.)
            for full,marg,assoc in zip(q['full_four_state_vector'],q['product_marginal_vector'],q['joint_association_vector']):
                self.assertAlmostEqual(full,marg+assoc,places=13)
            self.assertLessEqual(q['total_variation_joint'],q['total_variation_product_marginals']+q['total_variation_association_component']+1e-12)

    def test_total_variation_not_additive_no_false_causal_shares(self):
        x=observed_drift(self.table)['ndni']['year_difference_2008_minus_2007']
        self.assertNotAlmostEqual(x['total_variation_joint'],
                                  x['total_variation_product_marginals']+x['total_variation_association_component'])

    def test_ndni_observed_marginal_and_covariance_shifts(self):
        x=observed_drift(self.table)['ndni']['year_difference_2008_minus_2007']
        self.assertTrue(x['association_sign_changed'])
        self.assertAlmostEqual(x['apical_success'],-.6444444444444445)
        self.assertAlmostEqual(x['lateral_success'],-.6444444444444445)
        self.assertAlmostEqual(x['covariance'],.0982716049382716)
        self.assertAlmostEqual(x['total_variation_joint'],.8)
        self.assertAlmostEqual(x['total_variation_product_marginals'],.7017283950617285)

    def test_dni_has_no_association_drift_with_apical_never(self):
        x=observed_drift(self.table)['dni']['year_difference_2008_minus_2007']
        self.assertEqual(x['covariance'],0.)
        self.assertAlmostEqual(x['total_variation_association_component'],0.)
        self.assertGreater(x['total_variation_joint'],0.)

    def test_invalid_descriptive_cell_rejected(self):
        for x in ((0,0,0,0),(-1,1,1,1),(0.,1,2,3),(1,2,3),None):
            with self.assertRaises(SourceGateError):_cell_distribution(x)

    def test_training_lambda_selected_only_by_own_year(self):
        r=train_only_shrinkage_transfer(self.table)
        self.assertEqual(r['folds']['2007_to_2008']['model_selection']['selected_lambda'],.25)
        self.assertEqual(r['folds']['2008_to_2007']['model_selection']['selected_lambda'],.75)

    def test_training_year_2007_lower_loo_score_at_quarter(self):
        tr=_train_cells_from_one_year(self.table,'2007')
        result=choose_lambda_train_only(tr)
        x={z['lambda']:z['train_only_loo_logloss'] for z in result['grid']}
        self.assertLess(x[.25],x[0])
        self.assertLess(x[.25],x[1])
        self.assertEqual(result['training_mothers'],38)

    def test_training_year_2008_lower_loo_score_at_three_quarters(self):
        tr=_train_cells_from_one_year(self.table,'2008')
        result=choose_lambda_train_only(tr)
        x={z['lambda']:z['train_only_loo_logloss'] for z in result['grid']}
        self.assertLess(x[.75],x[0])
        self.assertLess(x[.75],x[1])
        self.assertEqual(result['training_mothers'],67)

    def test_target_year_mutation_cannot_affect_lambda_selection(self):
        original=choose_lambda_train_only(_train_cells_from_one_year(self.table,'2007'))
        manipulated=copy.deepcopy(self.table)
        manipulated['2008|ndni']=(0,0,0,15)
        modified=choose_lambda_train_only(_train_cells_from_one_year(manipulated,'2007'))
        self.assertEqual(original,modified)
        self.assertNotEqual(train_only_shrinkage_transfer(self.table)['folds']['2007_to_2008']['heldout_selected'],
                            train_only_shrinkage_transfer(manipulated)['folds']['2007_to_2008']['heldout_selected'])

    def test_lambda_grid_validation(self):
        tr=_train_cells_from_one_year(self.table,'2007')
        for grid in ((),(.2,.2),[0.,1.],(.2,1.2),(float('nan'),)):
            with self.assertRaises(SourceGateError):choose_lambda_train_only(tr,grid)

    def test_training_key_validation(self):
        train=_train_cells_from_one_year(self.table,'2008')
        del train['ndi']
        with self.assertRaises(SourceGateError):choose_lambda_train_only(train)

    def test_original_joint_same_as_v10_heldout(self):
        r=train_only_shrinkage_transfer(self.table)
        f=r['folds']['2007_to_2008']
        self.assertAlmostEqual(f['full_joint_minus_independent']['logloss_per_mother'],.057425013941451525)
        self.assertAlmostEqual(f['full_joint_minus_independent']['categorical_brier_per_mother'],.009949371004330848)
        f=r['folds']['2008_to_2007']
        self.assertAlmostEqual(f['full_joint_minus_independent']['logloss_per_mother'],-.11790584215462149)
        self.assertAlmostEqual(f['full_joint_minus_independent']['categorical_brier_per_mother'],.0009125754813731989)

    def test_shrunken_joint_still_fails_2007_to_2008_both_scores(self):
        f=train_only_shrinkage_transfer(self.table)['folds']['2007_to_2008']
        self.assertGreater(f['selected_minus_independent']['logloss_per_mother'],0.)
        self.assertGreater(f['selected_minus_independent']['categorical_brier_per_mother'],0.)

    def test_shrunken_joint_2008_to_2007_improves_both(self):
        f=train_only_shrinkage_transfer(self.table)['folds']['2008_to_2007']
        self.assertLess(f['selected_minus_independent']['logloss_per_mother'],0.)
        self.assertLess(f['selected_minus_independent']['categorical_brier_per_mother'],0.)

    def test_two_year_global_transfer_gate_rejects_one_direction(self):
        r=train_only_shrinkage_transfer(self.table)
        self.assertEqual(r['heldout_mother_years'],105)
        self.assertEqual(r['transport_gate'],'NO_REPLICATED_TRAIN_ONLY_SHRINKAGE_TRANSFER_GAIN')

    def test_no_leakage_tie_break_selects_smallest_lambda(self):
        tr={t:(5,5,5,5) for t in TREATMENTS}
        # All candidate distributions coincide when every treatment is exactly independent;
        # LOO creates a weak nonindependence for the single excluded category, so
        # verify a one-element grid instead of overclaiming a tie here.
        r=choose_lambda_train_only(tr,grid=(0.,))
        self.assertEqual(r['selected_lambda'],0.)

    def test_reproducible_report_and_no_original_mother_rows(self):
        summary=original_aggregate_summary()
        a=run_from_summary(summary,mother_entry=SOURCE)
        b=run_from_summary(copy.deepcopy(summary),mother_entry=SOURCE)
        self.assertEqual(a,b)
        self.assertEqual(a['scientific_decision'],
          'REAL_MARGINAL_AND_ASSOCIATION_DRIFT_NO_REPLICATED_SHRINKAGE_TRANSFER')
        self.assertGreater(len(json.dumps(a)),7000)
        self.assertNotIn('original_plant_rows',a)

    def test_bad_original_source_checksum_fails_closed(self):
        with self.assertRaises(SourceGateError):run_original_payload(b'this is not the source',mother_entry=SOURCE)

    def test_drift_is_not_defined_for_missing_strata(self):
        t=copy.deepcopy(self.table)
        t.pop('2008|ndni')
        with self.assertRaises(SourceGateError):observed_drift(t)


if __name__=='__main__':unittest.main()
