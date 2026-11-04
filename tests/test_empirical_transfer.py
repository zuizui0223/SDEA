"""v10: strict source, 4-state forecast, no-leakage and transfer negative tests."""
import copy
import math
import unittest

from sdea.empirical_admission import SourceGateError
from sdea.empirical_transfer import (
    ALPHA, CATEGORIES, TREATMENTS, YEARS, _input_table, _score,
    association_by_stratum, forecast, heldout_forecast,
    leave_one_mother_out_within_year, run_from_summary, run_original_payload,
)

# Bounded known source-reported four-state counts from the v9 original replay.
# Refit v10 on original pinned CSV in CI, NEVER treat these as new raw data.
KNOWN={
 '2007|ndni': (0,2,2,5), '2007|di': (0,1,6,4),
 '2007|ndi': (0,3,2,5), '2007|dni': (2,0,6,0),
 '2008|ndni': (12,1,1,1),'2008|di': (0,4,8,7),
 '2008|ndi': (0,7,1,13),'2008|dni': (5,0,7,0),
}
MOTHER={'repo':'aza3','commit':'locked-commit','path':'source.csv','blob_sha':'locked-blob','doi':'source-doi'}


def make_summary():
    return {'source_mothers':106,'admitted_coherent_mothers':105,
            'excluded_incoherent':[{'year':'2007','plant_id':'327'}],
            'year_treatment_strata':{key:{'n':sum(counts),**dict(zip(CATEGORIES,counts))}
                                     for key,counts in KNOWN.items()}}


class TransferTests(unittest.TestCase):
    def setUp(self):
        self.summary=make_summary()
        self.table=_input_table(self.summary)

    def test_source_counts_and_strata_locked(self):
        self.assertEqual(sum(sum(x) for x in self.table.values()),105)
        self.assertEqual(len(self.table),8)
        self.assertEqual(sum(sum(v) for k,v in self.table.items() if k.startswith('2007|')),38)
        self.assertEqual(sum(sum(v) for k,v in self.table.items() if k.startswith('2008|')),67)

    def test_source_bad_record_is_explicit(self):
        self.summary['excluded_incoherent']=[{'year':'2007','plant_id':'355'}]
        with self.assertRaises(SourceGateError): _input_table(self.summary)

    def test_source_admission_count_changed_fails(self):
        self.summary['admitted_coherent_mothers']=106
        with self.assertRaises(SourceGateError): _input_table(self.summary)

    def test_missing_or_extra_year_stratum_fails(self):
        self.summary['year_treatment_strata'].pop('2008|ndi')
        with self.assertRaises(SourceGateError): _input_table(self.summary)
        self.summary=make_summary()
        self.summary['year_treatment_strata']['2009|ndi']=self.summary['year_treatment_strata']['2008|ndi']
        with self.assertRaises(SourceGateError): _input_table(self.summary)

    def test_missing_or_incoherent_four_state_cell_fails(self):
        cell=self.summary['year_treatment_strata']['2008|ndni']
        del cell['both']
        with self.assertRaises(SourceGateError): _input_table(self.summary)
        self.summary=make_summary()
        self.summary['year_treatment_strata']['2008|ndni']['both']=4
        with self.assertRaises(SourceGateError): _input_table(self.summary)

    def test_fitted_marginals_exactly_identical(self):
        for x in ((0,1,6,4),(12,1,1,1),(0,0,0,0),(7,0,5,0),(0,4,8,7),(3,1,2,10)):
            p=forecast(x)
            self.assertAlmostEqual(p['joint'][1]+p['joint'][3],p['independent'][1]+p['independent'][3])
            self.assertAlmostEqual(p['joint'][2]+p['joint'][3],p['independent'][2]+p['independent'][3])
            self.assertAlmostEqual(sum(p['joint']),1.)
            self.assertAlmostEqual(sum(p['independent']),1.)

    def test_balanced_independence_exact_control(self):
        p=forecast((5,5,5,5))
        self.assertEqual(p['independent'],p['joint'])

    def test_structurally_nonindependent_same_marginals(self):
        p=forecast((20,0,0,20))
        self.assertAlmostEqual(p['joint'][1]+p['joint'][3],.5)
        self.assertGreater(p['joint'][3],p['independent'][3])

    def test_nonzero_smoothing_prevents_infinite_loss(self):
        p=forecast((12,0,0,0))
        for name in ('independent','joint'):
            self.assertTrue(all(x>0 for x in p[name]))
            self.assertTrue(math.isfinite(_score((0,0,0,1),p[name])['logloss_sum']))

    def test_malformed_counts_and_priors_fail(self):
        for x in ((1,2,3),(0,-1,1,0),(0.,1,0,0),None,('x',0,0,0)):
            with self.assertRaises(SourceGateError): forecast(x)
        for alpha in (0.,-1.,float('nan'),float('inf')):
            with self.assertRaises(SourceGateError): forecast((2,1,3,4),alpha=alpha)

    def test_score_manual_logloss_and_brier(self):
        p=(.25,)*4
        x=_score((1,0,0,0),p)
        self.assertAlmostEqual(x['logloss_sum'],math.log(4))
        self.assertAlmostEqual(x['brier_sum'],.75)

    def test_original_covariance_exact_identity_eight_cells(self):
        assoc=association_by_stratum(self.table)
        for z in assoc.values():
            self.assertAlmostEqual(z['observed_minus_independent'], -z['same_mother_covariance'],places=13)

    def test_ndni_covariance_reverses_year(self):
        assoc=association_by_stratum(self.table)
        self.assertLess(assoc['2007|ndni']['same_mother_covariance'],0.)
        self.assertGreater(assoc['2008|ndni']['same_mother_covariance'],0.)

    def test_year_holdout_no_training_leakage(self):
        a=heldout_forecast(self.table)
        changed=dict(self.table)
        changed['2008|ndni']=(0,0,0,15)  # mutate held-out labels, not training labels
        b=heldout_forecast(changed)
        ap=a['folds']['2007_to_2008']['by_treatment']
        bp=b['folds']['2007_to_2008']['by_treatment']
        for x,y in zip(ap,bp):
            self.assertEqual(x['prediction'],y['prediction'])
        self.assertNotEqual(a['folds']['2007_to_2008']['model_scores'],b['folds']['2007_to_2008']['model_scores'])

    def test_reversed_year_holdout_changes_training(self):
        o=heldout_forecast(self.table)
        a=o['folds']['2007_to_2008']['by_treatment'][0]['prediction']
        b=o['folds']['2008_to_2007']['by_treatment'][0]['prediction']
        self.assertNotEqual(a,b)

    def test_holdout_mother_denominators(self):
        r=heldout_forecast(self.table)
        self.assertEqual(r['folds']['2007_to_2008']['heldout_mothers'],67)
        self.assertEqual(r['folds']['2008_to_2007']['heldout_mothers'],38)
        self.assertEqual(r['tested_mother_year_units'],105)

    def test_joint_transfer_primary_decision_fails_preset(self):
        r=heldout_forecast(self.table)
        self.assertEqual(r['predeclared_transfer_gate'],'NO_REPLICATED_JOINT_FORECAST_GAIN')
        self.assertGreater(r['folds']['2007_to_2008']['joint_minus_independent']['logloss_per_mother'],0)
        self.assertLess(r['folds']['2008_to_2007']['joint_minus_independent']['logloss_per_mother'],0)
        self.assertGreater(r['folds']['2007_to_2008']['joint_minus_independent']['categorical_brier_per_mother'],0)
        self.assertGreater(r['folds']['2008_to_2007']['joint_minus_independent']['categorical_brier_per_mother'],0)

    def test_pooled_direction_cannot_overrule_both_folds(self):
        r=heldout_forecast(self.table)
        self.assertLess(r['pooled_joint_minus_independent']['logloss_per_mother'],0)
        self.assertGreater(r['pooled_joint_minus_independent']['categorical_brier_per_mother'],0)
        self.assertEqual(r['predeclared_transfer_gate'],'NO_REPLICATED_JOINT_FORECAST_GAIN')

    def test_secondary_same_year_loomo_different_estimand(self):
        r=leave_one_mother_out_within_year(self.table)
        self.assertEqual(r['2007']['heldout_mothers'],38)
        self.assertEqual(r['2008']['heldout_mothers'],67)
        self.assertGreater(r['2007']['joint_minus_independent']['logloss_per_mother'],0.)
        self.assertLess(r['2008']['joint_minus_independent']['logloss_per_mother'],0.)

    def test_integration_retains_biological_stop(self):
        r=run_from_summary(self.summary,mother_entry=MOTHER)
        self.assertEqual(r['biological_decision'],'OBSERVED_SAME_MOTHER_ROUTE_DEPENDENCE_BUT_NO_REPLICATED_YEAR_TRANSFER')
        self.assertEqual(r['source']['source_mothers'],106)
        self.assertEqual(r['source']['admitted_coherent_mothers'],105)
        self.assertTrue(r['design']['same_marginals_in_both_predictive_models'])

    def test_changed_original_source_sha_rejected(self):
        with self.assertRaises(SourceGateError): run_original_payload(b'non-original,maternal,source',mother_entry=MOTHER)

    def test_invalid_score_probabilities_rejected(self):
        for p in ((0.,.5,.2,.3),(float('nan'),.2,.3,.5),(.2,.3)):
            with self.assertRaises(SourceGateError):_score((1,0,0,0),p)

    def test_deterministic_report_and_json_serializable(self):
        import json
        x=run_from_summary(self.summary,mother_entry=MOTHER)
        y=run_from_summary(copy.deepcopy(self.summary),mother_entry=MOTHER)
        self.assertEqual(x,y)
        self.assertGreater(len(json.dumps(x)),1000)

if __name__=='__main__':unittest.main()
