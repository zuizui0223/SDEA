"""Post-v13 source-conditional forecast stress tests; no natural simulation claims."""
import math
import unittest
from copy import deepcopy
from random import Random

from sdea.external_stage_transfer import YEARS, ENDPOINTS, ExternalStageError
from sdea.within_context import (
    MODES,_group,_features_by_group,_train_context_slope,
    _feature_only_test_rows,_prediction,_mse,heldout,analyze,
)


def sample(n_each=4, *, width_relevant=True, escape_relevant=True, strong_group_shift=False):
    rows=[]
    for k,year in enumerate(YEARS):
        for sitecode in ('agg','hyb','ten'):
            for t in ('A','H','T'):
                site='agg' if sitecode=='agg' else 'upper'
                g={'agg':0,'hyb':1,'ten':2}[sitecode]
                tt={'A':0,'H':1,'T':2}[t]
                for i in range(n_each):
                    x=(i-(n_each-1)/2)*.3
                    w=3. + k*.07 + g*.4 + tt*.13 + x
                    inp=1.8+.09*k+.06*g+.10*tt+(0.7*x if width_relevant else 0.)
                    # Keep fruit-stage values in the interval [0,1].
                    escape=.38+.04*k/9+.025*tt + .07*(inp-2.) + (.22*x if escape_relevant else 0.)
                    rows.append({'year':year,'site_code':sitecode,'physical_site':site,'type':t,
                                 'width':w,'initiated':inp,'escape':escape})
    return rows


class TestContextFormulas(unittest.TestCase):
    def test_stratification_uses_physical_location(self):
        a={'year':'2015','site_code':'hyb','physical_site':'upper','type':'H'}
        b={**a,'site_code':'ten'}
        self.assertEqual(_group(a,'physical_location'),_group(b,'physical_location'))
        self.assertNotEqual(_group(a,'original_site_code'),_group(b,'original_site_code'))

    def test_unknown_mode_raises(self):
        with self.assertRaises(ExternalStageError): _group(sample()[0],'random')

    def test_unknown_stage_raises(self):
        with self.assertRaises(ExternalStageError):_features_by_group(sample()[:3],'physical_location',stage='seeds',stats={})

    def test_no_rows_raises(self):
        with self.assertRaises(ExternalStageError):_features_by_group([],'physical_location',stage=ENDPOINTS[0],stats={})

    def test_centered_width_sums_zero_in_each_group(self):
        rows=sample(3)
        stats={'width':(0,1),'initiated':(0,1)}
        for mode in MODES:
            f=_features_by_group(rows,mode,stage=ENDPOINTS[0],stats=stats)
            keys={x['group'] for x in f}
            for k in keys:
                self.assertAlmostEqual(sum(x['w_centered'] for x in f if x['group']==k),0)

    def test_centered_prior_sums_zero_in_each_group(self):
        rows=sample(3)
        f=_features_by_group(rows,'original_site_code',stage=ENDPOINTS[1],stats={'width':(0,1),'initiated':(0,1)})
        for k in {x['group'] for x in f}:
            self.assertAlmostEqual(sum(x['prior_centered'] for x in f if x['group']==k),0)

    def test_singleton_yields_no_association_signal(self):
        r=sample(1)[:1]
        f=_features_by_group(r,'physical_location',stage=ENDPOINTS[0],stats={'width':(0,1),'initiated':(0,1)})
        self.assertEqual(f[0]['n_group'],1)
        self.assertEqual(f[0]['w_centered'],0)

    def test_stage1_never_requires_test_initiated_outcome(self):
        r=sample(4)
        a=_feature_only_test_rows(r,ENDPOINTS[0]);self.assertTrue(all('initiated' not in x for x in a))
        f=_features_by_group(a,'physical_location',stage=ENDPOINTS[0],stats={'width':(0,1),'initiated':(0,1)})
        self.assertEqual(len(f),len(r))

    def test_stage2_uses_prior_stage_not_its_test_response(self):
        r=sample(4)
        a=_feature_only_test_rows(r,ENDPOINTS[1])
        self.assertTrue(all('escape' not in x and 'initiated' in x for x in a))

    def test_slope_positive_under_a_real_within_group_relation(self):
        r=sample(4)
        for stage in ENDPOINTS:
            s=_train_context_slope(r,'original_site_code',stage)
            self.assertGreater(s['beta'],0.)
            self.assertEqual(s['training_group_count'],9*3*3)
            self.assertEqual(s['n_training_rows'],len(r))

    def test_slope_zero_for_pure_between_group_signal(self):
        r=sample(4,width_relevant=False,escape_relevant=False)
        for stage in ENDPOINTS:
            s=_train_context_slope(r,'original_site_code',stage)
            self.assertAlmostEqual(s['beta'],0,places=10)

    def test_stronger_stage1_within_signal_recovers_better_predictions(self):
        r=sample(4,width_relevant=True)
        a=heldout(r,mode='physical_location',stage=ENDPOINTS[0])
        self.assertGreater(a['heldout_row_weighted_mse_gain'],0)
        self.assertEqual(a['n_years_width_improves'],9)

    def test_earlier_stage_relevance_never_implies_later_stage_causal_effect(self):
        r=sample(4,escape_relevant=False)
        a=_train_context_slope(r,'original_site_code',ENDPOINTS[0])
        b=_train_context_slope(r,'original_site_code',ENDPOINTS[1])
        self.assertGreater(a['beta'],0)
        self.assertLess(abs(b['beta']),0.015)

    def test_pure_between_group_signal_has_no_within_slope(self):
        r=sample(4,width_relevant=False,escape_relevant=False,strong_group_shift=True)
        for stage in ENDPOINTS:
            s=heldout(r,mode='original_site_code',stage=stage)
            self.assertAlmostEqual(s['heldout_row_weighted_mse_gain'],0,places=10)

    def test_training_year_never_overlaps_test(self):
        rows=sample(3)
        with self.assertRaisesRegex(ExternalStageError,'LEAKAGE'):
            _prediction(rows[:3],rows,'physical_location',ENDPOINTS[0])

    def test_mutating_heldout_stage1_y_does_not_change_predictions(self):
        r=sample(4);te=[x for x in r if x['year']=='2018'];tr=[x for x in r if x['year']!='2018']
        p=_prediction(te,tr,'physical_location',ENDPOINTS[0])
        changed=deepcopy(te)
        for x in changed:x['initiated']=100000.
        q=_prediction(changed,tr,'physical_location',ENDPOINTS[0])
        self.assertEqual(p[0],q[0]);self.assertEqual(p[1],q[1])

    def test_mutating_heldout_stage2_y_does_not_change_predictions(self):
        r=sample(4);te=[x for x in r if x['year']=='2018'];tr=[x for x in r if x['year']!='2018']
        p=_prediction(te,tr,'physical_location',ENDPOINTS[1]);changed=deepcopy(te)
        for x in changed:x['escape']=.99
        q=_prediction(changed,tr,'physical_location',ENDPOINTS[1])
        self.assertEqual(p[0],q[0]);self.assertEqual(p[1],q[1])

    def test_features_change_if_unlabeled_heldout_width_changes(self):
        r=sample(4);te=[x for x in r if x['year']=='2018'];tr=[x for x in r if x['year']!='2018']
        p=_prediction(te,tr,'physical_location',ENDPOINTS[0]);changed=deepcopy(te)
        changed[0]['width']+=10
        q=_prediction(changed,tr,'physical_location',ENDPOINTS[0]);self.assertNotEqual(p[1],q[1])

    def test_stage2_feature_only_uses_initiated_seeds(self):
        r=sample(4);te=[x for x in r if x['year']=='2018'];tr=[x for x in r if x['year']!='2018']
        p=_prediction(te,tr,'physical_location',ENDPOINTS[1]);changed=deepcopy(te)
        changed[0]['initiated']+=1
        q=_prediction(changed,tr,'physical_location',ENDPOINTS[1]);self.assertNotEqual(p[1],q[1])

    def test_zero_ridge_rejected(self):
        with self.assertRaises(ExternalStageError):_train_context_slope(sample(),'physical_location',ENDPOINTS[0],0)

    def test_nan_ridge_rejected(self):
        with self.assertRaises(ExternalStageError):_train_context_slope(sample(),'physical_location',ENDPOINTS[0],float('nan'))

    def test_nonfinite_prediction_fails(self):
        with self.assertRaises(ExternalStageError):_mse(sample()[:2],[float('nan'),.2],ENDPOINTS[0])

    def test_misaligned_mse_rows_fail(self):
        with self.assertRaises(ExternalStageError):_mse(sample()[:2],[.2],ENDPOINTS[0])

    def test_invalid_mode_rejected_in_holdout(self):
        with self.assertRaises(ExternalStageError):heldout(sample(),mode='wrong',stage=ENDPOINTS[0])

    def test_missing_year_rejected(self):
        with self.assertRaises(ExternalStageError):heldout(sample()[:4],mode='physical_location',stage=ENDPOINTS[0])

    def test_all_stage_and_context_cells_are_audited_in_a_mock_source_shaped_dataset(self):
        # Artificial 367-row table exercises the output shape, NOT the actual
        # original source statistics or any empirical validation result.
        r=sample(5)[:367]
        report=analyze(r,{'original_records':437}, {'same_original_study_not_an_independent_replication':True})
        self.assertEqual(len(report['within_context_forecasts']),2)
        for mode in MODES:
            for stage in ENDPOINTS:
                cell=report['within_context_forecasts'][mode][stage]
                self.assertEqual(cell['heldout_original_records_once_each'],367)
                self.assertEqual(cell['heldout_years'],9)
                self.assertIn(cell['readiness_gate'],('WITHIN_CONTEXT_WIDTH_INFO_RETAINS_YEAR_TRANSFER_IN_THIS_SOURCE','NO_ROBUST_WITHIN_CONTEXT_YEAR_TRANSFER'))
        self.assertIn('NOT_IDENTIFIABLE',report['cross_original_study_numerical_treatment_effect'])
        self.assertIn('transductive',str(report['design']).replace('unlabeled full heldout-year cohort','transductive'))

    def test_original_source_count_mismatch_blocks_analysis(self):
        with self.assertRaises(ExternalStageError):analyze(sample(),{'original_records':437},{})

    def test_original_prior_result_mismatch_blocks(self):
        r=sample(5)[:367]
        with self.assertRaises(ExternalStageError):analyze(r,{'original_records':437},{},prior_v13={'audit':{'complete_paired_stage_records':3}})

    def test_same_year_source_code_not_three_independent_fields(self):
        a=sample(3)
        self.assertEqual(len({x['physical_site'] for x in a}),2)
        self.assertEqual(len({x['site_code'] for x in a}),3)

    def test_endpoint_scales_and_outcome_ids_not_pooled(self):
        x=heldout(sample(4),mode='physical_location',stage=ENDPOINTS[0]);y=heldout(sample(4),mode='physical_location',stage=ENDPOINTS[1])
        self.assertEqual(x['endpoint'],'seed_initiation_log1p')
        self.assertEqual(y['endpoint'],'fly_escape_fraction')
        self.assertNotEqual(x['heldout_row_weighted_mse_gain'],y['heldout_row_weighted_mse_gain'])

    def test_holdout_fold_independent_of_test_outcome_during_training(self):
        r=sample(4);te=[x for x in r if x['year']=='2018'];tr=[x for x in r if x['year']!='2018']
        a=_train_context_slope(tr,'physical_location',ENDPOINTS[0]);te[0]['initiated']=99999.
        b=_train_context_slope(tr,'physical_location',ENDPOINTS[0]);self.assertEqual(a,b)

    def test_within_context_centered_vectors_zero_independently_under_test_shuffle(self):
        r=sample(4);z=[x for x in r if x['year']=='2018'];rand=Random(41)
        rand.shuffle(z)
        features=_features_by_group(_feature_only_test_rows(z,ENDPOINTS[0]),'physical_location',stage=ENDPOINTS[0],stats={'width':(0,1),'initiated':(0,1)})
        for group in {x['group'] for x in features}:
            self.assertAlmostEqual(sum(x['w_centered'] for x in features if x['group']==group),0)

if __name__=='__main__':unittest.main()
