"""Structural, original-source-gate and inference-boundary v9 tests (offline)."""
import csv
import hashlib
import io
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sdea.empirical_admission import (SourceGateError, blob_sha, pin_url,
  validate_source_bytes, load_inputs, parse_maternal_joint, summarize_source_receipts,
  run)


class TestSDEAEmpiricalAdmission(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = json.loads((ROOT/'data/v9_source_lock.json').read_text())

    @staticmethod
    def csv_bytes(rows):
        buf=io.StringIO(newline='')
        cols=['year','plt','trt','plt.gsd','ap.hd','o.hds']
        w=csv.DictWriter(buf,fieldnames=cols,lineterminator='\n')
        w.writeheader()
        for row in rows:
            w.writerow(dict(zip(cols,row)))
        return buf.getvalue().encode()

    def test_manifest_has_expected_six_unique_pins(self):
        self.assertEqual(len(self.manifest['sources']),6)
        for key,x in self.manifest['sources'].items():
            self.assertEqual(len(x['commit']),40)
            self.assertEqual(len(x['blob_sha']),40)
            self.assertTrue(all(ch in '0123456789abcdef' for ch in x['blob_sha']))
            self.assertEqual(x['repo'].count('/'),0)
            self.assertTrue(pin_url(x).startswith('https://raw.githubusercontent.com/zuizui0223/'+x['repo']+'/'+x['commit']+'/'))

    def test_git_blob_hash_known_reference(self):
        x=b'hello\n'
        self.assertEqual(blob_sha(x),hashlib.sha1(b'blob 6\0hello\n').hexdigest())
        self.assertEqual(blob_sha(b''),hashlib.sha1(b'blob 0\0').hexdigest())

    def test_blob_pinning_rejects_source_drift(self):
        fake={'path':'fixed.csv','blob_sha':blob_sha(b'expected')}
        self.assertEqual(validate_source_bytes(fake,b'expected'),b'expected')
        with self.assertRaisesRegex(SourceGateError,'SHA mismatch'):
            validate_source_bytes(fake,b'changed')

    def test_default_does_not_silently_access_network(self):
        with self.assertRaisesRegex(SourceGateError,'original local source root'):
            load_inputs(self.manifest)

    def test_grouped_same_mother_four_states(self):
        rows=[('2008','1','di','7','7','0'),('2008','2','di','2','0','2'),
              ('2008','3','di','5','2','3'),('2008','4','di','0','0','0')]
        result=parse_maternal_joint(self.csv_bytes(rows),check_original=False)
        v=result['year_treatment_strata']['2008|di']
        self.assertEqual((v['neither'],v['apical_only'],v['lateral_only'],v['both']),(1,1,1,1))
        self.assertEqual(v['either_success'],.75)
        self.assertEqual(v['apical_success'],.5)
        self.assertEqual(v['lateral_success'],.5)
        self.assertEqual(v['observed_two_route_added_coverage'],.25)
        self.assertEqual(v['observed_minus_independent'],0.)

    def test_independence_discrepancy_is_negative_covariance(self):
        cases=[('2008','1','di','1','1','0'),('2008','2','di','1','0','1'),
               ('2008','3','di','0','0','0')]
        x=parse_maternal_joint(self.csv_bytes(cases),check_original=False)['year_treatment_strata']['2008|di']
        self.assertAlmostEqual(x['same_mother_covariance'],-1/9)
        self.assertAlmostEqual(x['observed_minus_independent'],1/9)

    def test_false_mother_total_is_excluded_and_reported(self):
        cases=[('2007','327','dni','20','0','1'),('2007','328','dni','1','0','1')]
        result=parse_maternal_joint(self.csv_bytes(cases),check_original=False)
        self.assertEqual(result['source_mothers'],2)
        self.assertEqual(result['admitted_coherent_mothers'],1)
        self.assertEqual(result['excluded_incoherent'][0]['plant_id'],'327')

    def test_no_unexpected_year_treatment_or_repeated_mother(self):
        with self.assertRaises(SourceGateError):
            parse_maternal_joint(self.csv_bytes([('2008','3','bad','1','1','0')]),check_original=False)
        with self.assertRaisesRegex(SourceGateError,'duplicate original mother'):
            parse_maternal_joint(self.csv_bytes([('2008','3','di','1','1','0'),('2008','3','di','1','1','0')]),check_original=False)

    def test_no_invalid_or_missing_seed_counts(self):
        for invalid in ['NA','-3','1.25','nan','inf']:
            raw=self.csv_bytes([('2008','4','di',invalid,'0','0')])
            with self.subTest(invalid=invalid),self.assertRaises((SourceGateError,ValueError)):
                parse_maternal_joint(raw,check_original=False)

    def test_year_is_not_mother_identity(self):
        rows=[('2007','1','ndni','1','1','0'),('2008','1','ndni','1','0','1')]
        x=parse_maternal_joint(self.csv_bytes(rows),check_original=False)
        self.assertEqual(x['admitted_coherent_mothers'],2)

    def test_original_mode_must_see_106_rows(self):
        with self.assertRaisesRegex(SourceGateError,'expected 106'):
            parse_maternal_joint(self.csv_bytes([('2008','4','di','1','1','0')]),check_original=True)

    def test_no_real_files_altered_in_new_repository(self):
        # Original data remain in their separately versioned empirical repository.
        sources=self.manifest['sources']
        self.assertEqual(sources['mother_csv']['repo'],'aza3')
        self.assertIn('mother_plants.csv',sources['mother_csv']['path'])
        self.assertNotEqual(sources['mother_csv']['repo'],'SDEA')

    def test_no_model_only_source_counted_as_natural_study(self):
        self.assertEqual(len(self.manifest['admission_units']),5)
        self.assertNotIn('izu-core',str(self.manifest['admission_units']))
        self.assertEqual(len(self.manifest['full_contract']),5)

    def test_all_layers_required_without_imputation(self):
        assert 'same_organism_validated_joint_heritable_pathway' in self.manifest['full_contract']
        assert 'prospective_descendant_expression' in self.manifest['full_contract']
        assert 'trait_specific_demographic_establishment' in self.manifest['full_contract']

    def test_empty_or_extra_payloads_fail_closed(self):
        with self.assertRaisesRegex(SourceGateError,'payload set differ'):
            run(self.manifest,{})

    def test_no_real_inference_from_simulation_fields(self):
        self.assertIn('retrospective',self.manifest['status'])
        self.assertTrue(all(x['commit']!= 'main' for x in self.manifest['sources'].values()))

    def test_real_source_summary_receipt_shape(self):
        # This validates schema of the curated fixed receipts in a completely
        # offline synthetic-only unit test. Live pinned bytes are checked in CI.
        fake={
            'ryukyu_readiness':{'source_summary':{'public_run_count':25,'floral_RNA_runs_in_this_source':0,'phaseA_collected_here':0,'young_leaf_runs':25,'phaseA_target_individuals':298},
              'per_taxon':[{'same_plant_verified_functional_toolkit':'NOT_MEASURED_IN_REGISTERED_SOURCE','heritable_repeated_reexpression':'NOT_MEASURED_IN_REGISTERED_SOURCE'} for _ in range(12)]},
            'izu_summary':{'analysis_lock':{'row_count':572,'site_count':8},'aggregate_metrics':{'C1':{'row_weighted_mse':1.0877352926660404},'C2':{'row_weighted_mse':1.1320924612331091}}},
            'oenothera_summary':{'schema':{'n_seed_family_rows':60,'n_maternal_plants':23},'models':{'M0':{'loo_mse':0.11618622977542375},'M1':{'loo_mse':0.09186971001970597}}},
            'anto_registry':{'dataset_identity':'non_systematic_source-verified_mechanism_examples_NOT_random_species_sample','cases':[
              {'clade':'Petunia','case_role':'MOLECULARLY_SUPPORTED_PIGMENT_REGAIN','trait_specific_viable_seed_fitness':'NOT_QUALIFIED'},
              {'clade':'Petunia','case_role':'MOLECULARLY_SUPPORTED_ALTERNATIVE_HUE_REGAIN','trait_specific_viable_seed_fitness':'NOT_QUALIFIED'},
              *[{'clade':'OTHER','case_role':'OTHER','trait_specific_viable_seed_fitness':'NOT_QUALIFIED'} for _ in range(4)]]},
            'izu_core_boundary':'This is a finite diploid simulation study, not a direct historical estimate.'.encode()
        }
        result=summarize_source_receipts({k:json.dumps(v).encode() if k!='izu_core_boundary' else v for k,v in fake.items()})
        self.assertEqual(result['ryukyu_leaf_RNA']['public_flower_stage_RNA_runs'],0)
        self.assertAlmostEqual(result['honshu_izu_source_reported']['adding_distance_mse_change_percent'],4.077937791128323)
        self.assertAlmostEqual(result['oenothera_source_reported']['adding_isolation_mse_change_percent'],-20.928917138217813)
        self.assertEqual(result['anto_literature_source_reported']['independent_genera_of_these_two'],1)
        self.assertFalse(result['izu_core']['count_as_natural_replication'])

    def test_unexpected_positive_molecular_claim_does_not_pass(self):
        # Fail closed rather than declaring trait regain on unknown evidence.
        x={'source_summary':{'public_run_count':25,'floral_RNA_runs_in_this_source':0,'phaseA_collected_here':0},
          'per_taxon':[{'same_plant_verified_functional_toolkit':'CONFIRMED','heritable_repeated_reexpression':'NOT_MEASURED_IN_REGISTERED_SOURCE'} for _ in range(12)]}
        with self.assertRaises(SourceGateError):
            summarize_source_receipts({'ryukyu_readiness':json.dumps(x).encode(),
               'izu_summary':b'{}','oenothera_summary':b'{}','anto_registry':b'{}','izu_core_boundary':b''})


if __name__ == '__main__':
    unittest.main()
