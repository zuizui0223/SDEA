"""SDEA v12: synthetic source perturbations, no actual author raw data bundled."""
import json
import math
import unittest
from pathlib import Path
from sdea.site_transfer import (
    SiteTransferError, SITES, ARMS, MEASURES, EXPECTED_BLOBS,
    git_blob_sha, pinned_source_url, verified_input,
    parse_original_seeds, aggregate_complete_blocks, _avg,
    paired_contrast, heldout_site_evaluation, joint_sign_counts,
    validate_prior_independent_receipt, run_original,
)


def synthetic_csv():
    head=',Site_Code,Block_Code,Treatment_Code,Vial_Code,Viable_Seeds,NViable_Seeds,Damaged_Seeds,Missing_Seeds,Insects,Flowers'
    rows=[]
    for s in SITES:
        for block in ('B1','B2'):
            for arm in ARMS:
                viable={'Control':8,'Pollinator Exclusion':2,'Predator Exclusion':8,'Predator Model':8}[arm]
                ins={'Control':3,'Pollinator Exclusion':1,'Predator Exclusion':2,'Predator Model':2}[arm]
                rows.append(f'{len(rows)+1},{s},{block},{arm},1,{viable},0,0,{10-viable},{ins},10')
    return ('\n'.join([head,*rows])+'\n').encode()


def synthetic_data():
    rows,summary=parse_original_seeds(synthetic_csv(),strict=False)
    blocks,counters=aggregate_complete_blocks(rows,strict=False)
    return blocks


class TestV12SiteTransfer(unittest.TestCase):
    def test_git_blob_digest_exact_git_spec(self):
        import hashlib
        self.assertEqual(git_blob_sha(b'hello'),hashlib.sha1(b'blob 5\0hello').hexdigest())
    def test_pinned_author_original_url(self):
        entry={'repo':'Yogita25/Multitrophic_Interactions','commit':'8a151c62ed069f775a16b9bc5fbdce209cbb8691','path':'seeds.csv'}
        self.assertEqual(pinned_source_url(entry),'https://raw.githubusercontent.com/Yogita25/Multitrophic_Interactions/8a151c62ed069f775a16b9bc5fbdce209cbb8691/seeds.csv')
    def test_untrusted_repo_fails_closed(self):
        with self.assertRaises(SiteTransferError):pinned_source_url({'repo':'unknown','commit':'f'*40,'path':'seeds.csv'})
    def test_missing_source_requires_explicit_fetch(self):
        with self.assertRaises(SiteTransferError):verified_input({'blob_sha':'wrong'},fetch_pinned=False)
    def test_sha_mismatch_fails_even_for_valid_csv(self):
        with self.assertRaises(SiteTransferError):verified_input({'blob_sha':'0'*40},raw=synthetic_csv())
    def test_expected_github_blob_sha_keys(self):
        self.assertEqual(set(EXPECTED_BLOBS),{'original_author_seeds','aza3_independent_receipt'})
    def test_original_source_strict_615_gate(self):
        with self.assertRaises(SiteTransferError):parse_original_seeds(synthetic_csv(),strict=True)
    def test_synthetic_valid_csv_kept(self):
        rows,c=parse_original_seeds(synthetic_csv(),strict=False)
        self.assertEqual(len(rows),32)
        self.assertEqual(c['raw_head_rows'],32)
        self.assertEqual(c['eligible_head_rows'],32)
        self.assertEqual(c['ineligible_head_rows'],0)
    def test_synthetic_complete_blocks_and_proxy_grain(self):
        rows,_=parse_original_seeds(synthetic_csv(),strict=False)
        blocks,c=aggregate_complete_blocks(rows,strict=False)
        self.assertEqual(len(blocks),8)
        self.assertEqual(c['inferred_site_block_arm_plant_proxies'],32)
        self.assertEqual(c['source_site_blocks'],8)
        self.assertEqual(c['complete_site_blocks'],{s:2 for s in SITES})
    def test_incomplete_four_arm_block_is_not_counted(self):
        rows,_=parse_original_seeds(synthetic_csv(),strict=False)
        rows=[r for r in rows if not(r['Site_Code']=='campus' and r['Block_Code']=='B2' and r['Treatment_Code']=='Control')]
        blocks,c=aggregate_complete_blocks(rows,strict=False)
        self.assertEqual(len(blocks),7)
    def test_missing_numeric_rejected(self):
        raw=synthetic_csv().replace(b',8,0,0,2,3,10',b',BAD,0,0,2,3,10',1)
        with self.assertRaises(SiteTransferError):parse_original_seeds(raw,strict=False)
    def test_negative_non_missing_rejected(self):
        raw=synthetic_csv().replace(b',8,0,0,2,3,10',b',-8,0,0,18,3,10',1)
        with self.assertRaises(SiteTransferError):parse_original_seeds(raw,strict=False)
    def test_signed_missing_negative_is_admitted(self):
        raw=synthetic_csv().replace(b',8,0,0,2,3,10',b',11,0,0,-1,3,10',1)
        rows,c=parse_original_seeds(raw,strict=False)
        self.assertEqual(len(rows),32)
        self.assertEqual(c['negative_signed_missing_raw_heads'],1)
    def test_mass_balance_error_fails_closed(self):
        raw=synthetic_csv().replace(b',8,0,0,2,3,10',b',8,0,0,3,3,10',1)
        with self.assertRaises(SiteTransferError):parse_original_seeds(raw,strict=False)
    def test_zero_exposure_rejected(self):
        raw=synthetic_csv().replace(b',8,0,0,2,3,10',b',0,0,0,0,3,0',1)
        with self.assertRaises(SiteTransferError):parse_original_seeds(raw,strict=False)
    def test_invalid_site_and_treatment_excluded(self):
        raw=synthetic_csv()+b'100,other,B3,Control,1,8,0,0,2,3,10\n'
        rows,c=parse_original_seeds(raw,strict=False)
        self.assertEqual(len(rows),32)
        self.assertEqual(c['ineligible_head_rows'],1)
    def test_paired_seed_and_insect_differences(self):
        b=synthetic_data()[0]
        self.assertAlmostEqual(paired_contrast(b,intervention='Pollinator Exclusion',metric='viable_fraction'),-.6)
        self.assertAlmostEqual(paired_contrast(b,intervention='Pollinator Exclusion',metric='emerged_insects_per_head'),-2)
    def test_unsupported_endpoint_or_intervention_rejected(self):
        b=synthetic_data()[0]
        with self.assertRaises(SiteTransferError):paired_contrast(b,intervention='Control',metric='viable_fraction')
        with self.assertRaises(SiteTransferError):paired_contrast(b,intervention='Pollinator Exclusion',metric='ancestral_regain')
    def test_leave_site_out_reference_oracle_constant_effect(self):
        b=synthetic_data()
        out=heldout_site_evaluation(b,intervention='Pollinator Exclusion',metric='viable_fraction')
        self.assertAlmostEqual(out['site_mean_mse_zero_reference'],.36)
        self.assertAlmostEqual(out['site_mean_mse_train_other_sites'],0)
        self.assertEqual(out['site_benefit_count'],4)
        self.assertEqual(out['site_sign_match_count'],4)
        self.assertTrue(all(x['training_blocks']==6 and x['heldout_blocks']==2 for x in out['folds']))
    def test_target_site_perturbation_cannot_leak_into_own_prediction(self):
        import copy
        b=synthetic_data()
        base=heldout_site_evaluation(b,intervention='Pollinator Exclusion',metric='viable_fraction')
        changed=copy.deepcopy(b)
        for x in changed:
            if x['site']=='ditch':x['arms']['Pollinator Exclusion']['viable_fraction']=.8
        post=heldout_site_evaluation(changed,intervention='Pollinator Exclusion',metric='viable_fraction')
        self.assertEqual(next(x['predicted_from_other_sites'] for x in base['folds'] if x['heldout_site']=='ditch'),next(x['predicted_from_other_sites'] for x in post['folds'] if x['heldout_site']=='ditch'))
        self.assertNotEqual(next(x['observed_site_mean'] for x in base['folds'] if x['heldout_site']=='ditch'),next(x['observed_site_mean'] for x in post['folds'] if x['heldout_site']=='ditch'))
    def test_site_units_not_counted_as_independent_plant_experiments(self):
        b=synthetic_data()
        out=heldout_site_evaluation(b,intervention='Pollinator Exclusion',metric='viable_fraction')
        self.assertEqual(len(out['folds']),4)
        self.assertEqual(sum(x['heldout_blocks'] for x in out['folds']),8)
        self.assertEqual(out['holdout_unit'],'one whole original study site')
    def test_joint_signs_descriptive_no_causal_interpretation(self):
        signs=joint_sign_counts(synthetic_data(),'Pollinator Exclusion')
        self.assertEqual(signs,{'viable_down__emerged_down':8})
    def test_4_sites_required_not_a_single_site_pseudoreplication(self):
        b=[r for r in synthetic_data() if r['site']=='campus']
        with self.assertRaises(SiteTransferError):heldout_site_evaluation(b,intervention='Pollinator Exclusion',metric='viable_fraction')
    def test_previous_receipt_version_mismatch_fails(self):
        with self.assertRaises(SiteTransferError):validate_prior_independent_receipt({'version':'invented'}, {})
    def test_incomplete_source_manifest_fails(self):
        with self.assertRaises(SiteTransferError):run_original(b'x',b'{}',{'sources':{}})
    def test_site_mean_vs_block_mse_under_unbalanced_counts(self):
        b=synthetic_data()
        out=heldout_site_evaluation(b,intervention='Pollinator Exclusion',metric='viable_fraction')
        self.assertAlmostEqual(out['block_mse_zero_secondary'],.36)
        self.assertAlmostEqual(out['block_mse_train_other_sites_secondary'],0)
    def test_strict_eligible_source_proxy_control(self):
        rows,_=parse_original_seeds(synthetic_csv(),strict=False)
        with self.assertRaises(SiteTransferError):aggregate_complete_blocks(rows,strict=True)
    def test_approved_source_lock_expected_digests(self):
        manifest=json.loads((Path(__file__).resolve().parents[1]/'data/v12_source_lock.json').read_text())
        self.assertEqual({k:v['blob_sha'] for k,v in manifest['sources'].items()},EXPECTED_BLOBS)


if __name__=='__main__':unittest.main()
