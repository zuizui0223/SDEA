"""SDEA v13: adverse, source/measurement and genuine year-heldout design tests."""
import csv
import io
import json
import math
import unittest
from copy import deepcopy
from unittest.mock import patch

from sdea.external_stage_transfer import (
    ORIGINAL_URL, ORIGINAL_MD5, RECEIPT_BLOB_SHA,
    ExternalStageError, HEADER, YEARS, ENDPOINTS,
    git_blob_sha, ismissing, read_original, validate_prior_receipt,
    _solve, _stats, fit, predict, heldout_year_score, analyze,
)


def rows_synthetic(n_per_year=10, *, width_gain=True, stage_gain=True):
    results=[]
    for k,year in enumerate(YEARS):
        for i in range(n_per_year):
            width=2. + i*.3 + (k%2)*.12
            begin=1. + .35*width + .01*k if width_gain else 2. + .01*k
            escape=.18+.08*width-.01*begin if stage_gain else .45+.01*k
            results.append({'year':year,'site_code':'agg' if i%2==0 else 'hyb',
                            'physical_site':'agg' if i%2==0 else 'upper',
                            'type':('A','H','T')[i%3],'width':width,
                            'initiated':begin,'escape':max(0.,min(1.,escape))})
    return results


def source_text(rows_per_year=2):
    buf=io.StringIO(newline='')
    writer=csv.writer(buf,lineterminator='\n')
    writer.writerow(HEADER)
    for k,year in enumerate(YEARS):
        for i in range(rows_per_year):
            writer.writerow(['agg' if i%2==0 else 'hyb',year,'A' if i%2==0 else 'H',
                             f'{2.+.2*i}', '.', '.', '.', '.', '.', '3','2','0.7'])
    return buf.getvalue().encode()


class V13Tests(unittest.TestCase):
    def test_public_source_original_record_locked(self):
        self.assertIn('zenodo.org/records/5083399',ORIGINAL_URL)
        self.assertEqual(ORIGINAL_MD5,'12ed3dbe6d269c5f808f9c6b230d468c')
        self.assertEqual(len(RECEIPT_BLOB_SHA),40)

    def test_git_blob_hash_is_git_object_not_file_content_sha(self):
        raw=b'ab\n'
        self.assertNotEqual(git_blob_sha(raw),__import__('hashlib').sha1(raw).hexdigest())
        self.assertEqual(len(git_blob_sha(raw)),40)

    def test_missing_source_code_dot_vs_numeric_zero(self):
        self.assertTrue(ismissing('.'))
        self.assertTrue(ismissing(''))
        self.assertTrue(ismissing('NA'))
        self.assertFalse(ismissing('0'))
        self.assertFalse(ismissing('0.0'))

    def test_synthetic_csv_all_years_read_without_original_pin(self):
        records,audit=read_original(source_text(),strict=False)
        self.assertEqual(len(records),18)
        self.assertEqual(audit['complete_paired_stage_records'],18)
        self.assertEqual(audit['missing_per_source_column']['color'],18)
        self.assertEqual(audit['original_columns'],12)

    def test_altered_source_md5_blocked(self):
        with self.assertRaises(ExternalStageError):read_original(source_text(),strict=True)

    def test_required_header_must_be_exact(self):
        x=source_text().replace(b'Site,Year',b'Site2,Year')
        with self.assertRaisesRegex(ExternalStageError,'HEADER'):
            read_original(x,strict=False)

    def test_wrong_site_type_or_year_blocked(self):
        for old,new in ((b'2014',b'1999'),(b'agg',b'bad'),(b',A,',b',Q,')):
            txt=source_text().replace(old,new,1)
            with self.assertRaises(ExternalStageError):read_original(txt,strict=False)

    def test_synthetic_source_missing_required_year_blocked(self):
        x=source_text().replace(b'2010',b'2011')
        with self.assertRaises(ExternalStageError):read_original(x,strict=False)

    def test_negative_seed_count_fail_closed(self):
        x=source_text().replace(b',3,2,',b',-3,2,',1)
        with self.assertRaisesRegex(ExternalStageError,'NEGATIVE'):
            read_original(x,strict=False)

    def test_positive_width_required(self):
        x=source_text().replace(b'2.0',b'0.0',1)
        with self.assertRaisesRegex(ExternalStageError,'NONPOSITIVE'):
            read_original(x,strict=False)

    def test_fly_escape_must_be_fraction(self):
        x=source_text().replace(b',0.7',b',1.2',1)
        with self.assertRaisesRegex(ExternalStageError,'INVALID_FLYESCAPE'):
            read_original(x,strict=False)

    def test_nonfinite_numeric_fails_closed(self):
        x=source_text().replace(b',0.7',b',nan',1)
        with self.assertRaisesRegex(ExternalStageError,'NONFINITE'):
            read_original(x,strict=False)

    def test_earlier_receipt_wrong_hash_fails_closed(self):
        with self.assertRaisesRegex(ExternalStageError,'GIT_BLOB_MISMATCH'):
            validate_prior_receipt(b'{}')

    def test_earlier_receipt_arithmetic_not_recycled_as_independent_study(self):
        obj={'source':{'metadata_author_md5':ORIGINAL_MD5},'audit':{'n_original_data_records':437,'n_original_fields':12,'seeds_greater_than_seedsinit':41},'analyses':{'width':{'summary':{'n_source_complete_records':367}}}}
        txt=json.dumps(obj).encode()
        with patch('sdea.external_stage_transfer.RECEIPT_BLOB_SHA',git_blob_sha(txt)):
            receipt=validate_prior_receipt(txt)
        self.assertEqual(receipt['complete_two_stage_n'],367)
        self.assertTrue(receipt['same_original_study_not_an_independent_replication'])

    def test_prior_audit_wrong_counts_blocked_even_if_hash_is_trusted(self):
        obj={'source':{'metadata_author_md5':ORIGINAL_MD5},'audit':{'n_original_data_records':436,'n_original_fields':12,'seeds_greater_than_seedsinit':41},'analyses':{'width':{'summary':{'n_source_complete_records':367}}}}
        txt=json.dumps(obj).encode()
        with patch('sdea.external_stage_transfer.RECEIPT_BLOB_SHA',git_blob_sha(txt)):
            with self.assertRaisesRegex(ExternalStageError,'PROFILE_DRIFT'):
                validate_prior_receipt(txt)

    def test_gauss_jordan_exact_two_variable(self):
        beta=_solve([[2.,0.],[0.,3.]],[4.,9.])
        self.assertEqual(beta,[2.,3.])

    def test_gauss_jordan_singular_fails_closed(self):
        with self.assertRaisesRegex(ExternalStageError,'SINGULAR'):
            _solve([[1.,1.],[2.,2.]],[3.,6.])

    def test_ridge_positive_and_invalid_penalty(self):
        for bad in (0,-1,float('nan')):
            with self.assertRaises(ExternalStageError):
                fit(rows_synthetic(),endpoint=ENDPOINTS[0],width=True,penalty=bad)

    def test_train_standardization_ignores_heldout_extreme(self):
        rows=rows_synthetic();train=[r for r in rows if r['year']!='2017']
        m=fit(train,endpoint=ENDPOINTS[0],width=True)
        changed=deepcopy(rows)
        for r in changed:
            if r['year']=='2017':r['width']=1e8
        other=fit([r for r in changed if r['year']!='2017'],endpoint=ENDPOINTS[0],width=True)
        self.assertEqual(m,other)

    def test_stage2_already_conditions_on_initiated(self):
        x=fit(rows_synthetic(),endpoint=ENDPOINTS[1],width=False)
        y=fit(rows_synthetic(),endpoint=ENDPOINTS[0],width=False)
        self.assertEqual(len(x['beta']),len(y['beta'])+1)

    def test_both_models_share_background_and_only_one_width_column(self):
        rows=rows_synthetic()
        for ep in ENDPOINTS:
            a=fit(rows,width=False,endpoint=ep)
            b=fit(rows,width=True,endpoint=ep)
            self.assertEqual(len(a['beta'])+1,len(b['beta']))
            self.assertEqual(a['stats'],b['stats'])

    def test_constant_width_has_zero_added_prediction_gain(self):
        rows=rows_synthetic()
        for r in rows:r['width']=3.
        for ep in ENDPOINTS:
            fold=heldout_year_score(rows,year='2017',endpoint=ep)
            self.assertAlmostEqual(fold['width_information_gain_mse'],0.,places=11)

    def test_real_structured_synthetic_stage1_gain(self):
        rows=rows_synthetic(width_gain=True)
        result=heldout_year_score(rows,year='2017',endpoint=ENDPOINTS[0])
        self.assertGreater(result['width_information_gain_mse'],0.)

    def test_null_width_predictor_not_forced_positive(self):
        rows=rows_synthetic(width_gain=False,stage_gain=False)
        report=analyze(rows,{'synthetic':True},{'synthetic_receipt':True})
        self.assertTrue(all('heldout_weighted_mse_gain' in report['outcomes'][ep] for ep in ENDPOINTS))
        self.assertEqual(report['external_v12_transport']['decision'],
                         'NOT_IDENTIFIABLE_NO_COMMON_INTERVENTION_OR_ENDPOINT')

    def test_fraction_prediction_is_clipped_to_unit_interval(self):
        rows=rows_synthetic()
        m=fit(rows,width=True,endpoint=ENDPOINTS[1])
        extreme=deepcopy(rows[0]);extreme['width']=1000
        val=predict(m,extreme)
        self.assertGreaterEqual(val,0.)
        self.assertLessEqual(val,1.)

    def test_year_holdout_one_prediction_per_row(self):
        rows=rows_synthetic()
        out=analyze(rows,{'synthetic':True},{'synthetic_receipt':True})
        for ep in ENDPOINTS:
            d=out['outcomes'][ep]
            self.assertEqual(d['heldout_original_records_once_each'],len(rows))
            self.assertEqual(d['heldout_years'],9)
            self.assertTrue(all(x['train_rows']+x['heldout_rows']==len(rows) for x in d['folds']))

    def test_year_labels_cannot_be_computed_as_individual_IDs(self):
        x=analyze(rows_synthetic(),{'synthetic':True},{'synthetic_receipt':True})
        self.assertIn('individual ID absent',x['design']['unit'])
        self.assertTrue(x['design']['not_an_original_pre_registration'])

    def test_missing_year_in_heldout_fails_closed(self):
        x=[r for r in rows_synthetic() if r['year']!='2018']
        with self.assertRaisesRegex(ExternalStageError,'EMPTY_TEST_YEAR'):
            heldout_year_score(x,year='2018',endpoint=ENDPOINTS[0])

    def test_zero_width_gain_never_mislabels_cirsium_bagging_replication(self):
        x=analyze(rows_synthetic(),{'synthetic':True},{'synthetic_receipt':True})
        self.assertTrue(x['external_v12_transport']['no_pooled_cirsium_ipomopsis_treatment_effect'])
        self.assertNotIn('bagging',x['design']['candidate_added_variable'])


if __name__=='__main__':unittest.main()
