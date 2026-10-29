"""SDEA v15 offline tests: source identity, join, peer-context, no-leakage."""
from __future__ import annotations
import csv
import hashlib
import io
import math
import random
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sdea import ulex_partial as m


def mkcsv(header,rows):
    out=io.StringIO(newline='');w=csv.writer(out,lineterminator='\n');w.writerow(header)
    w.writerows(rows);return out.getvalue().encode()


def corpus(n=210, *, genotypes=225):
    sites=[f'S{i}' for i in range(6)]
    fruit=[];traits=[];poll=[]
    for i in range(n):
        site=sites[i%6];ident=f'P{i:03}'
        scars=30+(i%9);fruits=round(scars*(.2+.025*(i%5)))
        area=10.+i%17;weight=2.+(i%11)/3.
        fruit.append([site,ident,450+i%13,weight,area,scars,fruits,scars-fruits]);traits.append([site,ident,1,weight,area,450+i%13]);
    poll.append([2014,'xyz',9,2,1,0,1,.22,1])
    labels=[f'P{i:03}' for i in range(genotypes)]
    kin=[]
    for i in range(genotypes):
        kin.append([labels[i]]+[1.0 if i==j else (.08 if i%6==j%6 else -.012) for j in range(genotypes)])
    return {
      'fruits.and.mean.floral.traits.csv':mkcsv(m.EXPECTED_FRUIT_HEADER,fruit),
      'floral.traits.csv':mkcsv(m.EXPECTED_TRAIT_HEADER,traits),
      'pollinator.census.csv':mkcsv(m.EXPECTED_POLLINATOR_HEADER,poll),
      'relatedness.matrix.csv':mkcsv(('ind',*labels),kin),
      'README.md':b'original README synthetic testing only',
    }


def rows_for_model(n=60):
    out=[]
    for i in range(n):
        area=2+i/10;weight=1+(i%7)/3; elev=200+i/4
        out.append({'site':f'S{i%6}','ind':f'P{i:03}','elev':elev,'area':area,'weight':weight,
                    'scars':100,'fruits':round(100*(.1+.002*i)),'fruit_set':.1+.002*i})
    return out

class SourceAndGrain(unittest.TestCase):
    def test_registered_doi_and_md5_shapes(self):
        self.assertEqual(m.RECORD_ID,7761289)
        self.assertEqual(len(m.EXPECTED_MD5),5)
        self.assertTrue(all(len(x)==32 for x in m.EXPECTED_MD5.values()))

    def test_source_gate_requires_exact_five_files(self):
        c=corpus()
        with self.assertRaisesRegex(m.UlexGateError,'MISSING_OR_EXTRA'):
            m.verify_payloads({k:v for k,v in c.items() if k!='README.md'})

    def test_modified_source_bytes_fail_closed(self):
        c=corpus()
        with self.assertRaisesRegex(m.UlexGateError,'MD5_MISMATCH'):
            m.verify_payloads(c)

    def test_url_disallows_untrusted_domains_and_credentials(self):
        self.assertTrue(m.valid_url('https://zenodo.org/api/records/7761289'))
        for s in ['http://zenodo.org','https://example.org/a','https://zenodo.org.evil.com','https://abc@zenodo.org/a']:
            self.assertFalse(m.valid_url(s))

    def test_source_header_matches_probe_v15_contract(self):
        self.assertEqual(m.EXPECTED_FRUIT_HEADER[-3:],('scars','fruits','nofruit'))
        self.assertNotIn('ind',m.EXPECTED_POLLINATOR_HEADER)

    def test_parse_synthetic_cohort_six_sites(self):
        data=corpus()
        r,a=m.parse_fruit_and_traits(data)
        self.assertEqual(len(r),210)
        self.assertEqual(len(a['eligible_by_site']),6)
        self.assertEqual(sum(a['eligible_by_site'].values()),210)
        self.assertFalse(a['pollinator_plant_id_in_original_header'])
        self.assertEqual(a['fruit_records_with_any_repeated_flower_key'],210)

    def test_fruit_set_is_fruits_over_original_scars(self):
        rows,_=m.parse_fruit_and_traits(corpus())
        self.assertAlmostEqual(rows[0]['fruit_set'],rows[0]['fruits']/rows[0]['scars'])

    def test_duplicate_site_individual_is_blocked(self):
        d=corpus();d['fruits.and.mean.floral.traits.csv']+=d['fruits.and.mean.floral.traits.csv'].splitlines()[1]+b'\n'
        with self.assertRaisesRegex(m.UlexGateError,'DUPLICATE_SITE_INDIVIDUAL'):
            m.parse_fruit_and_traits(d)

    def test_fruit_scars_arithmetic_fails_closed(self):
        d=corpus(); s=d['fruits.and.mean.floral.traits.csv'].decode();s=s.replace('S0,P000,450,2.0,10.0,30,6,24','S0,P000,450,2.0,10.0,30,6,20',1)
        self.assertNotEqual(s,d['fruits.and.mean.floral.traits.csv'].decode())
        d['fruits.and.mean.floral.traits.csv']=s.encode()
        with self.assertRaisesRegex(m.UlexGateError,'ACCOUNTING_INCONSISTENCY'):
            m.parse_fruit_and_traits(d)

    def test_original_spreadsheet_error_counted_missing_not_zero(self):
        d=corpus()
        raw=d['fruits.and.mean.floral.traits.csv'].decode()
        self.assertIn('S0,P000,450,2.0,10.0,30,6,24',raw)
        d['fruits.and.mean.floral.traits.csv']=raw.replace(
            'S0,P000,450,2.0,10.0,30,6,24',
            'S0,P000,450,#VALUE!,10.0,30,6,24',1).encode()
        rows,a=m.parse_fruit_and_traits(d)
        self.assertEqual(len(rows),209)
        self.assertEqual(a['ineligible_reasons']['source_formula_error_weight'],1)
        self.assertEqual(a['ineligible_reasons']['missing_fruit_or_trait_field'],1)

    def test_arbitrary_nonnumeric_cannot_sneak_through_as_missing(self):
        d=corpus()
        raw=d['fruits.and.mean.floral.traits.csv'].decode()
        d['fruits.and.mean.floral.traits.csv']=raw.replace(
            'S0,P000,450,2.0,10.0,30,6,24',
            'S0,P000,450,UNKNOWN,10.0,30,6,24',1).encode()
        with self.assertRaisesRegex(m.UlexGateError,'NONNUMERIC'):
            m.parse_fruit_and_traits(d)

    def test_noninteger_fruit_bad(self):
        d=corpus();s=d['fruits.and.mean.floral.traits.csv'].decode().replace(',30,6,24',',30,6.1,24',1)
        d['fruits.and.mean.floral.traits.csv']=s.encode()
        with self.assertRaises(m.UlexGateError):m.parse_fruit_and_traits(d)

    def test_negative_trait_bad(self):
        d=corpus();s=d['fruits.and.mean.floral.traits.csv'].decode().replace('S0,P000,450,2.0','S0,P000,450,-2.0',1)
        d['fruits.and.mean.floral.traits.csv']=s.encode()
        with self.assertRaises(m.UlexGateError):m.parse_fruit_and_traits(d)

    def test_pollinator_schema_with_individual_id_must_stop(self):
        d=corpus();p=d['pollinator.census.csv'].decode().replace('year,locality,','ind,year,locality,',1)
        d['pollinator.census.csv']=p.encode()
        with self.assertRaises(m.UlexGateError):m.parse_fruit_and_traits(d)

    def test_two_code_plant_id_in_two_sites_cannot_join(self):
        data=corpus();s=data['fruits.and.mean.floral.traits.csv'].decode().replace('S1,P001','S1,P000',1)
        data['fruits.and.mean.floral.traits.csv']=s.encode()
        with patch.object(m,'verify_payloads',lambda _:None):
            with self.assertRaisesRegex(m.UlexGateError,'REUSED_GENOTYPE_KEY'):
                m.analyze(data)

    def test_no_genotype_match_fails_closed(self):
        d=corpus(210,genotypes=225)
        text=d['relatedness.matrix.csv'].decode()
        for i in range(225):text=text.replace(f'P{i:03}',f'Q{i:03}')
        d['relatedness.matrix.csv']=text.encode()
        with patch.object(m,'verify_payloads',lambda _:None):
            with self.assertRaisesRegex(m.UlexGateError,'SAME_PLANT'):
                m.analyze(d)

class Kinship(unittest.TestCase):
    def test_225_genomic_labels_and_negative_genomic_relatedness(self):
        K,a=m.parse_relatedness(corpus())
        self.assertEqual(len(K),225)
        self.assertLess(K['P000']['P001'],0)
        self.assertGreater(K['P000']['P006'],0)
        self.assertEqual(a['relatedness_observed_cells'],225*225)

    def test_changed_genomic_header_fails(self):
        d=corpus();d['relatedness.matrix.csv']=d['relatedness.matrix.csv'].replace(b'ind,P000',b'whatever,P000',1)
        with self.assertRaisesRegex(m.UlexGateError,'EXPECTED_225'):
            m.parse_relatedness(d)

    def test_genetic_relatedness_asymmetric_fails(self):
        d=corpus();s=d['relatedness.matrix.csv'].decode().splitlines();
        # replace the first row's second non-diagonal value with an asymmetric entry
        a=s[1].split(',');a[2]='0.8';s[1]=','.join(a);d['relatedness.matrix.csv']=('\n'.join(s)+'\n').encode()
        with self.assertRaisesRegex(m.UlexGateError,'ASYMMETRIC'):
            m.parse_relatedness(d)

    def test_no_kinship_does_not_become_zero_effect(self):
        d=corpus();cols=d['relatedness.matrix.csv'].decode().splitlines();d['relatedness.matrix.csv']=('\n'.join(cols[:-1])+'\n').encode()
        with self.assertRaises(m.UlexGateError):m.parse_relatedness(d)

class ModelControls(unittest.TestCase):
    def test_ridge_solution_matches_scalar_known(self):
        self.assertEqual(m._solve([[2.0]],[4.0]),[2.0])
        a=m._solve([[2.,0.],[0.,3.]],[4.,6.])
        self.assertEqual(a,[2.,2.])

    def test_predict_restricted_to_probability_scale(self):
        r=rows_for_model()
        f=m.fit(r,traits=True)
        self.assertTrue(all(0<=m.predict(f,x)<=1 for x in r))

    def test_fit_is_insensitive_to_heldout_response(self):
        r=rows_for_model();tr=r[:-6];te=r[-1].copy()
        f=m.fit(tr,traits=True)
        p=m.predict(f,te)
        te['fruit_set']=1000
        self.assertEqual(p,m.predict(f,te))

    def test_relatedness_residual_test_target_not_used(self):
        r=rows_for_model();tr=r[:-6];te=r[-1].copy();f=m.fit(tr,traits=True)
        K={te['ind']:{x['ind']:.1 for x in tr}}
        v=m.related_residual_prediction(te,tr,f,K)
        te['fruit_set']=-4
        self.assertEqual(v,m.related_residual_prediction(te,tr,f,K))

    def test_genetic_shrinkage_infinite_gives_traits_only(self):
        r=rows_for_model();tr=r[:-6];te=r[-1];f=m.fit(tr,traits=True)
        K={te['ind']:{x['ind']:.5 for x in tr}}
        v=m.related_residual_prediction(te,tr,f,K,shrinkage=1e12)
        self.assertAlmostEqual(v,m.predict(f,te),places=7)

    def test_nonpositive_shrinkage_fails(self):
        r=rows_for_model();f=m.fit(r[:-6],traits=True)
        with self.assertRaises(m.UlexGateError):m.related_residual_prediction(r[-1],r[:-6],f,{},shrinkage=0)

    def test_source_gating_malformed_site_coverage(self):
        d=corpus(n=60)
        with self.assertRaises(m.UlexGateError):m.parse_fruit_and_traits(d)

    def test_full_synthetic_analysis_includes_six_heldout_sites(self):
        d=corpus()
        with patch.object(m,'verify_payloads',lambda _:None):
            a=m.analyze(d)
        self.assertEqual(a['audit']['matched_genomic_flower_fruit_individuals'],210)
        self.assertEqual(len(a['heldout_sites']),6)
        self.assertEqual(sum(x['heldout_individuals'] for x in a['heldout_sites']),210)
        self.assertFalse(a['full_ANTO_EGWE_izu_model_admitted'])

    def test_no_pseudo_independent_pollinator_genomic_match(self):
        d=corpus()
        with patch.object(m,'verify_payloads',lambda _:None):
            a=m.analyze(d)
        self.assertFalse(a['audit']['pollinator_census_individual_link'])
        self.assertEqual(a['audit']['pollinator_census_rows_quarantined_unjoined'],1)

    def test_holdout_site_not_used_in_training(self):
        d=corpus()
        with patch.object(m,'verify_payloads',lambda _:None):a=m.analyze(d)
        n=a['audit']['matched_genomic_flower_fruit_individuals']
        self.assertTrue(all(x['heldout_individuals']+x['train_individuals']==n for x in a['heldout_sites']))

    def test_report_has_two_distinct_nested_increments(self):
        d=corpus()
        with patch.object(m,'verify_payloads',lambda _:None):a=m.analyze(d)
        self.assertEqual(a['trait_increment']['site_count'],6)
        self.assertEqual(a['genetic_increment_conditional_on_traits']['site_count'],6)
        self.assertNotIn('significant',a)

if __name__=='__main__':unittest.main()
