"""Adversarial checks for source identity, plant IDs, observation effort, LOSO, nulls."""
import csv
import hashlib
import io
import math
import random
import unittest
from copy import deepcopy
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sdea import sinapis_process as m

def fake_csv(cols,rows):
    s=io.StringIO(newline='')
    w=csv.writer(s,delimiter=';',lineterminator='\n')
    w.writerow(cols)
    w.writerows(rows)
    return s.getvalue().encode()

def corpus(*,nper=4,pure_background=False,none_visits=False):
    plants=[];floral=[];visits=[]
    for shelter in range(1,6):
        for treatment in m.TREATMENTS:
            for i in range(nper):
                plant=f'P{shelter}_{treatment}_{i}'
                family=f'Fam{i}'
                plot=f'{treatment}{shelter}'
                meta=[plant,family,treatment,str(shelter),plot]
                plants.append(meta)
                width=3.+i*.25+(shelter%2)*.3
                display=10.+i*.9+shelter*.2
                flowers=100.+i*14+shelter*7
                counts=((i+2)*2,(i%3+1)*2) if not none_visits else (0,0)
                visited=sum(counts)
                if pure_background:
                    seeds=(flowers*3)*(1.5 if treatment=='Control' else 1.)
                else:
                    seeds=flowers*(1.1+.15*width+visited*.023)
                r={key:'2.0' for key in m.FRUIT_FIELDS}
                r.update({'plantID':plant,'plotID':plot,'shelterID':str(shelter),
                          'motherID':family,'treatment':treatment,'sumflower':str(flowers),
                          'seedsperplant':str(seeds),'meandisplayflower':str(display),
                          'meanpetalwidth':str(width)})
                floral.append([r[k] for k in m.FRUIT_FIELDS])
                for j,n in enumerate(counts):
                    visits.append([plant,('Apis','Andrena')[j],family,treatment,str(shelter),plot,
                                  '42','100','5',str(n),'1','0.01'])
    return {'Plantlist.csv':fake_csv(m.PLANT_FIELDS,plants),
            'Floral_traits.csv':fake_csv(m.FRUIT_FIELDS,floral),
            'Int_data.csv':fake_csv(m.VISIT_FIELDS,visits),
            'README.md':b'Synthetic README, never the original source'}

def mutate(src,filename,old,new):
    x=dict(src)
    assert old in src[filename]
    x[filename]=src[filename].replace(old,new,1)
    return x

class OriginalAndSource(unittest.TestCase):
    def test_manifest_names_and_digests_frozen(self):
        self.assertEqual(m.RECORD_ID,7335086)
        self.assertEqual(len(m.ORIGINAL_MD5),4)
        self.assertTrue(all(len(x)==32 for x in m.ORIGINAL_MD5.values()))
        self.assertEqual(m.N_ORIGINAL_PLANTS,53)
        self.assertEqual(m.N_ORIGINAL_VISIT_GROUP_ROWS,427)

    def test_only_original_zenodo_https_is_trusted(self):
        self.assertTrue(m._safe_url('https://zenodo.org/api/records/7335086/files/README.md/content'))
        for u in ['http://zenodo.org/a','https://zenodo.org.evil.org/a','https://user:pw@zenodo.org/a','file:///etc/passwd']:
            self.assertFalse(m._safe_url(u))

    def test_changed_or_missing_source_fails_closed(self):
        d=corpus()
        with self.assertRaisesRegex(m.SinapisGateError,'MD5'):
            m.validate_bytes(d)
        d.pop('README.md')
        with self.assertRaisesRegex(m.SinapisGateError,'FILE_SET'):
            m.validate_bytes(d)

    def test_source_original_row_lock_only_strict(self):
        x=corpus()
        rows,a=m.source_join(x,strict=False)
        self.assertEqual(len(rows),60)
        with self.assertRaisesRegex(m.SinapisGateError,'ORIGINAL_FILE_MD5'):
            m.source_join(x,strict=True)

    def test_unexpected_header_rejected(self):
        d=mutate(corpus(),'Plantlist.csv',b'plantID',b'plant id')
        with self.assertRaisesRegex(m.SinapisGateError,'HEADER|COLUMN'):
            m.source_join(d,strict=False)

    def test_duplicate_plant_id_rejected(self):
        d=corpus()
        row=d['Plantlist.csv'].decode().splitlines()[1]
        d['Plantlist.csv']+= (row+'\n').encode()
        with self.assertRaisesRegex(m.SinapisGateError,'DUPLICATE'):
            m.source_join(d,strict=False)

    def test_mother_identity_mismatch_rejected(self):
        d=mutate(corpus(),'Int_data.csv',b'Fam0',b'UNKNOWN')
        with self.assertRaisesRegex(m.SinapisGateError,'MISMATCH'):
            m.source_join(d,strict=False)

    def test_rainfall_treatment_mismatch_rejected(self):
        d=mutate(corpus(),'Floral_traits.csv',b';Control;',b';Drought;')
        with self.assertRaisesRegex(m.SinapisGateError,'MISMATCH'):
            m.source_join(d,strict=False)

    def test_plot_and_shelter_mismatch_rejected(self):
        d=mutate(corpus(),'Int_data.csv',b';1;Control1;',b';5;Control1;')
        with self.assertRaisesRegex(m.SinapisGateError,'MISMATCH'):
            m.source_join(d,strict=False)

    def test_unknown_nonnumeric_fails_closed(self):
        d=corpus()
        rs=list(csv.DictReader(io.StringIO(d['Floral_traits.csv'].decode()),delimiter=';'))
        rs[0]['seedsperplant']='not_a_seed_count'
        d['Floral_traits.csv']=fake_csv(m.FRUIT_FIELDS,[[r[k] for k in m.FRUIT_FIELDS] for r in rs])
        with self.assertRaisesRegex(m.SinapisGateError,'NONNUMERIC'):
            m.source_join(d,strict=False)

    def test_missing_numeric_not_zero_imputation(self):
        d=corpus()
        # Fixed source schema; replace a numeric only in a known field by row dict rewrite.
        l=list(csv.DictReader(io.StringIO(d['Floral_traits.csv'].decode()),delimiter=';'))
        l[0]['meanpetalwidth']='NA'
        d['Floral_traits.csv']=fake_csv(m.FRUIT_FIELDS,[[r[k] for k in m.FRUIT_FIELDS] for r in l])
        records,a=m.source_join(d,strict=False)
        self.assertEqual(len(records),59)
        self.assertEqual(a['ineligibility_counts']['missing_meanpetalwidth'],1)
        self.assertEqual(a['ineligibility_counts']['ineligible_complete_source_outcomes_and_trait'],1)

    def test_original_required_nonmissing_plant_key(self):
        d=mutate(corpus(),'Plantlist.csv',b'P1_Control_0',b'')
        with self.assertRaises(m.SinapisGateError):m.source_join(d,strict=False)

    def test_missing_visit_plant_not_no_interactions(self):
        d=corpus()
        l=[r for r in csv.DictReader(io.StringIO(d['Int_data.csv'].decode()),delimiter=';') if r['plantID']!='P1_Control_0']
        d['Int_data.csv']=fake_csv(m.VISIT_FIELDS,[[r[k] for k in m.VISIT_FIELDS] for r in l])
        _,a=m.source_join(d,strict=False)
        self.assertEqual(a['ineligibility_counts']['no_direct_visitation_record'],1)
        self.assertEqual(a['individuals_with_all_model_fields'],59)

    def test_visit_time_repeated_across_groups_not_added(self):
        d=corpus()
        r,a=m.source_join(d,strict=False)
        self.assertEqual(a['individuals_with_all_model_fields'],60)
        row=next(x for x in r if x['treatment']=='Control' and x['shelter']=='1')
        # First plant has group counts 4+2 and *single* observation time 100 min.
        self.assertAlmostEqual(row['log_visit_rate'],math.log1p(6/(100/60)))
        self.assertTrue(a['visitor_observation_time_reused_once_per_plant'])

    def test_inconsistent_visit_effort_blocks(self):
        d=corpus()
        text=d['Int_data.csv'].decode().splitlines()
        # Modify exactly one observation time in second visitor group of first plant.
        for i,line in enumerate(text):
            if 'P1_Control_0;Andrena;' in line:
                text[i]=line.replace(';100;',';125;');break
        d['Int_data.csv']=('\n'.join(text)+'\n').encode()
        with self.assertRaisesRegex(m.SinapisGateError,'EXPOSURE'):
            m.source_join(d,strict=False)

    def test_visitor_group_rows_are_not_plant_replicates(self):
        d=corpus()
        records,a=m.source_join(d,strict=False)
        self.assertEqual(len(records),60)
        self.assertEqual(a['visitor_group_source_rows'],120)
        self.assertEqual(a['plants_with_visit_groups'],60)

    def test_fruit_negative_seeds_rejected(self):
        d=corpus();rows=list(csv.DictReader(io.StringIO(d['Floral_traits.csv'].decode()),delimiter=';'))
        rows[0]['seedsperplant']='-4'
        d['Floral_traits.csv']=fake_csv(m.FRUIT_FIELDS,[[r[k] for k in m.FRUIT_FIELDS] for r in rows])
        with self.assertRaisesRegex(m.SinapisGateError,'INVALID'):
            m.source_join(d,strict=False)

    def test_plant_counts_per_shelter_preserved(self):
        rs,a=m.source_join(corpus(),strict=False)
        self.assertEqual(set(a['observed_shelter_counts']),set(m.SHELTERS))
        self.assertTrue(all(z==12 for z in a['observed_shelter_counts'].values()))

class ForecastAndControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows,cls.audit=m.source_join(corpus(),strict=False)

    def test_training_ridge_rejects_invalid_penalty(self):
        with self.assertRaises(m.SinapisGateError):m.train_model(self.rows,penalty=0)
        with self.assertRaises(m.SinapisGateError):m.train_model(self.rows,model='M9')

    def test_training_only_normalization(self):
        a=m.train_model([r for r in self.rows if r['shelter']!='1'],model='M2')
        te=[r for r in self.rows if r['shelter']=='1'][0]
        changed=dict(te,log_seed=1000000)
        self.assertEqual(m.predict(a,te),m.predict(a,changed))

    def test_model_feature_stage_nested(self):
        self.assertEqual(m.MODEL_FEATURES['M0'],('log_flower',))
        self.assertEqual(m.MODEL_FEATURES['M1'],('log_flower','display','petal_width'))
        self.assertEqual(m.MODEL_FEATURES['M2'][-1],'log_visit_rate')
        self.assertEqual(m.MODEL_FEATURES['M3'],('log_flower','log_visit_rate'))

    def test_all_five_shelters_held_out_exactly_once(self):
        out=m.score_shelters(self.rows)
        self.assertEqual(out['total_heldout_individuals_once_each'],60)
        self.assertEqual(len(out['folds']),5)
        self.assertEqual({r['shelter'] for r in out['folds']},set(m.SHELTERS))
        self.assertTrue(all(r['train_plants']+r['heldout_plants']==60 for r in out['folds']))

    def test_both_process_increments_are_reported(self):
        a=m.score_shelters(self.rows)
        self.assertEqual(set(a['gains']),{'trait_gain_m0_to_m1','process_gain_m1_to_m2','process_gain_m0_to_m3'})
        for k,v in a['gains'].items():
            self.assertEqual(0<=v['improved_shelters']<=5,True)
            self.assertTrue(math.isfinite(v['equal_shelter']))

    def test_score_rejects_missing_original_shelter(self):
        with self.assertRaisesRegex(m.SinapisGateError,'SHELTERS'):
            m.score_shelters([r for r in self.rows if r['shelter']!='4'])

    def test_shuffle_reference_seed_reproducible(self):
        a=m.shuffle_placebo(self.rows,reps=12,seed=20261105)
        b=m.shuffle_placebo(self.rows,reps=12,seed=20261105)
        self.assertEqual(a,b)

    def test_visitor_placebo_exposes_correct_original_units(self):
        p=m.shuffle_placebo(self.rows,reps=12,seed=123)
        self.assertEqual(p['unit'].split()[0],'shuffle')
        self.assertEqual(p['replicates'],12)
        self.assertFalse(p['is_randomized_causal_p_value'])

    def test_nonconstant_visits_respond_to_whole_site_permutation(self):
        y=m.score_shelters(self.rows)
        x=[dict(r) for r in self.rows]
        rng=random.Random(12)
        vals=[r['log_visit_rate'] for r in x]
        rng.shuffle(vals)
        for r,v in zip(x,vals):r['log_visit_rate']=v
        a=m.score_shelters(x)
        self.assertNotEqual(y['gains']['process_gain_m1_to_m2'],a['gains']['process_gain_m1_to_m2'])

    def test_no_visitation_variation_has_no_process_gain(self):
        r,_=m.source_join(corpus(none_visits=True),strict=False)
        out=m.score_shelters(r)
        self.assertAlmostEqual(out['gains']['process_gain_m1_to_m2']['equal_shelter'],0,places=9)

    def test_ridge_protection_singular_backdrops(self):
        data=[dict(x,log_flower=3.,display=2.,petal_width=1.) for x in self.rows]
        out=m.score_shelters(data)
        self.assertTrue(math.isfinite(out['gains']['process_gain_m1_to_m2']['equal_shelter']))

    def test_signed_negative_results_are_not_removed(self):
        r,_=m.source_join(corpus(pure_background=True),strict=False)
        out=m.score_shelters(r)
        self.assertIn('process_gain_m1_to_m2',out['gains'])
        self.assertEqual(len(out['folds']),5)

    def test_placebo_parameter_budgets_validated(self):
        with self.assertRaisesRegex(m.SinapisGateError,'INVALID_PLACEBO'):
            m.shuffle_placebo(self.rows,reps=1)
        with self.assertRaisesRegex(m.SinapisGateError,'INVALID_PLACEBO'):
            m.shuffle_placebo(self.rows,reps=1001)

    def test_exact_helper_quantiles(self):
        self.assertEqual(m._quantile([0.,2.,4.],.5),2.)
        self.assertEqual(m._quantile([0.,2.,4.],.25),1.)

    def test_original_data_not_needed_for_synthetic_testing(self):
        d=corpus()
        self.assertEqual(set(d),set(m.ORIGINAL_MD5))
        self.assertFalse(any(b'genotype' in x for x in d.values()))

if __name__=='__main__':unittest.main()
