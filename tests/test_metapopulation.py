"""Adversarial, conservation and temporal tests for the v5 synthetic process."""
from dataclasses import replace
from math import isclose
import unittest

from sdea.metapopulation import (
    Meta, cross, GENOTYPES, slot_distribution, weights,
    expression_probs, ecological_next, generation, ensemble,
    run_trajectory, allele_copies, exact_recolonization_probability_empty, exact_first_step_contrasts, monte_carlo_difference,
)
from random import Random

O=(0,0); A=(2,0); B=(0,2)


class MetapopulationTests(unittest.TestCase):
    def test_all_genotypes_mendelian_cross_rows_normalize(self):
        for m in GENOTYPES:
            for f in GENOTYPES:
                probs=cross(m,f)
                self.assertAlmostEqual(sum(probs),1,places=14)
                self.assertTrue(all(0<=v<=1 for v in probs))

    def test_mendelian_expected_offspring_copies(self):
        for m,f in ((A,B),((0,2),(2,0)),((1,0),(2,2)),((0,0),(1,1))):
            p=cross(m,f)
            for j in (0,1):
                expected=sum(g[j]*pr for g,pr in zip(GENOTYPES,p))
                self.assertAlmostEqual(expected,(m[j]+f[j])/2)

    def test_two_route_hybrid_from_distinct_parents(self):
        p=cross(A,B)
        self.assertAlmostEqual(p[GENOTYPES.index((1,1))],1.)

    def test_source_seed_is_meiotic_not_clone(self):
        par=replace(Meta(),seed_rates=(1.,0.),assurance=.5)
        # Mixed donor produces Mendelian segregating seeds, not adult clones.
        d=slot_distribution((),((1,1),(1,1)),0,1,0,par)
        imported={g:v for src,g,v in d if src=='seed_from_other'}
        self.assertGreater(imported.get((0,0),0),0)
        self.assertGreater(imported.get((2,2),0),0)
        self.assertGreater(imported.get((1,1),0),0)

    def test_slot_kernel_normalization_many_cells(self):
        for recipient in ((),(O,), (A,A),(A,B)):
            for donor in ((),(B,), (B,B),(A,B)):
                for i in (0,1):
                    for k in (0,1):
                        p=Meta(pollen_rates=(.45,.35),seed_rates=(.5,.15),assurance=.7)
                        for target in (0,1):
                            dist=slot_distribution(recipient,donor,i,k,target,p)
                            self.assertAlmostEqual(sum(x[2] for x in dist),1.,places=12)
                            self.assertTrue(all(0<=v<=1 for _,_,v in dist))

    def test_empty_cannot_make_local_zygote(self):
        p=Meta(pollen_rates=(1.,1.),seed_rates=(0.,0.))
        d=slot_distribution((),(A,A),0,1,0,p)
        self.assertEqual(d,(('empty',None,1.),))

    def test_seed_import_recolonizes_even_without_mother(self):
        p=Meta(seed_rates=(.5,0.))
        d=slot_distribution((),(A,A),0,1,0,p)
        self.assertGreater(sum(w for o,g,w in d if o=='seed_from_other'),0)
        self.assertFalse(any(o=='pollen_from_other' for o,g,w in d))

    def test_seed_recolonization_exact_formula_matches_slot(self):
        p=Meta(seed_rates=(.5,0.), catastrophe_rates=(.15,0.))
        d=slot_distribution((),(A,A),0,1,0,p)
        per=sum(w for o,g,w in d if o=='seed_from_other')
        expected=(1-(1-per)**p.capacity)*(1-p.catastrophe_rates[0])
        self.assertAlmostEqual(expected,exact_recolonization_probability_empty((A,A),1,.5,p))

    def test_donor_extinct_cannot_send_even_with_positive_rates(self):
        p=Meta(seed_rates=(1.,1.),pollen_rates=(1.,1.))
        d=slot_distribution((),(),0,1,0,p)
        self.assertEqual(d,(('empty',None,1.),))

    def test_both_empty_remain_empty_regardless_of_flow(self):
        p=Meta(seed_rates=(1.,1.),pollen_rates=(1.,1.))
        out=run_trajectory(((),()),(1,1),8,p,seed=88)
        self.assertEqual(out['populations'],((),()))
        self.assertFalse(any(out['entry'].values()))

    def test_no_innovation_without_initial_allele(self):
        p=Meta(seed_rates=(.8,.8),pollen_rates=(.8,.8))
        initial=((O,O),(B,B))
        for seed in range(50):
            r=run_trajectory(initial,(1,1),8,p,seed)
            self.assertTrue(all(t['allele_A']==(0,0) for t in r['path']))
            self.assertFalse(r['ever_A'][0] or r['ever_A'][1])

    def test_one_generation_pollen_can_build_joint_seed_cannot(self):
        pollen=Meta(pollen_rates=(.6,.6))
        seed=Meta(seed_rates=(.6,.6))
        pp=slot_distribution((A,A),(B,B),1,1,0,pollen)
        sp=slot_distribution((A,A),(B,B),1,1,0,seed)
        pollen_joint=sum(w for o,g,w in pp if g is not None and g[0] and g[1])
        seed_joint=sum(w for o,g,w in sp if g is not None and g[0] and g[1])
        self.assertGreater(pollen_joint,0.)
        self.assertEqual(seed_joint,0.)

    def test_neutral_no_movement_cannot_build_composite(self):
        p=Meta()
        initial=((A,A),(B,B))
        for i in range(40):
            r=run_trajectory(initial,(0,1),7,p,i)
            self.assertFalse(any(r['ever_joint_individual']))

    def test_seed_can_erode_local_A_with_B_only_donor_at_one_step(self):
        closed=slot_distribution((A,A),(B,B),1,1,0,Meta())
        imported=slot_distribution((A,A),(B,B),1,1,0,Meta(seed_rates=(.8,0.)))
        def per_slot_A(d):
            return sum(w for _,g,w in d if g is not None and g[0]>0)
        self.assertLess(per_slot_A(imported),per_slot_A(closed))
        def per_slot_viable(d):
            return sum(w for _,g,w in d if g is not None)
        self.assertGreater(per_slot_viable(imported),per_slot_viable(closed))

    def test_symmetric_dynamics_are_not_fixed_source(self):
        p=Meta(seed_rates=(.4,.4),pollen_rates=(.3,.3))
        x=run_trajectory(((A,A),(B,B)),(1,1),7,p,0)
        self.assertEqual(len(x['path']),7)
        self.assertTrue(any(row['allele_A'][1]!=0 or row['allele_B'][0]!=0 for row in x['path']))
        self.assertTrue(any(row['allele_A'][0]!=4 or row['allele_B'][1]!=4 for row in x['path']))

    def test_every_recolonization_starts_from_seed(self):
        p=Meta(seed_rates=(.5,.5),pollen_rates=(.5,.5),catastrophe_rates=(.15,.15))
        for seed in range(80):
            r=run_trajectory(((),(A,A)),(0,1),8,p,seed)
            self.assertGreaterEqual(r['recolon_count'][0],0)
            self.assertEqual(r['recolonized_after_extinction'][0],r['recolon_count'][0]>0)
            self.assertEqual(r['colonized_initial_vacancy'][0],r['colonization_count'][0]>0)

    def test_q_feedback_and_exclusive_route_probabilities(self):
        p=Meta(a_expression_q_feedback=.15,b_expression_q_feedback=-.15)
        self.assertAlmostEqual(sum(expression_probs((2,2),.5,p)),1.)
        self.assertEqual(expression_probs((0,0),.5,p),(1.,0.,0.))
        self.assertEqual(ecological_next(0,('A',),replace(p,low_becomes_high=1.),Random(1)),1)
        self.assertEqual(ecological_next(0,('B',),replace(p,low_becomes_high=0.),Random(1)),0)

    def test_q_zero_blocks_pollen_even_with_occupied_donor(self):
        p=Meta(q_low=0.,pollen_rates=(1.,0.),assurance=1.,inbreeding_depression=0.)
        x=slot_distribution((A,A),(B,B),0,1,0,p)
        self.assertEqual(sum(w for o,g,w in x if o=='pollen_from_other'),0.)
        self.assertGreater(sum(w for o,g,w in x if o=='self_local'),0.)

    def test_force_catastrophe_extinction_cannot_self_recolonize(self):
        p=Meta(seed_rates=(0.,0.), catastrophe_rates=(1.,1.))
        x=run_trajectory(((A,A),(B,B)),(1,1),4,p,seed=7)
        self.assertEqual(x['populations'],((),()))
        self.assertEqual(x['extinction_count'],(1,1))

    def test_generation_parentage_conservation_and_gene_entry(self):
        p=Meta(pollen_rates=(.6,.6),seed_rates=(.6,.6),catastrophe_rates=(0.,0.))
        for n in range(50):
            pops,q,vis,origins,entries=generation(((),(A,A)),(0,1),((),()),p,Random(n))
            self.assertLessEqual(len(pops[0]),2)
            for route,origin in entries[0]:
                self.assertEqual(route,'A')
                self.assertEqual(origin,'seed_from_other')

    def test_pollen_cannot_recolonize_empty_in_ensemble(self):
        p=Meta(pollen_rates=(1.,1.),seed_rates=(0.,0.))
        x=ensemble(((),(A,A)),p,horizon=4,replicates=100,seed=7)
        self.assertEqual(x['event_counts']['patch0_colonized_initial_vacancy'],0)
        self.assertEqual(x['event_counts']['patch0_occupied'],0)

    def test_positive_seed_recolonization_in_ensemble(self):
        p=Meta(seed_rates=(.6,.6))
        x=ensemble(((),(A,A)),p,horizon=4,replicates=160,seed=19)
        self.assertGreater(x['event_counts']['patch0_colonized_initial_vacancy'],0)
        self.assertEqual(x['event_counts']['patch0_ancestral_gene_entered_via_pollen'],0)
        self.assertGreater(x['event_counts']['patch0_ancestral_gene_entered_via_seed'],0)

    def test_reproducible_monte_carlo_and_valid_wilson_intervals(self):
        p=Meta(pollen_rates=(.2,.2),seed_rates=(.2,.2))
        x=ensemble(((A,A),(B,B)),p,replicates=100,horizon=4,seed=14)
        y=ensemble(((A,A),(B,B)),p,replicates=100,horizon=4,seed=14)
        self.assertEqual(x,y)
        for metric,v in x['probabilities'].items():
            a,b=x['wilson_95'][metric]
            self.assertTrue(0<=a<=v<=b<=1,(metric,v,a,b))

    def test_closed_world_route_histories_do_not_prove_ancestry(self):
        p=Meta(seed_rates=(.5,.5),pollen_rates=(.5,.5))
        r=run_trajectory(((A,A),(B,B)),(1,0),7,p,seed=36)
        self.assertEqual(r['horizon'],7)
        self.assertTrue(set(r['first_visible']).issubset({'none','A','B','both'}))
        for t in r['path']:
            self.assertTrue(all(0<=n<=4 for n in t['allele_A']))

    def test_first_step_exact_pollen_composite_seed_displacement(self):
        result=exact_first_step_contrasts()
        self.assertGreater(result['pollen_only']['any_joint_competence_generation1'],0.)
        self.assertEqual(result['seed_only']['any_joint_competence_generation1'],0.)
        self.assertEqual(result['no_flow']['any_joint_competence_generation1'],0.)
        self.assertLess(result['seed_only']['any_ancestral_support_generation1'],
                        result['no_flow']['any_ancestral_support_generation1'])
        self.assertGreater(result['seed_only']['any_occupancy_generation1'],
                           result['no_flow']['any_occupancy_generation1'])

    def test_monte_carlo_difference_is_second_minus_first(self):
        x=ensemble(((),(A,A)),Meta(seed_rates=(.25,.25)),horizon=3,replicates=120,seed=2)
        y=ensemble(((),(A,A)),Meta(seed_rates=(.45,.45)),horizon=3,replicates=120,seed=3)
        c=monte_carlo_difference({'x':x,'y':y},'x','y','patch0_colonized_initial_vacancy')
        self.assertAlmostEqual(c['difference_second_minus_first'],
             y['probabilities'][c['endpoint']]-x['probabilities'][c['endpoint']])
        self.assertGreaterEqual(c['mc_standard_error'],0.)

    def test_invalid_parameters_fail_closed(self):
        for p in (dict(capacity=0),dict(pollen_rates=(.5,)),dict(seed_rates=(.1,1.2)),
                  dict(high_stays_high=-.01),dict(q_low=.9,q_high=.1),
                  dict(catastrophe_rates=(0.1, float('nan'))),
                  dict(a_expression_q_feedback=1.1)):
            with self.assertRaises(ValueError):Meta(**p)
        with self.assertRaises(ValueError):slot_distribution((A,A,A), (B,B),0,1,0)
        with self.assertRaises(ValueError):run_trajectory(((A,A),(B,B)),horizon=0)
        with self.assertRaises(ValueError):cross((3,0),A)

if __name__=='__main__':unittest.main()
