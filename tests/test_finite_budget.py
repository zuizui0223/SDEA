"""SDEA v6: exact budget controls, causal boundaries and MC integrity tests."""
import json
import unittest
from dataclasses import replace
from pathlib import Path
from random import Random
from math import isclose
from sdea.finite_budget import (
    A, B, O, Budget, GENOTYPES, IDENTITY, _cross, _father, _population,
    _settle, allele_count, successful_seed_probability,
    exact_colonization_from_one_source, exact_shared_budget_colonization,
    generation, run_trajectory, ensemble, difference, benchmark,
)

STRONG=((.2,.4,.4),(.4,.2,.4),(.4,.4,.2))
MODERATE=((.65,.175,.175),(.175,.65,.175),(.175,.175,.65))


class FiniteBudgetTests(unittest.TestCase):
    def test_validated_routing_rows_and_capacity(self):
        for kw in (dict(seed_routing=((1.,0.,0.),)),
                   dict(pollen_routing=((.9,.1,0.),(0.,1.,0.),(0.,0.,.8))),
                   dict(seed_routing=((.5,-.1,.6),IDENTITY[1],IDENTITY[2])),
                   dict(seed_attempts=-1),dict(seed_attempts=1.5),
                   dict(capacity=0),dict(capacity=True),
                   dict(assurance=1.1),dict(q_low=.9,q_high=.7),
                   dict(selection_a=float('nan')),
                   dict(catastrophe_rates=(0.,0.)),
                   dict(expression_feedback_a=1.2)):
            with self.subTest(kw=kw),self.assertRaises(ValueError):
                Budget(**kw)

    def test_population_invalid_and_allele_counts(self):
        _population((A,(1,1)),2)
        self.assertEqual(allele_count((A,(1,1)),0),3)
        self.assertEqual(allele_count((A,(1,1)),1),1)
        for state in (((3,0),),((True,0),), (A,A,A), [A]):
            with self.subTest(state=state), self.assertRaises(ValueError):
                _population(state,2)
        with self.assertRaises(ValueError):
            allele_count((A,),2)

    def test_mendelian_AxB_combines_routes_in_one_hybrid(self):
        for seed in range(30):
            self.assertEqual(_cross(A,B,Random(seed)),(1,1))
            self.assertEqual(_cross(B,A,Random(seed)),(1,1))

    def test_single_mother_cannot_local_outcross(self):
        j,dad=_father(((A,),(),()),0,0,IDENTITY,Random(2))
        self.assertIsNone(j)
        self.assertIsNone(dad)
        p=Budget(assurance=0.)
        self.assertEqual(successful_seed_probability(((A,),(),()),(1,0,0),0,p),0.)

    def test_incoming_pollen_gamete_requires_living_mother(self):
        pp=((0.,1.,0.),(0.,1.,0.),(0.,0.,1.))
        p=Budget(pollen_routing=pp,assurance=0.,q_high=1.,q_low=0.)
        pops=((A,),(B,B),())
        self.assertAlmostEqual(successful_seed_probability(pops,(1,1,0),0,p),1.)
        new,q,v,ledger=generation(pops,(1,1,0),((),(),()),p,Random(8))
        self.assertTrue(new[0])
        self.assertTrue(all(g==(1,1) for g in new[0]))
        self.assertEqual(new[2],())
        self.assertTrue(all(i==0 for i in ledger['attempts'][2:]))

    def test_pollen_cannot_colonize_two_empty_destinations(self):
        p=Budget(seed_routing=IDENTITY,pollen_routing=STRONG,seed_attempts=8)
        for k in range(40):
            r=run_trajectory(((A,A),(),()),par=p,horizon=5,seed=k)
            self.assertEqual(r['populations'][1],())
            self.assertEqual(r['populations'][2],())
            self.assertFalse(any(r['colonized_initial_vacancy'][1:]))

    def test_budget_zero_has_no_recruits_and_no_seed_reentry(self):
        p=Budget(seed_attempts=0,seed_routing=STRONG,pollen_routing=STRONG)
        x=run_trajectory(((A,A),(),()),par=p,horizon=4,seed=22)
        self.assertEqual(x['populations'],((),(),()))
        self.assertEqual(x['total_viable_seed_produced'],0)
        self.assertEqual(x['total_seed_exports'],0)

    def test_exact_viability_one_source_q_high(self):
        p=Budget()
        s=.45*(1-.75)
        expected=s*.8+(1-s)*.75
        self.assertAlmostEqual(successful_seed_probability(((A,A),(),()),(1,0,0),0,p),expected,places=14)
        self.assertAlmostEqual(expected,.755625)

    def test_exact_single_source_colonization_formula(self):
        p=Budget(seed_attempts=3,seed_routing=STRONG)
        x=exact_colonization_from_one_source(p,0,1)
        S=x['source_viability_per_attempt']
        self.assertAlmostEqual(x['colonized_after_catastrophe'],1-(1-S*.4)**3)
        y=exact_colonization_from_one_source(replace(p,catastrophe_rates=(0.,.3,0.)),0,1)
        self.assertAlmostEqual(y['colonized_after_catastrophe'],x['colonized_after_catastrophe']*.7)
        with self.assertRaises(ValueError):
            exact_colonization_from_one_source(p,0,0)
        with self.assertRaises(ValueError):
            exact_colonization_from_one_source(p,0,1,count_donor=3)

    def test_shared_source_negative_colonization_covariance(self):
        for b in (1,2,3,5,8):
            r=exact_shared_budget_colonization(Budget(seed_attempts=b,seed_routing=STRONG))
            self.assertLessEqual(r['joint_colonization'],r['naive_independent_joint']+1e-12)
            self.assertGreaterEqual(r['joint_colonization'],-1e-12)
            self.assertTrue(all(0<=x<=1 for x in r['recipient_colonization']))
        r=exact_shared_budget_colonization(Budget(seed_attempts=3,seed_routing=STRONG))
        self.assertAlmostEqual(r['joint_colonization'],0.38245796915625,places=12)
        self.assertAlmostEqual(r['naive_independent_joint'],0.435991969923443,places=12)
        self.assertLess(r['joint_minus_independent'],-.05)

    def test_budget_one_cannot_colonize_two_patches(self):
        p=Budget(seed_attempts=1,seed_routing=STRONG)
        self.assertEqual(exact_shared_budget_colonization(p)['joint_colonization'],0.)
        for seed in range(200):
            pops,*_=generation(((A,A),(),()),(1,0,0),((),(),()),p,Random(seed))
            self.assertFalse(pops[1] and pops[2])

    def test_budget_two_can_split_two_seeds_with_probability_half(self):
        m=((0.,.5,.5),IDENTITY[1],IDENTITY[2])
        par=Budget(seed_attempts=2,seed_routing=m,assurance=0.,q_high=1.,q_low=0.)
        exact=exact_shared_budget_colonization(par)
        self.assertAlmostEqual(exact['joint_colonization'],.5)
        self.assertAlmostEqual(exact['recipient_colonization'][0],.75)
        self.assertLess(exact['joint_colonization'],.75**2)

    def test_single_source_monte_carlo_matches_joint_exact(self):
        p=Budget(seed_attempts=3,seed_routing=STRONG)
        n=8000
        counts=[0,0,0]
        for seed in range(n):
            pops,*_=generation(((A,A),(),()),(1,0,0),((),(),()),p,Random(seed))
            counts[0]+=bool(pops[1]);counts[1]+=bool(pops[2]);counts[2]+=bool(pops[1] and pops[2])
        exact=exact_shared_budget_colonization(p)
        for emp,target in zip(counts,[*exact['recipient_colonization'],exact['joint_colonization']]):
            self.assertLess(abs(emp/n-target),.023)

    def test_seed_conservation_and_no_clone_transport(self):
        p=Budget(seed_attempts=4,seed_routing=STRONG,pollen_routing=MODERATE)
        old=((A,(1,1)),(B,B),())
        for seed in range(50):
            nxt, q, vis, log=generation(old,(1,1,0),((),(),()),p,Random(seed))
            self.assertEqual(log['attempts'],(4,4,0))
            self.assertLessEqual(sum(log['viable_produced']),8)
            self.assertEqual(sum(log['viable_produced']),sum(log['pre_recruit_arrivals']))
            self.assertEqual(sum(log['viable_produced']),sum(log['retained'])+sum(log['exports']))
            self.assertLessEqual(sum(len(patch) for patch in nxt),sum(log['viable_produced']))
            for dst, admits in enumerate(log['after_recruit_admissions']):
                self.assertLessEqual(len(admits),p.capacity)
                for genotype,mother,father in admits:
                    self.assertIn(genotype,GENOTYPES)
                    self.assertIn(mother,range(3))
                    self.assertIn(father,range(3))

    def test_all_routes_absent_are_absorbing_without_mutation(self):
        p=Budget(seed_routing=STRONG,pollen_routing=STRONG,seed_attempts=10)
        initial=((O,O),((0,2),(0,2)),())
        for seed in range(40):
            outcome=run_trajectory(initial,par=p,horizon=6,seed=seed)
            self.assertTrue(all(row['A_copies']==(0,0,0) for row in outcome['path']))
            self.assertFalse(any(outcome['ever_visible_A']))

    def test_positive_seed_export_can_lower_local_first_step_retention(self):
        local=Budget(seed_routing=IDENTITY)
        export=Budget(seed_routing=STRONG)
        p0=1-(1-successful_seed_probability(((A,A),(),()),(1,0,0),0,local))**3
        v=exact_shared_budget_colonization(export)['source_retention']
        self.assertLess(v,p0)
        self.assertGreater(exact_shared_budget_colonization(export)['recipient_colonization'][0],0.)
        self.assertAlmostEqual(v,1-(1-.755625*.2)**3)

    def test_stronger_budget_increases_exact_source_retention(self):
        a=exact_shared_budget_colonization(Budget(seed_attempts=3,seed_routing=STRONG))
        b=exact_shared_budget_colonization(Budget(seed_attempts=8,seed_routing=STRONG))
        self.assertGreater(b['source_retention'],a['source_retention'])
        self.assertGreater(b['joint_colonization'],a['joint_colonization'])

    def test_positive_selection_changes_lottery_only_with_competition(self):
        pool=((A,0,0),(B,2,2))
        base=Budget(capacity=1,selection_a=0.)
        strong=replace(base,selection_a=6.)
        self.assertEqual(_settle(pool,.75,2,strong,Random(0)),pool)
        n=1000
        selected0=sum(_settle(pool,.75,1,base,Random(i))[0][0]==A for i in range(n))
        selected1=sum(_settle(pool,.75,1,strong,Random(i))[0][0]==A for i in range(n))
        self.assertLess(abs(selected0/n-.5),.06)
        self.assertGreater(selected1/n,.92)

    def test_catastrophe_erases_recruits_even_after_successful_seed_production(self):
        p=Budget(seed_routing=STRONG,seed_attempts=6,catastrophe_rates=(1.,1.,1.))
        new,_,_,ledger=generation(((A,A),(),()),(1,0,0),((),(),()),p,Random(2))
        self.assertEqual(new,((),(),()))
        self.assertTrue(all(len(x)==0 for x in ledger['after_recruit_admissions']))
        self.assertGreater(sum(ledger['viable_produced']),0.)

    def test_q_feedback_uses_previous_visible_not_future_endpoint(self):
        p=Budget(seed_routing=IDENTITY,high_stays_high=0.,low_becomes_high=0.,expression_feedback_a=1.)
        new,q,vis,_=generation(((A,A),(),()),(0,0,0),(('A',),(),()),p,Random(2))
        self.assertEqual(q[0],1)
        b=replace(p,expression_feedback_a=0.)
        _,q2,_,_=generation(((A,A),(),()),(0,0,0),(('A',),(),()),b,Random(2))
        self.assertEqual(q2[0],0)

    def test_first_colonization_different_from_recolonization_after_loss(self):
        initial=((A,A),(),())
        p=Budget(seed_routing=STRONG,catastrophe_rates=(.2,.2,.2),seed_attempts=8)
        r=run_trajectory(initial,par=p,horizon=10,seed=28)
        self.assertEqual(len(r['colonized_initial_vacancy']),3)
        self.assertEqual(len(r['recolonized_after_extinction']),3)
        self.assertFalse(r['colonized_initial_vacancy'][0])
        for row in r['path']:
            self.assertLessEqual(sum(row['viable_produced']),sum(row['seed_attempts']))

    def test_mc_seeded_replay_and_wilson_bounds(self):
        p=Budget(seed_attempts=4,seed_routing=STRONG)
        one=ensemble(((A,A),(),()),p,horizon=4,replicates=120,seed=51)
        two=ensemble(((A,A),(),()),p,horizon=4,replicates=120,seed=51)
        self.assertEqual(one,two)
        for key, probability in one['probabilities'].items():
            lo,hi=one['mc_wilson_95'][key]
            self.assertTrue(0<=lo<=probability<=hi<=1,(key,probability,lo,hi))

    def test_mc_difference_correct_orientation(self):
        low=ensemble(((A,A),(),()),Budget(seed_routing=IDENTITY),horizon=3,replicates=100,seed=7)
        high=ensemble(((A,A),(),()),Budget(seed_routing=STRONG),horizon=3,replicates=100,seed=8)
        r=difference({'low':low,'high':high},'low','high','patch1_first_colonization')
        self.assertEqual(r['p_first'],0.)
        self.assertGreater(r['difference_second_minus_first'],0.)
        self.assertAlmostEqual(r['difference_second_minus_first'],r['p_second']-r['p_first'])

    def test_generation_rejects_nonfinite_and_malformed_states(self):
        for pops,qi,visible in (
            (((A,A),(),()),(1,0),((),(),())),
            (((A,A),(),()),(1,2,0),((),(),())),
            (((A,A),(),()),(1,0,0),((),())),
            (((A,A),(),()),(1,0,0),('A',(),())),
            (((A,A),(),()),(1,0,0),(('unknown',),(),())),
        ):
            with self.subTest(pops=pops,qi=qi,visible=visible),self.assertRaises(ValueError):
                generation(pops,qi,visible,Budget(),Random(1))

    def test_exact_joint_fails_closed_for_invalid_recipient_pair(self):
        p=Budget()
        for recipients in ((0,1),(1,1),(1,3),(1,),(1,2,0)):
            with self.subTest(recipients=recipients),self.assertRaises(ValueError):
                exact_shared_budget_colonization(p,recipients=recipients)

    def test_one_source_q_zero_and_no_assurance_blocks_reproduction(self):
        p=Budget(q_low=0.,assurance=0.,seed_routing=STRONG)
        r=exact_shared_budget_colonization(p,donor_q='low')
        self.assertEqual(r['joint_colonization'],0.)
        self.assertEqual(r['source_retention'],0.)

    def test_fixture_match_and_reproducible_three_patch_counts(self):
        result=benchmark(replicates=80)
        result2=benchmark(replicates=80)
        self.assertEqual(result,result2)
        self.assertEqual(len(result['scenarios']),8)
        self.assertEqual(len(result['contrasts']),7)
        self.assertEqual(result['scenarios']['single_local_B3']['event_counts']['patch1_first_colonization'],0)
        self.assertEqual(result['scenarios']['single_local_B3']['event_counts']['patch2_first_colonization'],0)

if __name__=='__main__':
    unittest.main()
