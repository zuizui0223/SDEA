"""SDEA v7b: exact graph/frontier controls and independently seeded contrasts."""
import itertools
import json
from dataclasses import replace
from pathlib import Path
from random import Random
import unittest

from sdea.network_reachability import (
    A, B, O, IDENTITY, CYCLE, RETURN, REVERSE, BALANCED, Scenario,
    _choose, _zygote, _routing, viability, reachable_sources,
    step, exact_matched_two_step, exhaustive_two_step, operator_brier_gain,
    run_history, ensemble, mc_difference,
)

SINGLE=((A,A),(),())
DISTINCT=((A,A),(),(B,B))


class NetworkTopologyTests(unittest.TestCase):
    def test_budget_and_seed_retention_identical_across_networks(self):
        for network in (CYCLE,RETURN,REVERSE,BALANCED):
            for i,row in enumerate(network):
                self.assertAlmostEqual(sum(row),1.)
                self.assertAlmostEqual(row[i],.6)
                self.assertAlmostEqual(1-row[i],.4)
        self.assertEqual(CYCLE[0],RETURN[0])
        self.assertEqual(CYCLE[2],RETURN[2])

    def test_exact_mother_success(self):
        self.assertEqual(viability(0),0.)
        self.assertAlmostEqual(viability(1),.09)
        self.assertAlmostEqual(viability(2),.755625)
        self.assertAlmostEqual(viability(4),.755625)

    def test_t1_source_and_receiver_identical_by_formula(self):
        forward=exact_matched_two_step(Scenario(seed_routing=CYCLE))
        backward=exact_matched_two_step(Scenario(seed_routing=RETURN))
        self.assertEqual(forward['t1_patch0_occupied'],backward['t1_patch0_occupied'])
        self.assertEqual(forward['t1_patch1_colonized'],backward['t1_patch1_colonized'])
        self.assertAlmostEqual(forward['t1_patch1_colonized'],.660296880140625)
        self.assertEqual(forward['t1_patch2_colonized'],0.)

    def test_exact_t2_reachable_vs_inaccessible(self):
        forward=exact_matched_two_step(Scenario(seed_routing=CYCLE))
        backward=exact_matched_two_step(Scenario(seed_routing=RETURN))
        self.assertAlmostEqual(forward['t2_patch2_colonized'],.19048150974644207,places=13)
        self.assertEqual(backward['t2_patch2_colonized'],0.)
        self.assertEqual(forward['reachability'],[0,1,2])
        self.assertEqual(backward['reachability'],[0,1])

    def test_exhaustive_allocation_proof_independent(self):
        for route in (CYCLE,RETURN):
            q=Scenario(seed_routing=route)
            self.assertAlmostEqual(exact_matched_two_step(q)['t2_patch2_colonized'],
                                   exhaustive_two_step(q),places=14)

    def test_known_operator_strictly_improves_full_state_forecast(self):
        p=operator_brier_gain()
        self.assertAlmostEqual(p['operator_hidden_brier'],.07019017811220257)
        self.assertAlmostEqual(p['operator_observed_brier'],.045139601351184125)
        self.assertAlmostEqual(p['brier_gain_from_observing_operator'],.02505057676101845)
        self.assertGreater(p['brier_gain_from_observing_operator'],0.)
        self.assertAlmostEqual(sum(v['probability'] for v in p['conditional_states']),1.)

    def test_operator_brier_zero_when_same_future_law(self):
        p=operator_brier_gain(Scenario(seed_routing=CYCLE),Scenario(seed_routing=CYCLE))
        self.assertAlmostEqual(p['brier_gain_from_observing_operator'],0.,places=14)

    def test_operator_brier_rejects_changes_in_genetic_or_q_law(self):
        with self.assertRaises(ValueError):
            operator_brier_gain(Scenario(seed_routing=CYCLE),Scenario(seed_routing=RETURN,q=.8))
        with self.assertRaises(ValueError):
            operator_brier_gain(Scenario(seed_routing=CYCLE),Scenario(seed_routing=REVERSE))

    def test_short_budget_edges(self):
        for k in (0,1,2,3,4,8):
            q=Scenario(seed_routing=CYCLE,seed_attempts=k)
            a=exact_matched_two_step(q)
            self.assertTrue(0<=a['t2_patch2_colonized']<=1)
            self.assertAlmostEqual(a['t2_patch2_colonized'],exhaustive_two_step(q))
        self.assertEqual(exact_matched_two_step(Scenario(seed_attempts=0))['t2_patch2_colonized'],0.)

    def test_global_route_inaccessibility_with_source_pollen(self):
        # Even incoming paternal pollen cannot create a mother in an empty
        # patch without any directed incoming seed path.
        mixed=((1/3,1/3,1/3),)*3
        q=Scenario(seed_routing=RETURN,pollen_routing=mixed)
        for seed in range(100):
            x=run_history(SINGLE,q,7,seed)
            self.assertFalse(x['ever_occupied'][2])
            self.assertFalse(x['final'][2])

    def test_t1_full_population_state_equal_under_same_seed(self):
        # Row 1 is the sole changed routing row, and initially has no mother.
        p=Scenario(seed_routing=CYCLE)
        q=replace(p,seed_routing=RETURN)
        for seed in range(80):
            self.assertEqual(run_history(SINGLE,p,1,seed),run_history(SINGLE,q,1,seed))

    def test_t2_divergence_is_possible_under_same_seed(self):
        p=Scenario(seed_routing=CYCLE)
        q=replace(p,seed_routing=RETURN)
        witnesses=0
        for seed in range(220):
            a=run_history(SINGLE,p,2,seed)
            b=run_history(SINGLE,q,2,seed)
            witnesses+=bool(a['final'][2] and not b['final'][2])
        self.assertGreater(witnesses,0)

    def test_three_patch_reachability_uses_direction(self):
        self.assertEqual(reachable_sources(CYCLE,(0,)),(0,1,2))
        self.assertEqual(reachable_sources(RETURN,(0,)),(0,1))
        self.assertEqual(reachable_sources(BALANCED,(0,)),(0,1,2))
        self.assertEqual(reachable_sources(RETURN,(0,2)),(0,1,2))

    def test_reverse_is_permutation_of_cycle_for_symmetric_founders(self):
        swap=(2,1,0)
        self.assertEqual(REVERSE,tuple(tuple(CYCLE[swap[i]][swap[j]] for j in range(3)) for i in range(3)))
        self.assertEqual(Scenario(seed_routing=CYCLE).q,Scenario(seed_routing=REVERSE).q)

    def test_mendelian_gene_copy_expectations(self):
        rng=Random(70)
        for _ in range(800):
            child=_zygote(A,B,rng)
            self.assertEqual(child,(1,1))
        for _ in range(100):
            self.assertEqual(_zygote(O,O,rng),(0,0))

    def test_no_spontaneous_B_from_all_A_founders(self):
        for seed in range(40):
            history=run_history(SINGLE,Scenario(seed_routing=BALANCED),7,seed)
            self.assertTrue(all(all(x==0 for x in t['B']) for t in history['trajectory']))
            self.assertTrue(all(not (a and b) for p in history['final'] for a,b in p))

    def test_seed_budget_mass_ledger(self):
        for initial in (SINGLE,DISTINCT):
            for net in (CYCLE,RETURN,BALANCED):
                pops=initial
                par=Scenario(seed_routing=net)
                rng=Random(41)
                for _ in range(7):
                    parents=pops
                    pops,ledger=step(pops,par,rng)
                    self.assertEqual(sum(ledger['viable']),sum(ledger['arrived']))
                    self.assertTrue(all(0<=x<=par.seed_attempts for x in ledger['attempted']))
                    self.assertEqual(tuple(par.seed_attempts if p else 0 for p in parents),ledger['attempted'])
                    self.assertTrue(all(v<=a for v,a in zip(ledger['viable'],ledger['attempted'])))
                    self.assertTrue(all(len(p)<=par.capacity for p in pops))

    def test_no_seed_production_or_propagation_with_zero_budget(self):
        q=Scenario(seed_attempts=0)
        x=run_history(SINGLE,q,4,0)
        self.assertEqual(x['final'],((),(),()))
        self.assertFalse(x['ever_occupied'][1] or x['ever_occupied'][2])

    def test_all_initially_empty_is_absorbing(self):
        for seed in range(30):
            x=run_history(((),(),()),Scenario(seed_routing=BALANCED),5,seed)
            self.assertEqual(x['final'],((),(),()))

    def test_no_pollen_viability_when_q_zero_and_no_assurance(self):
        q=Scenario(q=0.,assurance=0.)
        self.assertEqual(viability(2,q),0.)
        for seed in range(30):
            x=run_history(SINGLE,q,3,seed)
            self.assertEqual(x['final'],((),(),()))

    def test_binomial_exhaustive_works_with_no_recruits(self):
        q=Scenario(assurance=0.,q=0.,seed_routing=CYCLE)
        self.assertEqual(exact_matched_two_step(q)['t2_patch2_colonized'],0.)
        self.assertEqual(exhaustive_two_step(q),0.)

    def test_mc_close_to_prospective_exact_two_step(self):
        exact=exact_matched_two_step(Scenario(seed_routing=CYCLE))
        n=2400
        x=ensemble(SINGLE,Scenario(seed_routing=CYCLE),horizon=2,replicates=n,seed=1212)
        self.assertLess(abs(x['probabilities']['t2_patch2_occupied']-exact['t2_patch2_colonized']),.035)
        y=ensemble(SINGLE,Scenario(seed_routing=RETURN),horizon=2,replicates=n,seed=9887)
        self.assertEqual(y['event_counts']['t2_patch2_occupied'],0)

    def test_replicates_are_whole_landscapes_and_seeded(self):
        par=Scenario(seed_routing=BALANCED)
        a=ensemble(DISTINCT,par,replicates=80,seed=99)
        b=ensemble(DISTINCT,par,replicates=80,seed=99)
        self.assertEqual(a,b)
        self.assertEqual(a['replicates'],80)
        for name,value in a['probabilities'].items():
            l,u=a['mc_wilson_95'][name]
            self.assertTrue(-1e-12<=l<=value<=u+1e-12<=1+1e-12)

    def test_simulation_contrast_sign_convention(self):
        c=ensemble(SINGLE,Scenario(seed_routing=CYCLE),horizon=2,replicates=80,seed=77)
        r=ensemble(SINGLE,Scenario(seed_routing=RETURN),horizon=2,replicates=80,seed=111)
        d=mc_difference({'r':r,'c':c},'r','c','t2_patch2_occupied')
        self.assertAlmostEqual(d['difference_right_minus_left'],c['probabilities']['t2_patch2_occupied'])
        self.assertTrue(d['simulation_95_interval'][0]<=d['difference_right_minus_left']<=d['simulation_95_interval'][1])

    def test_prespecified_negative_null_has_no_claim_of_equivalence(self):
        # Same two A sources but reversed cycle: distributionally equivalent
        # under patch0<->patch2 permutation; finite Monte Carlo differs.
        a=ensemble(((A,A),(),(A,A)),Scenario(seed_routing=CYCLE),replicates=250,seed=3)
        b=ensemble(((A,A),(),(A,A)),Scenario(seed_routing=REVERSE),replicates=250,seed=14)
        d=mc_difference({'a':a,'b':b},'a','b','patch1_ever_occupied')
        self.assertLess(abs(d['difference_right_minus_left']),.1)

    def test_bad_inputs_fail_closed(self):
        with self.assertRaises(ValueError): Scenario(capacity=0)
        with self.assertRaises(ValueError): Scenario(seed_attempts=41)
        with self.assertRaises(ValueError): Scenario(seed_attempts=2.0)
        with self.assertRaises(ValueError): Scenario(q=-.1)
        with self.assertRaises(ValueError): Scenario(seed_routing=((.5,.5,0),(.4,.4,.1),(.4,0,.6)))
        with self.assertRaises(ValueError): Scenario(seed_routing=((1,0,0),(0,1,0)))
        with self.assertRaises(ValueError): run_history(((A,A),(),(A,B)),Scenario(),0,1)
        with self.assertRaises(ValueError): reachable_sources(CYCLE,())
        with self.assertRaises(ValueError): exact_matched_two_step(Scenario(seed_routing=BALANCED))
        with self.assertRaises(ValueError): exhaustive_two_step(Scenario(seed_attempts=13))

    def test_no_cross_family_pooling_in_the_fixture(self):
        data=json.loads((Path(__file__).resolve().parents[1]/'results/v7b_network_reachability_mc.json').read_text())
        self.assertEqual(set(data['scenarios']),{
            'single_cycle','single_return','single_balanced',
            'distinct_cycle','distinct_return','distinct_balanced',
            'neutral_cycle','neutral_reverse'})
        self.assertEqual({v['replicates'] for v in data['scenarios'].values()},{4000})
        self.assertEqual(data['replicate_unit'],'whole 3-patch history; separate founder families not pooled')
        self.assertEqual(len(data['contrasts']),6)

if __name__=='__main__': unittest.main()
