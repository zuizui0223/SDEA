"""SDEA v8: exact pool proof, mechanism stages, adverse and simulation checks."""
import math
import unittest
from dataclasses import replace
from random import Random

from sdea.network_topology import A, B, INIT, Network, routing, _first_hybrid_exact
from sdea.selection_filter import (
    _paired_summary, benchmark, cases, exact_pool_joint_recruitment,
    exact_three_seed_witness, generation_selected, joint, production, settle,
    trajectory,
)


class TestSelectionFilter(unittest.TestCase):
    def test_pool_identity_exact_all_predeclared_weights(self):
        for w in (0., 1/3, 1., 3.):
            row=exact_three_seed_witness(w)
            self.assertAlmostEqual(row['conditional_joint_recruitment'],
                                   row['enumerated_cross_check'],places=14)
        self.assertEqual(exact_three_seed_witness(0.)['conditional_joint_recruitment'],0.)
        self.assertAlmostEqual(exact_three_seed_witness(1.)['conditional_joint_recruitment'],2/3)
        self.assertAlmostEqual(exact_three_seed_witness(3.)['conditional_joint_recruitment'],.9)
        self.assertAlmostEqual(exact_three_seed_witness(1/3)['conditional_joint_recruitment'],5/14)

    def test_exact_pool_enumeration_handles_replicated_genotypes(self):
        self.assertEqual(exact_pool_joint_recruitment((), 2, 3.),0.)
        self.assertEqual(exact_pool_joint_recruitment((A,B),2,3.),0.)
        self.assertEqual(exact_pool_joint_recruitment(((1,1),),2,0.),0.)
        self.assertEqual(exact_pool_joint_recruitment(((1,1),),2,3.),1.)
        self.assertEqual(exact_pool_joint_recruitment(((1,1),(1,1),A,B),4,1.),1.)
        self.assertAlmostEqual(exact_pool_joint_recruitment(((1,1),A,B),1,3.),3/5)

    def test_competitive_selection_changes_mixed_pool_in_declared_direction(self):
        p=[exact_pool_joint_recruitment(((1,1),A,B),2,w) for w in (0.,1/3,1.,3.)]
        self.assertEqual(p[0],0.)
        self.assertEqual(p,sorted(p))
        self.assertGreater(p[-1],p[-2])
        self.assertLess(p[1],p[2])

    def test_stochastic_settlement_matches_exact_probability(self):
        pool=((1,1), A, B)
        for w in (1/3,3.):
            expected=exact_pool_joint_recruitment(pool,2,w)
            n=2200
            p=sum(any(joint(g) for g in settle(pool,2,w,Random(12345+i*53))) for i in range(n))/n
            se=math.sqrt(expected*(1-expected)/n)
            self.assertLess(abs(p-expected),4.5*se+.004)

    def test_hard_lethal_gate_removes_joint_even_without_competition(self):
        self.assertEqual(settle(((1,1),),2,0.,Random(1)),())
        self.assertEqual(settle(((1,1),A),2,0.,Random(1)),(A,))
        self.assertEqual(settle(((1,1),A),2,3.,Random(1)),tuple(sorted(((1,1),A))))

    def test_fit_positive_no_overflow_is_inert(self):
        pool=(A,(1,1))
        for w in (1/3,1.,3.):
            self.assertEqual(settle(pool,2,w,Random(21)),tuple(sorted(pool)))

    def test_invalid_parameters_fail_closed(self):
        for w in (-1, float('nan'),float('inf'),31,True, '3'):
            with self.assertRaises(ValueError):exact_three_seed_witness(w)
        with self.assertRaises(ValueError):exact_pool_joint_recruitment((A,),0,1.)
        with self.assertRaises(ValueError):exact_pool_joint_recruitment(((4,0),),2,1.)
        with self.assertRaises(ValueError):settle(((3,0),),2,1.,Random(0))
        with self.assertRaises(ValueError):trajectory(Network(),1.,horizon=0)
        with self.assertRaises(ValueError):trajectory(Network(),1.,seed=1.2)
        with self.assertRaises(ValueError):benchmark(replicates=20)

    def test_model_routing_matching_is_not_graph_only(self):
        a,b=routing('clustered'),routing('alternating')
        self.assertNotEqual(a,b)
        for graph in (a,b,routing('components')):
            for i in range(6):
                self.assertAlmostEqual(sum(graph[i]),1.)
                self.assertAlmostEqual(sum(row[i] for row in graph),1.)
                self.assertEqual(sum(x>0 for x in graph[i]),3)
                self.assertAlmostEqual(graph[i][i],.25)

    def test_mendelian_one_cross_has_two_route_alleles(self):
        from sdea.network_topology import _cross
        for i in range(30):
            self.assertEqual(_cross(A,B,Random(i)),(1,1))

    def test_first_generation_seed_production_is_exactly_selection_invariant(self):
        par=Network(seed=routing('clustered'),pollen=routing('alternating'))
        for i in range(120):
            before=production(INIT,par,Random(i))
            # The downstream selector has no parameter in this production law.
            for w in (0.,1/3,1.,3.):
                after,ledger=generation_selected(INIT,par,w,i,71+i)
                self.assertEqual(before[1]['source_joint_seeds'],ledger['source_joint_seeds'])
                self.assertEqual(before[1]['source_viable_seeds'],ledger['source_viable_seeds'])
                self.assertLessEqual(sum(map(len,after)),sum(before[1]['source_viable_seeds']))

    def test_no_cross_links_prevents_joint_production(self):
        par=Network(seed=routing('components'),pollen=routing('components'))
        self.assertEqual(_first_hybrid_exact(par)['probability_any_hybrid_seed_produced'],0.)
        for w in (0.,1/3,1.,3.):
            for i in range(30):
                v=trajectory(par,w,horizon=6,seed=17+i)
                self.assertFalse(v['ever_joint_seed_produced'])
                self.assertFalse(v['ever_joint_adult_recruited'])

    def test_lethal_selection_does_not_undo_zigote_formation(self):
        par=Network(seed=routing('alternating'),pollen=routing('alternating'))
        self.assertGreater(_first_hybrid_exact(par)['probability_any_hybrid_seed_produced'],.99)
        for i in range(25):
            v=trajectory(par,0.,7,i)
            self.assertFalse(v['gen1_joint_adult_recruited'])
            self.assertFalse(v['ever_joint_adult_recruited'])
            self.assertFalse(v['final_joint_adult'])

    def test_positive_competition_does_not_create_missing_allele(self):
        par=Network(seed=routing('alternating'),pollen=routing('alternating'))
        init=((A,A),)*6
        for i in range(20):
            v=trajectory(par,3.,7,seed=i,initial=init)
            self.assertFalse(v['ever_joint_seed_produced'])
            self.assertFalse(v['final_both_alleles_global'])

    def test_production_seeds_are_conserved_and_not_multiplied(self):
        par=Network(seed=routing('alternating'),pollen=routing('clustered'))
        for s in range(50):
            pool, ledger=production(INIT,par,Random(s))
            self.assertEqual(sum(map(len,pool)),sum(ledger['source_viable_seeds']))
            self.assertEqual(ledger['source_attempts'],(3,)*6)
            self.assertLessEqual(sum(ledger['source_exported_seeds']),sum(ledger['source_viable_seeds']))
            self.assertLessEqual(sum(ledger['source_joint_seeds']),sum(ledger['source_viable_seeds']))

    def test_selection_is_inactive_when_no_recruitment_overcrowding(self):
        par=Network(seed=routing('alternating'),pollen=routing('clustered'))
        for s in range(45):
            pools,ledger=production(INIT,par,Random(s))
            x=tuple(settle(p,par.capacity,1/3,Random(s)) for p in pools)
            y=tuple(settle(p,par.capacity,3.,Random(s)) for p in pools)
            for i,p in enumerate(pools):
                if len(p)<=par.capacity:
                    self.assertEqual(x[i],y[i])

    def test_exact_first_gen_production_matches_independent_mc(self):
        for topology in ('clustered','alternating'):
            par=Network(pollen=routing(topology),seed=routing('alternating'))
            exact=_first_hybrid_exact(par)['probability_any_hybrid_seed_produced']
            n=1200
            seen=sum(sum(production(INIT,par,Random(i*191+100))[1]['source_joint_seeds'])>0 for i in range(n))
            ph=seen/n
            se=math.sqrt(exact*(1-exact)/n)
            self.assertLess(abs(ph-exact),4.2*se+.004)

    def test_stage_invariance_across_selection_with_same_replicate_seed(self):
        par=Network(seed=routing('clustered'),pollen=routing('clustered'))
        for key in range(18):
            values=[trajectory(par,w,horizon=1,seed=1000+key) for w in (0.,1/3,1.,3.)]
            self.assertEqual(len({v['gen1_joint_seed_count'] for v in values}),1)
            self.assertEqual(len({v['gen1_joint_seed_produced'] for v in values}),1)
            self.assertEqual(len({v['total_viable_seeds'] for v in values}),1)

    def test_fixed_pool_forecast_value_is_not_historical_memory(self):
        a=exact_three_seed_witness(1/3)['conditional_joint_recruitment']
        b=exact_three_seed_witness(3.)['conditional_joint_recruitment']
        pooled=(a+b)/2
        coarse=pooled*(1-pooled)
        complete=(a*(1-a)+b*(1-b))/2
        self.assertGreater(coarse-complete,0)
        self.assertAlmostEqual(coarse-complete,(b-a)**2/4.)

    def test_paired_mc_interval_is_centered_on_signed_difference(self):
        for data in ([0,0,0],[-1,0,1],[1,1,1],[-1,-1,0]):
            row=_paired_summary(data)
            self.assertAlmostEqual(row['difference'],sum(data)/len(data))
            self.assertLessEqual(row['paired_mc_95'][0],row['difference'])
            self.assertGreaterEqual(row['paired_mc_95'][1],row['difference'])
        with self.assertRaises(ValueError):_paired_summary([1])

    def test_predeclared_cases_keep_demographic_and_resource_controls_fixed(self):
        c=cases()
        self.assertEqual(len(c),12)
        self.assertEqual(len({p.capacity for p,w in c.values()}),1)
        self.assertEqual(len({p.seed_attempts for p,w in c.values()}),1)
        self.assertEqual(len({p.q for p,w in c.values()}),1)
        self.assertEqual(set(w for p,w in c.values()),{0.,1/3,1.,3.})

    def test_model_benchmark_is_seed_reproducible_and_falsifiable(self):
        x=benchmark(replicates=50,horizon=3,seed=777)
        y=benchmark(replicates=50,horizon=3,seed=777)
        self.assertEqual(x,y)
        self.assertEqual(x['design']['number_cells'],12)
        self.assertEqual(x['exact_first_gen_joint_seed_production']['components']['probability_any_hybrid_seed_produced'],0)
        self.assertEqual(x['scenarios']['clustered_lethal']['probabilities']['ever_joint_adult_recruited'],0)
        self.assertEqual(x['scenarios']['alternating_lethal']['probabilities']['final_joint_adult'],0)
        self.assertEqual(x['paired_mc_contrasts']['neutral_alt_minus_cluster']['gen1_joint_seed_produced']['difference'],
                         x['paired_mc_contrasts']['favored_alt_minus_cluster']['gen1_joint_seed_produced']['difference'])


if __name__=='__main__':
    unittest.main()
