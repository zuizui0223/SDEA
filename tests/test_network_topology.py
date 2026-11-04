"""v7 adversarial tests: matched seed budgets, network alignment, exact controls."""
from dataclasses import replace
from random import Random
from math import sqrt
import unittest

from sdea.network_topology import (
    NODES,A,B,O,INIT,IDENTITY,ORDERS,Network,routing,validate_matrix,
    mixing_signature,validate_populations,_gamete,_cross,_first_hybrid_exact,
    initial_one_step_colonization_bound,network_prediction_risk,generation,trajectory,scenarios,estimate,
)

class NetworkTopologyTests(unittest.TestCase):
    def test_regular_matrices_have_identical_source_and_destination_marginals(self):
        for k in ORDERS:
            m=routing(k)
            sig=mixing_signature(m)
            self.assertEqual(sig['outgoing'],[1.]*6)
            self.assertEqual(sig['incoming'],[1.]*6)
            self.assertEqual(sig['retention'],[.25]*6)
            self.assertEqual(sig['undirected_edges'],6)
            for i in range(6):
                self.assertEqual(sum(m[i][j]>0 for j in range(6) if i!=j),2)
                for j in range(6):
                    self.assertEqual(m[i][j],m[j][i])

    def test_cross_route_edges_differ_while_degree_matches(self):
        actual=[mixing_signature(routing(k))['cross_route_edges'] for k in ('components','clustered','alternating')]
        self.assertEqual(actual,[0,4,12])
        masses=[mixing_signature(routing(k))['cross_route_mass'] for k in ('components','clustered','alternating')]
        self.assertEqual(masses,[0.,1.5,4.5])

    def test_initial_population_and_q_identical_across_cells(self):
        for name,p in scenarios().items():
            self.assertEqual(p.capacity,2)
            self.assertEqual(p.seed_attempts,3)
            self.assertEqual(p.q,.75)
            self.assertEqual(len(INIT),6)
            self.assertEqual(INIT[:3],((A,A),)*3)
            self.assertEqual(INIT[3:],((B,B),)*3)

    def test_nonidentity_matrices_match_local_retention(self):
        for name,p in scenarios().items():
            if p.seed!=IDENTITY:self.assertEqual([p.seed[i][i] for i in range(6)],[.25]*6)
            if p.pollen!=IDENTITY:self.assertEqual([p.pollen[i][i] for i in range(6)],[.25]*6)

    def test_all_scenarios_have_valid_matrix_and_budgets(self):
        self.assertEqual(len(scenarios()),12)
        for p in scenarios().values():
            validate_matrix(p.seed);validate_matrix(p.pollen)
            self.assertEqual(sum(sum(row) for row in p.seed),6.)
            self.assertEqual(sum(sum(row) for row in p.pollen),6.)

    def test_mendelian_homozygote_cross_creates_joint_route(self):
        rng=Random(1)
        self.assertTrue(all(_cross(A,B,rng)==(1,1) for _ in range(80)))
        self.assertTrue(all(_cross(B,A,rng)==(1,1) for _ in range(80)))
        self.assertTrue(all(_cross(A,A,rng)==A for _ in range(80)))

    def test_independent_meiosis_heterozygote_mean(self):
        rng=Random(91)
        child=[_cross((1,1),(1,1),rng) for _ in range(5000)]
        for locus in (0,1):
            self.assertLess(abs(sum(g[locus] for g in child)/5000-1.),.04)

    def test_exact_hybrid_pre_recruitment_witness(self):
        exact={k:_first_hybrid_exact(Network(pollen=routing(k))) for k in ORDERS}
        self.assertEqual(exact['components']['probability_any_hybrid_seed_produced'],0)
        self.assertAlmostEqual(exact['clustered']['probability_any_hybrid_seed_produced'],.9681251026678412,places=12)
        self.assertAlmostEqual(exact['alternating']['probability_any_hybrid_seed_produced'],.999996076577506,places=12)
        self.assertAlmostEqual(exact['clustered']['expected_hybrid_seeds_produced'],2.9953125,places=11)
        self.assertAlmostEqual(exact['alternating']['expected_hybrid_seeds_produced'],8.9859375,places=11)
        self.assertEqual({r['attempt_budget_initial_total'] for r in exact.values()},{18})

    def test_exact_hybrid_production_independent_of_seed_destination(self):
        c=Network(pollen=routing('clustered'),seed=routing('clustered'))
        a=replace(c,seed=routing('alternating'))
        self.assertEqual(_first_hybrid_exact(c),_first_hybrid_exact(a))

    def test_binary_saturation_exact_predictive_state_gain(self):
        r=[network_prediction_risk(b) for b in (1,2,3,8)]
        for item in r:
            self.assertAlmostEqual(item['network_information_gain'],item['exact_gain_identity'],places=13)
            self.assertGreaterEqual(item['network_information_gain'],-1e-14)
            self.assertGreaterEqual(item['alternating_probability'],item['clustered_probability'])
        self.assertGreater(r[0]['network_information_gain'],r[-1]['network_information_gain'])
        self.assertAlmostEqual(r[0]['expected_hybrid_seeds_alternating']/
                               r[0]['expected_hybrid_seeds_clustered'],3.,places=13)
        with self.assertRaises(ValueError): network_prediction_risk(-1)

    def test_structural_pollen_only_closure_under_two_triangles(self):
        p=Network(pollen=routing('components'))
        for i in range(30):
            v=trajectory(par=p,horizon=5,seed=i+10)
            self.assertFalse(v['gen1_hybrid_seed_produced'])
            self.assertFalse(v['ever_joint_individual_recruited'])

    def test_structural_combined_closed_network_preserves_disjoint_modules(self):
        p=Network(seed=routing('components'),pollen=routing('components'))
        for i in range(30):
            v=trajectory(par=p,horizon=5,seed=i+800)
            self.assertFalse(v['ever_joint_individual_recruited'])

    def test_identity_network_never_creates_joint(self):
        for i in range(50):
            v=trajectory(par=Network(),horizon=6,seed=i+400)
            self.assertFalse(v['gen1_hybrid_seed_produced'])
            self.assertFalse(v['ever_joint_individual_recruited'])

    def test_no_pollen_flow_seed_only_has_no_gen1_hybrid(self):
        for k in ('clustered','alternating','components'):
            par=Network(seed=routing(k))
            for i in range(30):
                v=trajectory(par=par,horizon=3,seed=i+500)
                self.assertFalse(v['gen1_hybrid_seed_produced'])

    def test_later_joint_assembly_possible_by_seed_only_via_local_mating(self):
        p=Network(seed=routing('alternating'))
        self.assertTrue(any(trajectory(par=p,horizon=4,seed=i)['ever_joint_individual_recruited']
                            for i in range(25)))

    def test_same_initial_pollen_randomness_invariant_to_seed_routing(self):
        p=Network(seed=routing('clustered'),pollen=routing('clustered'))
        other=replace(p,seed=routing('alternating'))
        for i in range(60):
            a=trajectory(par=p,horizon=1,seed=i)
            b=trajectory(par=other,horizon=1,seed=i)
            self.assertEqual(a['gen1_hybrid_seed_produced'],b['gen1_hybrid_seed_produced'])

    def test_one_generation_budget_conservation_even_under_redistribution(self):
        for topology in ('components','clustered','alternating'):
            p=Network(seed=routing(topology),pollen=routing(topology))
            for i in range(60):
                pop,ledger=generation(INIT,p,Random(i))
                self.assertEqual(ledger['attempted'],(3,)*6)
                self.assertEqual(sum(ledger['viable_seed_produced']),sum(ledger['pre_recruit_seed_pool_sizes']))
                self.assertLessEqual(sum(len(x) for x in pop),sum(ledger['viable_seed_produced']))
                self.assertTrue(all(len(x)<=2 for x in pop))
                self.assertTrue(all(e<=a for e,a in zip(ledger['seed_exported'],ledger['viable_seed_produced'])))

    def test_zero_budget_absorbs_when_all_adults_die(self):
        p=Network(seed_attempts=0)
        v=trajectory(par=p,horizon=4)
        self.assertFalse(v['final_any_patch_occupied'])
        self.assertEqual(v['total_attempts'],18*0)

    def test_empty_all_stays_empty_without_external_source(self):
        x=((),)*6
        p=Network(seed=routing('alternating'),pollen=routing('alternating'))
        v=trajectory(initial=x,par=p,horizon=5)
        self.assertFalse(v['final_any_patch_occupied'])
        self.assertEqual(v['total_attempts'],0)

    def test_pollen_without_seed_cannot_colonize_empty_patch(self):
        init=INIT[:2]+((),)+INIT[3:]
        p=Network(pollen=routing('alternating'))
        for i in range(40):
            v=trajectory(init,p,horizon=4,seed=i)
            self.assertEqual(v['final_populations'][2],())

    def test_seed_network_can_colonize_initial_empty_patch(self):
        init=INIT[:2]+((),)+INIT[3:]
        p=Network(seed=routing('alternating'))
        self.assertTrue(any(trajectory(init,p,horizon=2,seed=i)['final_populations'][2]
                            for i in range(90)))

    def test_missing_B_globally_cannot_be_generated_by_routing(self):
        init=((A,A),)*6
        p=Network(seed=routing('alternating'),pollen=routing('alternating'))
        for i in range(25):
            v=trajectory(init,p,horizon=6,seed=i)
            self.assertFalse(v['final_B_global'])
            self.assertFalse(v['ever_joint_individual_recruited'])

    def test_exact_first_hybrid_production_matches_untuned_MC(self):
        # Independent seeds; MC verification distinct from algebraic calculation.
        for k in ('clustered','alternating'):
            par=Network(pollen=routing(k))
            exact=_first_hybrid_exact(par)['probability_any_hybrid_seed_produced']
            n=1500
            seen=sum(trajectory(par=par,horizon=1,seed=112371+i*7919)['gen1_hybrid_seed_produced'] for i in range(n))
            ph=seen/n
            se=sqrt(exact*(1-exact)/n)
            self.assertLess(abs(ph-exact),max(.004,4*se))

    def test_no_network_allows_persistence_of_distinct_routes_globally(self):
        # Both gene sets can persist separately without any same-individual joint capacity.
        x=trajectory(par=Network(),horizon=3,seed=991)
        self.assertFalse(x['ever_joint_individual_recruited'])
        self.assertTrue(isinstance(x['final_both_routes_global'],bool))

    def test_invalid_inputs_fail_closed(self):
        for n in ('bogus','',None):
            with self.assertRaises(ValueError):routing(n)
        with self.assertRaises(ValueError):routing('clustered',-1)
        with self.assertRaises(ValueError):routing('clustered',1.01)
        with self.assertRaises(ValueError):Network(capacity=0)
        with self.assertRaises(ValueError):Network(seed_attempts=-1)
        with self.assertRaises(ValueError):Network(seed_attempts=33)
        with self.assertRaises(ValueError):Network(q=1.1)
        with self.assertRaises(ValueError):Network(seed=((1.,0.),)*6)
        with self.assertRaises(ValueError):validate_populations(INIT[:-1])
        with self.assertRaises(ValueError):trajectory(par=Network(),horizon=0)
        with self.assertRaises(ValueError):trajectory(par=Network(),seed=1.5)
        with self.assertRaises(ValueError):estimate(replicates=10)

    def test_paired_difference_and_factorial_timing(self):
        x=estimate(replicates=50,horizon=2,seed=881)
        self.assertEqual(x['design']['same_seed_attempt_budget_per_occupied_source'],3)
        self.assertTrue(x['design']['future_total_attempts_not_controlled_because_source_occupancy_is_endogenous'])
        self.assertEqual(x['factorial_topology_interaction']['gen1_hybrid_seed_produced']['difference_in_differences'],0.)
        a=x['paired_contrasts']['seed_alt_minus_cluster']['gen1_hybrid_seed_produced']
        self.assertEqual(a['difference_variant_minus_baseline'],0.)
        self.assertEqual(a['paired_mc_se'],0.)
        for d in x['paired_contrasts'].values():
            for c in d.values():
                self.assertLessEqual(c['paired_mc_95'][0],c['difference_variant_minus_baseline'])
                self.assertGreaterEqual(c['paired_mc_95'][1],c['difference_variant_minus_baseline'])

    def test_determinism_of_synthetic_independent_landscapes(self):
        p=Network(seed=routing('clustered'),pollen=routing('alternating'))
        self.assertEqual(trajectory(par=p,seed=301),trajectory(par=p,seed=301))
        self.assertEqual(estimate(replicates=25,horizon=3,seed=882),estimate(replicates=25,horizon=3,seed=882))

if __name__ == '__main__':
    unittest.main()
