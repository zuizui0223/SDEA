"""v2 exact-state normalization, non-equivalence and adverse controls."""
from dataclasses import replace
from math import isclose
import unittest

from sdea.multigeneration import Dynamics, evolve_exact, history_masks, patch_alignment, run_multigeneration, _transition, horizon_brier_risks, multigen_cost_crossover
from scripts.verify_fixture import assert_same


class TestMultigeneration(unittest.TestCase):
    def test_mass_conservation_each_generation(self):
        p = Dynamics(module_retention=.82,module_reacquisition=.10,mutation=.03,
                     immigrant_fraction=.04,post_expression_damage=.15,
                     expression_recoupling=.12)
        trajectory = evolve_exact(history_masks("odd"),2,"high",8,p)["trajectory"]
        for row in trajectory:
            self.assertAlmostEqual(row["alive"]+row["extinct"],1.,places=11)
            for key,value in row.items():
                if key!="generation":
                    self.assertTrue(-1e-12 <= value <= 1+1e-12,(key,value))

    def test_allele_transition_is_normalized(self):
        p = Dynamics(module_retention=.96,module_reacquisition=.08,
                     immigrant_fraction=.12,post_expression_damage=.3)
        for m in (0,3,7):
            for copies in range(2*p.diploids+1):
                for qi in (0,1):
                    for expression in (False,True):
                        self.assertAlmostEqual(sum(w for _,w in _transition(m,copies,qi,expression,p)),1.,places=12)

    def test_first_is_not_three_or_alive(self):
        f = evolve_exact({7:1.}, horizon=8)["final"]
        self.assertGreater(f["first_expression_ever"],f["three_consecutive_survived_ever"])
        self.assertGreaterEqual(f["three_consecutive_survived_ever"],f["alive_and_ever_three"])
        self.assertGreaterEqual(f["alive_and_ever_three"],f["alive_and_current_three"])

    def test_short_horizon_cannot_produce_three_survived(self):
        for h in (1,2):
            f=evolve_exact({7:1.},horizon=h)["final"]
            self.assertEqual(f["three_consecutive_survived_ever"],0.)
            self.assertGreater(f["first_expression_ever"],0.)

    def test_history_even_closed_odd_open_when_no_module_input(self):
        a=evolve_exact(history_masks("even"),horizon=8)["final"]
        b=evolve_exact(history_masks("odd"),horizon=8)["final"]
        self.assertEqual(a["first_expression_ever"],0.)
        self.assertGreater(b["three_consecutive_survived_ever"],0.)

    def test_reacquisition_breaks_absolute_module_closure(self):
        p=Dynamics(module_reacquisition=.2)
        self.assertEqual(evolve_exact({0:1.},horizon=8)["final"]["first_expression_ever"],0.)
        self.assertGreater(evolve_exact({0:1.},horizon=8,par=p)["final"]["first_expression_ever"],0.)

    def test_allele_absence_closed_and_mutation_reopens(self):
        x=evolve_exact({7:1.},initial_copies=0,horizon=8)["final"]
        y=evolve_exact({7:1.},initial_copies=0,horizon=8,par=Dynamics(mutation=.04))["final"]
        self.assertEqual(x["first_expression_ever"],0.)
        self.assertEqual(x["alive_and_allele_supported"],0.)
        self.assertGreater(y["first_expression_ever"],0.)

    def test_certain_post_expression_damage_first_still_possible_run_not(self):
        p=Dynamics(post_expression_damage=1.,module_reacquisition=0.)
        f=evolve_exact({7:1.},horizon=8,par=p)["final"]
        self.assertGreater(f["first_expression_ever"],0.)
        self.assertEqual(f["three_consecutive_survived_ever"],0.)

    def test_alignment_and_negative_control(self):
        a=patch_alignment(horizon=7)
        self.assertGreater(a["aligned"],a["misaligned"])
        equal=patch_alignment(horizon=7,par=Dynamics(q_low=.6,q_high=.6))
        self.assertAlmostEqual(equal["difference"],0.,places=13)

    def test_alleles_expected_to_be_neutral_given_survival(self):
        p=Dynamics(selection_strength=0.,assurance_cost=0.,mutation=0.,immigrant_fraction=0.)
        events=_transition(7,2,1,False,p)
        mean=sum(k*weight for (m,k,q),weight in events)
        self.assertAlmostEqual(mean,2.,places=12)

    def test_survival_increases_with_assurance_cost_zero(self):
        lo=evolve_exact({7:1.},horizon=8,par=Dynamics(assurance=.1,assurance_cost=0.))["final"]
        hi=evolve_exact({7:1.},horizon=8,par=Dynamics(assurance=.9,assurance_cost=0.))["final"]
        self.assertGreater(hi["alive"],lo["alive"])

    def test_time_monotonic_cumulative_but_not_occupancy(self):
        rows=evolve_exact({7:1.},horizon=10)["trajectory"]
        for a,b in zip(rows,rows[1:]):
            self.assertLessEqual(a["first_expression_ever"],b["first_expression_ever"]+1e-12)
            self.assertLessEqual(a["three_consecutive_survived_ever"],b["three_consecutive_survived_ever"]+1e-12)
            self.assertGreaterEqual(a["alive"]+1e-12,b["alive"])

    def test_fail_closed_parameters(self):
        with self.assertRaises(ValueError): Dynamics(diploids=0)
        with self.assertRaises(ValueError): Dynamics(module_reacquisition=1.1)
        with self.assertRaises(ValueError): Dynamics(q_low=.9,q_high=.4)
        with self.assertRaises(ValueError): Dynamics(consecutive_required=1)
        with self.assertRaises(ValueError): evolve_exact({9:1.})
        with self.assertRaises(ValueError): evolve_exact({7:.7})
        with self.assertRaises(ValueError): evolve_exact({7:1.},initial_copies=7)
        with self.assertRaises(ValueError): evolve_exact({7:1.},start_q="middle")

    def test_history_proxy_redundant_given_full_state(self):
        risks=horizon_brier_risks(horizon=8)
        self.assertLess(risks["history"],risks["visible"])
        self.assertLess(risks["assay_history"],risks["assay"])
        self.assertAlmostEqual(risks["full_history"],risks["full"],places=13)
        null=horizon_brier_risks(null_history=True)
        self.assertAlmostEqual(null["history"],null["visible"],places=13)
        perfect=horizon_brier_risks(sensitivity=1.,specificity=1.)
        self.assertAlmostEqual(perfect["assay"],perfect["full"],places=13)

    def test_crossover_endpoint_depends_on_measure(self):
        a=multigen_cost_crossover("alive_and_ever_three")
        b=multigen_cost_crossover("three_consecutive_survived_ever")
        c=multigen_cost_crossover("alive_and_joint_toolkit_supported")
        for x in (a,b,c):
            self.assertLess(abs(x["root_residual"]),1e-10)
            self.assertGreater(x["diff_at_cost_zero"],0.)
            self.assertLess(x["diff_at_upper"],0.)
        self.assertNotAlmostEqual(a["cost_star"],b["cost_star"],places=3)
        with self.assertRaises(ValueError):
            multigen_cost_crossover("alive")
        with self.assertRaises(ValueError):
            multigen_cost_crossover(par=Dynamics(mutation=.1))

    def test_expression_recoupling_is_nonnegative_in_this_example(self):
        base=evolve_exact({7:1.},horizon=8)["final"]
        feedback=evolve_exact({7:1.},horizon=8,
             par=Dynamics(expression_recoupling=.10))["final"]
        self.assertGreater(feedback["three_consecutive_survived_ever"],
                           base["three_consecutive_survived_ever"])

    def test_fixture_recomputation(self):
        import json
        from pathlib import Path
        path=Path(__file__).resolve().parents[1]/"results/v2_multigeneration_exact.json"
        expected=json.loads(path.read_text())
        assert_same(expected,run_multigeneration())


if __name__=="__main__": unittest.main()
