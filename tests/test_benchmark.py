"""Deterministic, negative-control and inference-boundary checks for SDEA v1."""

from dataclasses import replace
import math
import unittest

from scripts.verify_fixture import assert_same

from sdea.benchmark import (
    Parameters, HISTORIES, aligned_patch_test, assay_likelihood, assurance_contrast,
    assurance_cost_crossover,
    brier_risks, module_marginals, outcome_probability, run_benchmarks,
    transmitted_frequency,
)


class TestSDEA(unittest.TestCase):
    def test_matched_single_and_pairwise_but_different_joint(self):
        e, o = module_marginals("even"), module_marginals("odd")
        self.assertEqual(e["individual"], [0.5]*3)
        self.assertEqual(e["pairs"], [0.25]*3)
        self.assertEqual(e["individual"], o["individual"])
        self.assertEqual(e["pairs"], o["pairs"])
        self.assertEqual((e["all_three"], o["all_three"]), (0.0, .25))

    def test_history_label_informative_when_toolkit_is_hidden(self):
        x = brier_risks()
        self.assertLess(x["history"], x["visible"])
        self.assertLess(x["assay_history"], x["assay"])
        self.assertLess(x["full"], x["history"])
        self.assertAlmostEqual(x["full"], x["full_history"], places=13)

    def test_uninformative_history_negative_control(self):
        x = brier_risks(null_history=True)
        self.assertAlmostEqual(x["visible"], x["history"], places=13)
        self.assertAlmostEqual(x["assay"], x["assay_history"], places=13)

    def test_imperfect_and_perfect_assay(self):
        none = brier_risks(sensitivity=.5, specificity=.5)
        perfect = brier_risks(sensitivity=1., specificity=1.)
        self.assertAlmostEqual(none["assay"], none["visible"], places=13)
        self.assertAlmostEqual(perfect["assay"], perfect["full"], places=13)
        self.assertAlmostEqual(perfect["assay_history"], perfect["assay"], places=13)

    def test_assay_normalization(self):
        for mask in range(8):
            self.assertAlmostEqual(sum(assay_likelihood(z, mask, .9, .95) for z in range(8)), 1.)

    def test_same_coarse_marginals_different_transition(self):
        x = aligned_patch_test()
        self.assertGreater(x["aligned"], x["misaligned"])
        self.assertGreater(x["difference"], 0)
        control = aligned_patch_test(.6, .6)
        self.assertAlmostEqual(control["difference"], 0., places=13)

    def test_survival_and_support_can_diverge(self):
        x = assurance_contrast(q=.8, cost=2.)["difference_high_minus_low"]
        self.assertGreater(x["survival"], 0)
        self.assertLess(x["genetic_support"], 0)
        self.assertLess(x["joint"], 0)

    def test_cost_free_assurance_control(self):
        for q in (.2, .5, .8):
            x = assurance_contrast(q, 0.)["difference_high_minus_low"]
            self.assertGreater(x["survival"], 0.)
            self.assertAlmostEqual(x["genetic_support"], 0., places=13)
            self.assertGreater(x["joint"], 0.)

    def test_allele_support_closed_without_input_and_reopened(self):
        p = replace(Parameters(), allele_frequency=0.)
        self.assertEqual(outcome_probability(7, .6, p, "genetic_support"), 0.)
        self.assertGreater(outcome_probability(7, .6, replace(p, mutation_rate=.02), "genetic_support"), 0.)
        self.assertGreater(outcome_probability(7, .6, replace(p, immigrant_fraction=.02), "genetic_support"), 0.)

    def test_expression_eligibility_not_survival(self):
        p = Parameters()
        self.assertEqual(outcome_probability(0, .8, p, "expression"), 0)
        self.assertGreater(outcome_probability(0, .8, p, "survival"), 0)
        self.assertGreater(outcome_probability(7, .8, p, "expression"), 0)
        self.assertLess(outcome_probability(7, .8, p, "joint"), outcome_probability(7, .8, p, "expression"))

    def test_finite_sample_support_probability(self):
        p = Parameters(allele_frequency=.2, assurance_cost=0., selection_strength=0.)
        for n in (1, 2, 10):
            k = replace(p, offspring_diploids=n)
            self.assertAlmostEqual(outcome_probability(7, .5, k, "genetic_support"), 1 - .8**(2*n))

    def test_endpoints_probability_bounds(self):
        p = Parameters()
        for mask in HISTORIES["even"] + HISTORIES["odd"]:
            for q in (0., .01, .5, 1.):
                for endpoint in ("expression", "survival", "genetic_support", "joint"):
                    v = outcome_probability(mask, q, p, endpoint)
                    self.assertTrue(math.isfinite(v) and 0 <= v <= 1)

    def test_invalid_inputs_fail_closed(self):
        with self.assertRaises(ValueError):
            Parameters(assurance=1.1)
        with self.assertRaises(ValueError):
            Parameters(offspring_diploids=0)
        with self.assertRaises(ValueError):
            Parameters(assurance_cost=-1.)
        with self.assertRaises(ValueError):
            outcome_probability(8, .5)
        with self.assertRaises(ValueError):
            outcome_probability(7, -0.2)
        with self.assertRaises(ValueError):
            outcome_probability(7, .5, endpoint="true_establishment")

    def test_assurance_cost_frontier_crosses(self):
        roots = [assurance_cost_crossover(q)["cost_star"] for q in (.2, .5, .8)]
        self.assertGreater(roots[0], roots[1])
        self.assertGreater(roots[1], roots[2])
        for q, root in zip((.2, .5, .8), roots):
            self.assertGreater(assurance_contrast(q, root - .02)["difference_high_minus_low"]["joint"], 0)
            self.assertLess(assurance_contrast(q, root + .02)["difference_high_minus_low"]["joint"], 0)
            self.assertAlmostEqual(assurance_contrast(q, root)["difference_high_minus_low"]["joint"], 0., places=12)

    def test_cost_frontier_domain_and_input_bounds(self):
        with self.assertRaises(ValueError):
            assurance_cost_crossover(0.)
        with self.assertRaises(ValueError):
            assurance_cost_crossover(1.)
        with self.assertRaises(ValueError):
            assurance_cost_crossover(.5, par=Parameters(mutation_rate=.02))
        with self.assertRaises(ValueError):
            assurance_cost_crossover(.5, par=Parameters(allele_frequency=0.))
        with self.assertRaises(ValueError):
            assurance_cost_crossover(.5, upper=.2)

    def test_fixture_allows_machine_roundoff(self):
        assert_same({"score": [0.1, {"p": 0.2}]}, {"score": [0.1 + 1e-15, {"p": 0.2}]})
        with self.assertRaises(AssertionError):
            assert_same({"score": 0.1}, {"score": 0.10001})

    def test_fixture_rejects_structure_changes(self):
        with self.assertRaises(AssertionError):
            assert_same({"a": 2}, {"b": 2})
        with self.assertRaises(AssertionError):
            assert_same([True], [1])
        with self.assertRaises(AssertionError):
            assert_same([1, 2], [1])

    def test_deterministic_protocol(self):
        self.assertEqual(run_benchmarks(), run_benchmarks())
        self.assertEqual(len(run_benchmarks()["assurance_tradeoff"]), 12)
        self.assertEqual(len(run_benchmarks()["conditional_cost_frontier"]), 3)


if __name__ == "__main__":
    unittest.main()
