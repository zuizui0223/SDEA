"""SDEA v18 source-grain, negative-gate and retrospective risk audits."""
import copy
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sdea import evidence_ladder as e


class EvidenceLadderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sources = e.load_locked()
        cls.out = e.build_ladder(cls.sources)

    def test_original_unit_and_historical_source_count(self):
        self.assertEqual(self.out["independent_original_studies"], 3)
        self.assertEqual(self.out["comparisons"]["ulex_flower_traits_over_elevation"]["original_heldout_records"], 202)
        self.assertEqual(self.out["comparisons"]["sinapis_visits_over_flower_supply"]["original_heldout_records"], 52)
        self.assertEqual(self.out["adversarial_controls"]["Ipomopsis_context"]["years_total"], 9)
        self.assertEqual(self.out["joint_G_interaction_phenotype_fitness_same_individual"], "NOT_IDENTIFIABLE")

    def test_source_locked_quantitative_anchors(self):
        c = self.out["comparisons"]
        self.assertAlmostEqual(c["ulex_flower_traits_over_elevation"]["equal_block_mse_gain"], .0010207813409493843, places=12)
        self.assertAlmostEqual(c["ulex_flower_traits_over_elevation"]["relative_gain_ratio_of_equal_block_means"], .11687039006760376, places=12)
        self.assertEqual(c["ulex_flower_traits_over_elevation"]["blocks_with_positive_gain"], 5)
        self.assertGreater(c["ulex_flower_traits_over_elevation"]["min_omit_one_block_gain"], 0)
        self.assertAlmostEqual(c["ulex_genomic_similarity_over_flower_traits"]["equal_block_mse_gain"], .00001722707955827491, places=12)
        self.assertEqual(c["ulex_genomic_similarity_over_flower_traits"]["blocks_with_positive_gain"], 4)
        self.assertAlmostEqual(c["sinapis_flower_shape_over_flower_supply"]["equal_block_mse_gain"], -.023034472099674543, places=12)
        self.assertEqual(c["sinapis_flower_shape_over_flower_supply"]["blocks_with_positive_gain"], 1)
        self.assertAlmostEqual(c["sinapis_visits_over_flower_supply"]["equal_block_mse_gain"], .058002058310772577, places=12)
        self.assertAlmostEqual(c["sinapis_visits_over_flower_supply"]["relative_gain_ratio_of_equal_block_means"], .24596901754624959, places=12)
        self.assertEqual(c["sinapis_visits_over_flower_supply"]["blocks_with_positive_gain"], 5)
        self.assertGreater(c["sinapis_visits_over_flower_supply"]["min_omit_one_block_gain"], 0)
        self.assertEqual(c["sinapis_visits_over_supply_and_shape"]["blocks_with_positive_gain"], 4)

    def test_no_pooled_meta_or_genetic_reexpression_claim(self):
        self.assertTrue(self.out["cross_study_numeric_effect_pooling"].startswith("PROHIBITED"))
        self.assertEqual(self.out["trait_regain_across_generations"], "NOT_MEASURED")
        self.assertFalse(self.out["adversarial_controls"]["Ulex_genomic_identity"]["individual_genome_identity_gate_passed"])
        self.assertTrue(self.out["adversarial_controls"]["Sinapis_visitor_identity"]["descriptive_within_source_gate_passed"])
        self.assertFalse(self.out["adversarial_controls"]["Ipomopsis_context"]["strict_within_context_gate_passed"])

    def test_duplicate_source_block_fails(self):
        d = copy.deepcopy(self.sources)
        d["ulex"]["heldout_sites"][1]["site"] = d["ulex"]["heldout_sites"][0]["site"]
        with self.assertRaisesRegex(e.EvidenceGateError, "DUPLICATE_SOURCE_BLOCK"):
            e.build_ladder(d)

    def test_missing_original_individual_fails(self):
        d = copy.deepcopy(self.sources)
        d["sinapis"]["heldout_results"]["folds"][0]["heldout_plants"] -= 1
        with self.assertRaisesRegex(e.EvidenceGateError, "SOURCE_RECORD_DENOMINATOR_DRIFT"):
            e.build_ladder(d)

    def test_nonfinite_or_negative_risk_fails(self):
        for wrong in (float("nan"), -1.0):
            d = copy.deepcopy(self.sources)
            d["sinapis"]["heldout_results"]["folds"][0]["mse_models"]["M3"] = wrong
            with self.assertRaisesRegex(e.EvidenceGateError, "INVALID_PREDICTION_RISK"):
                e.build_ladder(d)

    def test_changed_original_result_or_null_fails(self):
        d = copy.deepcopy(self.sources)
        d["ulex"]["heldout_sites"][0]["plus_relatedness_residual_mse"] += .001
        with self.assertRaisesRegex(e.EvidenceGateError, "RECEIPT_DRIFT"):
            e.build_ladder(d)
        d = copy.deepcopy(self.sources)
        d["ulex_placebo"]["within_site_genotype_identity_placebo"]["equal_site_genomic_gain"]["q975"] = -1
        with self.assertRaisesRegex(e.EvidenceGateError, "ULEX_GENOMIC_PLACEBO_DRIFT"):
            e.build_ladder(d)

    def test_same_original_study_cannot_be_counted_twice(self):
        d = copy.deepcopy(self.sources)
        d["ipomopsis"]["v13_previous_same_source_record"]["seed_initiation_log1p"]["v13_row_weighted_gain"] = -1
        with self.assertRaisesRegex(e.EvidenceGateError, "IPOMOPSIS_CONTEXT_CHALLENGE_DRIFT"):
            e.build_ladder(d)

    def test_missing_genetic_or_visitor_link_must_not_be_promoted(self):
        d = copy.deepcopy(self.sources)
        d["ulex"]["audit"]["pollinator_census_individual_link"] = True
        with self.assertRaisesRegex(e.EvidenceGateError, "ORIGINAL_UNIT_OR_OBSERVABILITY_DRIFT"):
            e.build_ladder(d)
        d = copy.deepcopy(self.sources)
        d["sinapis"]["audit"]["plant_genetic_alleles_or_ancestors_measured"] = True
        with self.assertRaisesRegex(e.EvidenceGateError, "ORIGINAL_UNIT_OR_OBSERVABILITY_DRIFT"):
            e.build_ladder(d)

    def test_omission_does_not_equal_independent_refitting(self):
        c = self.out["comparisons"]["sinapis_visits_over_flower_supply"]
        self.assertGreater(c["min_omit_one_block_gain"], 0)
        self.assertIn("not independent re-fitting", self.out["inference_limit"])


if __name__ == "__main__":
    unittest.main()
