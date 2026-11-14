"""Read-only, source-locked cross-study predictive evidence ladder (SDEA v18).

This is a retrospective audit of previously frozen aggregate result receipts.
There is no source-row access, biological refit, effect-size pooling or meta-analysis.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FILES = {
    "ulex": "results/v15_ulex_partial.json",
    "ulex_placebo": "results/v16_ulex_genomic_placebo.json",
    "sinapis": "results/v17_sinapis_process.json",
    "ipomopsis": "results/v14_within_context.json",
}


class EvidenceGateError(ValueError):
    """Malformed or scientifically inadmissible source-level comparison."""


def load_locked(root=ROOT):
    data = {}
    for key, path in FILES.items():
        p = Path(root) / path
        data[key] = json.loads(p.read_text(encoding="utf-8"))
    return data


def _assert_close(a, b, what, tol=1e-11):
    if not (math.isfinite(a) and math.isfinite(b) and abs(a - b) <= tol):
        raise EvidenceGateError(f"RECEIPT_DRIFT:{what}")


def summarize_folds(folds, *, key, count_key, baseline_key, model_key, expected_n):
    """Equal-source-unit and original-record weighted performance, not independent re-fits."""
    if len(folds) < 3:
        raise EvidenceGateError("TOO_FEW_SOURCE_BLOCKS")
    seen = set()
    rows = []
    for f in folds:
        name = str(f[key])
        if name in seen:
            raise EvidenceGateError("DUPLICATE_SOURCE_BLOCK")
        seen.add(name)
        n = f[count_key]
        a, b = f[baseline_key], f[model_key]
        if not isinstance(n, int) or isinstance(n, bool) or n <= 0:
            raise EvidenceGateError("INVALID_BLOCK_DENOMINATOR")
        if not all(isinstance(x, (int, float)) and math.isfinite(x) and x >= 0 for x in (a, b)):
            raise EvidenceGateError("INVALID_PREDICTION_RISK")
        if a <= 0:
            raise EvidenceGateError("ZERO_BASELINE_RISK")
        rows.append((name, n, float(a), float(b)))
    total_n = sum(x[1] for x in rows)
    if total_n != expected_n:
        raise EvidenceGateError("SOURCE_RECORD_DENOMINATOR_DRIFT")
    gains = [a - b for _, _, a, b in rows]
    total_abs = math.fsum(abs(x) for x in gains)
    base = math.fsum(a for _, _, a, _ in rows) / len(rows)
    improved = math.fsum(b for _, _, _, b in rows) / len(rows)
    mean_gain = math.fsum(gains) / len(rows)
    one_block_out = [(math.fsum(gains) - g) / (len(rows)-1) for g in gains]
    return {
        "block_count": len(rows),
        "original_heldout_records": total_n,
        "blocks_with_positive_gain": sum(x > 0 for x in gains),
        "equal_block_baseline_mse": base,
        "equal_block_candidate_mse": improved,
        "equal_block_mse_gain": mean_gain,
        "record_weighted_mse_gain": math.fsum(n*(a-b) for _, n, a, b in rows) / total_n,
        "relative_gain_ratio_of_equal_block_means": mean_gain/base,
        "mean_blockwise_relative_gain": math.fsum((a-b)/a for _, _, a, b in rows)/len(rows),
        "min_omit_one_block_gain": min(one_block_out),
        "max_omit_one_block_gain": max(one_block_out),
        "max_absolute_block_gain_share": max(abs(x) for x in gains)/total_abs if total_abs else 0.0,
        "signed_block_gain": {name: a-b for name, _, a, b in rows},
    }


def build_ladder(sources):
    if set(sources) != set(FILES):
        raise EvidenceGateError("SOURCE_RECEIPT_SET_DRIFT")
    u, up, s, i = (sources[k] for k in ("ulex", "ulex_placebo", "sinapis", "ipomopsis"))
    if (u["audit"]["matched_genomic_flower_fruit_individuals"] != 202
            or u["audit"]["pollinator_census_individual_link"] is not False
            or up["matched_genotypes"] != 202
            or up["decision"] != "NO_ROBUST_INDIVIDUAL_GENOMIC_IDENTITY_INFORMATION_BEYOND_SITE_STRUCTURE"
            or s["audit"]["individuals_with_all_model_fields"] != 52
            or s["audit"]["plant_genetic_alleles_or_ancestors_measured"] is not False
            or i["audit"]["complete_paired_stage_records"] != 367):
        raise EvidenceGateError("ORIGINAL_UNIT_OR_OBSERVABILITY_DRIFT")
    if len(u["heldout_sites"]) != 6 or len(s["heldout_results"]["folds"]) != 5:
        raise EvidenceGateError("HELDOUT_BLOCK_COUNT_DRIFT")
    sf = [
        {"shelter": f["shelter"], "n": f["heldout_plants"], **f["mse_models"]}
        for f in s["heldout_results"]["folds"]
    ]
    comparisons = {
        "ulex_flower_traits_over_elevation": summarize_folds(
            u["heldout_sites"], key="site", count_key="heldout_individuals",
            baseline_key="baseline_elevation_mse",
            model_key="plus_floral_area_weight_quadratic_mse", expected_n=202),
        "ulex_genomic_similarity_over_flower_traits": summarize_folds(
            u["heldout_sites"], key="site", count_key="heldout_individuals",
            baseline_key="plus_floral_area_weight_quadratic_mse",
            model_key="plus_relatedness_residual_mse", expected_n=202),
        "sinapis_flower_shape_over_flower_supply": summarize_folds(
            sf, key="shelter", count_key="n", baseline_key="M0",
            model_key="M1", expected_n=52),
        "sinapis_visits_over_supply_and_shape": summarize_folds(
            sf, key="shelter", count_key="n", baseline_key="M1",
            model_key="M2", expected_n=52),
        "sinapis_visits_over_flower_supply": summarize_folds(
            sf, key="shelter", count_key="n", baseline_key="M0",
            model_key="M3", expected_n=52),
    }
    _assert_close(
        comparisons["ulex_genomic_similarity_over_flower_traits"]["equal_block_mse_gain"],
        u["genetic_increment_conditional_on_traits"]["equal_site_mean_mse_gain"],
        "Ulex v15")
    _assert_close(
        comparisons["ulex_genomic_similarity_over_flower_traits"]["equal_block_mse_gain"],
        up["real_genomic_gain"]["equal_site"], "Ulex v16")
    _assert_close(
        comparisons["ulex_flower_traits_over_elevation"]["equal_block_mse_gain"],
        u["trait_increment"]["equal_site_mean_mse_gain"], "Ulex floral")
    for k, receipt_key in (
        ("sinapis_flower_shape_over_flower_supply", "trait_gain_m0_to_m1"),
        ("sinapis_visits_over_supply_and_shape", "process_gain_m1_to_m2"),
        ("sinapis_visits_over_flower_supply", "process_gain_m0_to_m3")
    ):
        _assert_close(comparisons[k]["equal_block_mse_gain"],
                      s["heldout_results"]["gains"][receipt_key]["equal_shelter"], k)
    placebo = up["within_site_genotype_identity_placebo"]
    if not (placebo["null_replicates_ge_observed_equal_site_gain"] == 312
            and up["null_design"]["replicates"] == 1000
            and placebo["equal_site_genomic_gain"]["q975"]
                > up["real_genomic_gain"]["equal_site"]):
        raise EvidenceGateError("ULEX_GENOMIC_PLACEBO_DRIFT")
    p = s["within_stratum_visitation_placebo"]
    if not (p["replicates"] == 200 and
            p["observed_equal_shelter_increment"] > p["placebo_distribution"]["q975"] and
            p["fraction_placebo_greater_or_equal_observed_descriptive"] == 0):
        raise EvidenceGateError("SINAPIS_VISIT_PLACEBO_DRIFT")
    physical = i["within_context_forecasts"]["physical_location"]["seed_initiation_log1p"]
    if not (i["v13_previous_same_source_record"]["seed_initiation_log1p"]["v13_row_weighted_gain"] > 0
            and physical["heldout_row_weighted_mse_gain"] < 0
            and physical["n_years_width_improves"] == 6
            and physical["readiness_gate"] == "NO_ROBUST_WITHIN_CONTEXT_YEAR_TRANSFER"):
        raise EvidenceGateError("IPOMOPSIS_CONTEXT_CHALLENGE_DRIFT")
    return {
        "schema": "SDEA_v18_cross_study_source_block_information_audit",
        "status": "RETROSPECTIVE_AGGREGATE_RECEIPT_AUDIT_NOT_META_ANALYSIS",
        "inputs": FILES,
        "independent_original_studies": 3,
        "comparison_scales": {
            "Ulex": "squared fruit fraction, six original field sites",
            "Sinapis": "squared log1p(estimated seeds per plant), five rain shelters",
            "Ipomopsis": "squared log1p(seed initiation), nine years (same source as v13)",
        },
        "comparisons": comparisons,
        "adversarial_controls": {
            "Ulex_genomic_identity": {
                "actual_gain": up["real_genomic_gain"]["equal_site"],
                "site_preserving_label_975": placebo["equal_site_genomic_gain"]["q975"],
                "placebo_relabels_ge_real": 312,
                "placebo_relabels_total": 1000,
                "individual_genome_identity_gate_passed": False,
            },
            "Sinapis_visitor_identity": {
                "actual_gain": p["observed_equal_shelter_increment"],
                "within_shelter_treatment_placebo_975": p["placebo_distribution"]["q975"],
                "placebo_relabels_ge_real": 0,
                "placebo_relabels_total": 200,
                "descriptive_within_source_gate_passed": True,
            },
            "Ipomopsis_context": {
                "original_v13_record_weighted_gain": i["v13_previous_same_source_record"]
                    ["seed_initiation_log1p"]["v13_row_weighted_gain"],
                "v14_within_context_record_weighted_gain": physical["heldout_row_weighted_mse_gain"],
                "years_with_positive_within_context_gain": physical["n_years_width_improves"],
                "years_total": 9,
                "strict_within_context_gate_passed": False,
            },
        },
        "joint_G_interaction_phenotype_fitness_same_individual": "NOT_IDENTIFIABLE",
        "trait_regain_across_generations": "NOT_MEASURED",
        "cross_study_numeric_effect_pooling": "PROHIBITED_DIFFERENT_OUTCOMES_CONTRASTS_AND_UNITS",
        "inference_limit": (
            "Delete-one-fold calculations summarize influence of existing overlapping-training "
            "cross-validation folds, not independent re-fitting, p-values, causal effects, "
            "meta-analysis, or evolutionary predictions."
        ),
    }


def main():
    print(json.dumps(build_ladder(load_locked()), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
