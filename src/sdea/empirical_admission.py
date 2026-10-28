"""SDEA v9: genuinely source-pinned *partial* natural-data evidence gate.

CC0 USDA mother-plant records are read from the original aza3-controlled mirror
at an immutable Git commit. Original source data are never copied into SDEA.
EGWE and ANTO locked source summaries are cross-checked but NOT re-analysed as
if they were independently replicated. Izu-core remains a synthetic programme.
No complete ancestral floral reexpression plus fitness sequence is inferred.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from pathlib import Path
from urllib.request import Request, urlopen

SOURCE_FIELDS = ("plt", "year", "trt", "plt.gsd", "ap.hd", "o.hds")
ALL_TREATMENTS = ("di", "dni", "ndi", "ndni")
CORRUPT_ORIGINAL_ROW = ("2007", "327")


class SourceGateError(ValueError):
    """Fail closed on input/schema/source corruption; never call it ecological absence."""


def blob_sha(payload: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(payload)).encode() + b"\0" + payload).hexdigest()


def pin_url(entry: dict) -> str:
    return ("https://raw.githubusercontent.com/zuizui0223/" + entry["repo"] + "/"
            + entry["commit"] + "/" + entry["path"])


def validate_source_bytes(entry: dict, raw: bytes) -> bytes:
    if blob_sha(raw) != entry["blob_sha"]:
        raise SourceGateError("source Git blob SHA mismatch for " + entry["path"])
    return raw


def load_inputs(lock: dict, *, root: Path | None = None, allow_network: bool = False) -> dict[str, bytes]:
    if root is None and not allow_network:
        raise SourceGateError("provide an original local source root or --fetch-pinned")
    result = {}
    for key, entry in lock["sources"].items():
        if root is not None:
            file = root / entry["repo"] / entry["commit"] / entry["path"]
            if not file.is_file():
                raise SourceGateError("pinned local source is missing: " + str(file))
            raw = file.read_bytes()
        else:
            try:
                req = Request(pin_url(entry), headers={"User-Agent": "SDEA-v9-pinned-source-validation/1.0"})
                with urlopen(req, timeout=45) as response:
                    raw = response.read()
            except Exception as exc:
                raise SourceGateError("SOURCE_UNAVAILABLE: " + key + ": " + str(exc)) from exc
        result[key] = validate_source_bytes(entry, raw)
    return result


def parse_maternal_joint(raw: bytes, *, check_original: bool = True) -> dict:
    try:
        table = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
        if table.fieldnames is None or any(k not in table.fieldnames for k in SOURCE_FIELDS):
            raise SourceGateError("mother-plant CSV lacks a required source column")
        rows = list(table)
    except UnicodeDecodeError as exc:
        raise SourceGateError("mother-plant source not UTF-8") from exc
    if check_original and len(rows) != 106:
        raise SourceGateError("expected 106 original USDA maternal records")
    seen = set()
    omitted = []
    valid = []
    for row in rows:
        year, plant, trt = row["year"].strip(), row["plt"].strip(), row["trt"].strip()
        if year not in ("2007", "2008") or trt not in ALL_TREATMENTS or not plant:
            raise SourceGateError("unexpected original year/treatment/maternal ID")
        key = (year, plant)
        if key in seen:
            raise SourceGateError("duplicate original mother key: " + str(key))
        seen.add(key)
        try:
            vals = []
            for col in ("plt.gsd", "ap.hd", "o.hds"):
                x = float(row[col])
                if not math.isfinite(x) or x < 0 or not x.is_integer():
                    raise SourceGateError("invalid original integer good-seed count")
                vals.append(int(x))
        except (TypeError, ValueError) as exc:
            raise SourceGateError("invalid original maternal count") from exc
        total, apical, lateral = vals
        if apical + lateral != total:
            omitted.append({"year": year, "plant_id": plant, "apical": apical,
                            "lateral": lateral, "total": total,
                            "reason": "apical_plus_lateral_inconsistent_with_maternal_total"})
            continue
        valid.append((year, trt, apical > 0, lateral > 0))
    if check_original:
        if len(valid) != 105 or len(omitted) != 1 or (omitted[0]["year"], omitted[0]["plant_id"]) != CORRUPT_ORIGINAL_ROW:
            raise SourceGateError("unexpected coherence exclusions in original source")
    grouped: dict[str, dict] = {}
    for year, treatment, a, l in valid:
        key = year + "|" + treatment
        group = grouped.setdefault(key, {"n": 0, "neither": 0, "apical_only": 0,
                                         "lateral_only": 0, "both": 0})
        group["n"] += 1
        group["both" if a and l else "apical_only" if a else "lateral_only" if l else "neither"] += 1
    if check_original and set(grouped) != {y + "|" + t for y in ("2007", "2008") for t in ALL_TREATMENTS}:
        raise SourceGateError("eight expected original year-treatment strata not represented")
    out = {}
    for group, c in sorted(grouped.items()):
        n = c["n"]
        ap = (c["apical_only"] + c["both"]) / n
        lat = (c["lateral_only"] + c["both"]) / n
        observed = 1 - c["neither"] / n
        expected_independent = ap + lat - ap * lat
        overlap = c["both"] / n
        out[group] = {**c, "apical_success": ap, "lateral_success": lat,
                      "either_success": observed,
                      "independent_marginal_either": expected_independent,
                      "observed_minus_independent": observed - expected_independent,
                      "same_mother_covariance": overlap - ap * lat,
                      "observed_two_route_added_coverage": observed - max(ap, lat)}
    return {"source_mothers": len(rows), "admitted_coherent_mothers": len(valid),
            "excluded_incoherent": omitted, "year_treatment_strata": out,
            "unit": "one original year-specific maternal plant",
            "claim_ceiling": "two simultaneously observed anatomical good-seed routes; not latent genotype/reexpression or causal removal insurance"}


def summarize_source_receipts(inputs: dict[str, bytes]) -> dict:
    def j(key: str) -> dict:
        try:
            return json.loads(inputs[key].decode("utf-8-sig"))
        except (UnicodeDecodeError, ValueError, KeyError) as exc:
            raise SourceGateError("invalid locked JSON receipt: " + key) from exc
    ryuk = j("ryukyu_readiness")
    izu = j("izu_summary")
    oe = j("oenothera_summary")
    reg = j("anto_registry")
    model_doc = inputs["izu_core_boundary"].decode("utf-8-sig")

    ryuk_info = ryuk["source_summary"]
    if (ryuk_info["public_run_count"] != 25 or
        ryuk_info["floral_RNA_runs_in_this_source"] != 0 or
        ryuk_info["phaseA_collected_here"] != 0 or
        len(ryuk["per_taxon"]) != 12):
        raise SourceGateError("unexpected Ryukyu readiness source schema")
    if any(t["same_plant_verified_functional_toolkit"] != "NOT_MEASURED_IN_REGISTERED_SOURCE" or
           t["heritable_repeated_reexpression"] != "NOT_MEASURED_IN_REGISTERED_SOURCE"
           for t in ryuk["per_taxon"]):
        raise SourceGateError("unexpected same-plant molecular or descendant claim in Ryukyu registry")
    izu_c1 = izu["aggregate_metrics"]["C1"]["row_weighted_mse"]
    izu_c2 = izu["aggregate_metrics"]["C2"]["row_weighted_mse"]
    if izu["analysis_lock"]["row_count"] != 572 or izu["analysis_lock"]["site_count"] != 8:
        raise SourceGateError("EGWE Honshu-Izu source cohort has changed")
    if oe["schema"]["n_seed_family_rows"] != 60 or oe["schema"]["n_maternal_plants"] != 23:
        raise SourceGateError("EGWE Oenothera source cohort has changed")
    mating_m0 = oe["models"]["M0"]["loo_mse"]
    mating_m1 = oe["models"]["M1"]["loo_mse"]
    cases = reg["cases"]
    if len(cases) != 6 or reg["dataset_identity"] != "non_systematic_source-verified_mechanism_examples_NOT_random_species_sample":
        raise SourceGateError("ANTO case source differs from nonrandom fixed registry")
    petunia_supported = [c for c in cases if c["clade"] == "Petunia" and c["case_role"].startswith("MOLECULARLY_SUPPORTED")]
    if len(petunia_supported) != 2 or any(c["trait_specific_viable_seed_fitness"] != "NOT_QUALIFIED" for c in petunia_supported):
        raise SourceGateError("unexpected Petunia literature mechanism/fecundity qualification")
    if "finite diploid simulation study" not in model_doc or "not a direct historical estimate" not in model_doc:
        raise SourceGateError("izu-core model versus natural-result boundary not found")
    return {
        "ryukyu_leaf_RNA": {"public_young_leaf_runs": ryuk_info["young_leaf_runs"],
                            "public_flower_stage_RNA_runs": ryuk_info["floral_RNA_runs_in_this_source"],
                            "collected_phaseA_original_individuals": ryuk_info["phaseA_collected_here"],
                            "prospective_phaseA_target_n_not_real_samples": ryuk_info["phaseA_target_individuals"],
                            "same_individual_flower_toolkit_and_descendant_regain_test_ready": False},
        "honshu_izu_source_reported": {"holdout_sites": izu["analysis_lock"]["site_count"],
                          "pollen_observation_rows": izu["analysis_lock"]["row_count"],
                          "functional_partial_state_mse": izu_c1,
                          "plus_mainland_distance_mse": izu_c2,
                          "adding_distance_mse_change_percent": 100*(izu_c2/izu_c1-1),
                          "independently_refit_by_SDEA": False,
                          "future_trait_reexpression_endpoint": False},
        "oenothera_source_reported": {"maternal_plants": oe["schema"]["n_maternal_plants"],
                            "seed_family_rows": oe["schema"]["n_seed_family_rows"],
                            "treatment_only_LOO_MSE": mating_m0,
                            "plus_isolation_LOO_MSE": mating_m1,
                            "adding_isolation_mse_change_percent": 100*(mating_m1/mating_m0-1),
                            "independently_refit_by_SDEA": False,
                            "future_trait_reexpression_endpoint": False},
        "anto_literature_source_reported": {"nonrandom_cases": len(cases),
                            "petunia_molecular_mechanism_cases": len(petunia_supported),
                            "independent_genera_of_these_two": len({c["clade"] for c in petunia_supported}),
                            "petunia_cases_with_qualified_trait_fitness": 0,
                            "population_reexpression_rate_estimable": False},
        "izu_core": {"model_only": True, "count_as_natural_replication": False}
    }


def run(lock: dict, original_payloads: dict[str, bytes]) -> dict:
    manifest_keys = set(lock["sources"])
    if set(original_payloads) != manifest_keys:
        raise SourceGateError("source manifest and payload set differ")
    for key, raw in original_payloads.items():
        validate_source_bytes(lock["sources"][key], raw)
    maternals = parse_maternal_joint(original_payloads["mother_csv"])
    summaries = summarize_source_receipts(original_payloads)
    statuses = {
        "Cirsium_undulatum_original_mothers": "DIRECT_ORIGINAL_CONTEMPORARY_ANATOMICAL_SEED_FUNCTION_ONLY",
        "Ryukyu_Cirsium_leaf_RNA": "PUBLIC_LEAF_METADATA_NO_SAME_INDIVIDUAL_FLOWER_FUNCTION",
        "Honshu_Izu_ecological_holdout": "EGWE_SOURCE_REPORTED_CURRENT_ECOLOGICAL_PARTIAL_STATE",
        "Oenothera_mating_holdout": "EGWE_SOURCE_REPORTED_CONTEMPORARY_MATING_STATE",
        "Petunia_mechanism_examples": "ANTO_LITERATURE_SELECTED_MECHANISMS_NO_FITNESS_SEQUENCE"
    }
    # The public names deliberately differ from internal stable record keys.
    if len(statuses) != len(lock["admission_units"]):
        raise SourceGateError("admission unit count changed")
    return {"protocol": lock["protocol"],
            "scope": "one independent raw-data replay plus four source-reported bounded evidence units; never pooled as taxon prevalence",
            "source_pins": {k:{"repo":v["repo"],"commit":v["commit"],
                               "path":v["path"],"git_blob_sha":v["blob_sha"],
                               "verified": True} for k,v in lock["sources"].items()},
            "original_reanalysis": {"usda_maternal_good_seed_routes": maternals},
            "source_reported_partial_results": summaries,
            "natural_evidence_admission": statuses,
            "full_cross_layer_contract_admitted_units": 0,
            "denominator_is_convenience_heterogeneous_evidence_units": len(statuses),
            "full_contract_decision": "NOT_IDENTIFIABLE_FROM_THESE_SOURCES",
            "disallowed_inferences": ["ancestral_flower_reexpression_proven_in_Cirsium",
              "maternal_apical_lateral_functions_are_distinct_genetic_ancestral_modules",
              "EGWE_representations_independently_refit_or_cross_validated_by_SDEA",
              "two_Petunia_cases_are_independent_plant_genera",
              "five_evidence_units_are_randomly_sampled_species",
              "izu_core_synthetic_Model_3_counts_as_natural_replication",
              "repeated_flowering_equals_demographic_genetic_fixation"]}
