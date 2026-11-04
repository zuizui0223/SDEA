"""SDEA v12: original 2022 Cirsium multi-site, held-out-site transfer.

This is a *retrospective* original-source predictive audit, not a replicated
organism, a controlled pollen/seed mediation effect, or genetic re-evolution.
The source is always pinned by exact Git commit and Git blob of raw bytes.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
import math
from urllib.request import Request, urlopen

SITES = ("campus", "construction", "ditch", "rapseed")
ARMS = ("Control", "Pollinator Exclusion", "Predator Exclusion", "Predator Model")
MEASURES = ("viable_fraction", "emerged_insects_per_head")
SOURCE_NUMERIC = ("Flowers", "Viable_Seeds", "NViable_Seeds", "Damaged_Seeds", "Missing_Seeds", "Insects")
EXPECTED_BLOBS = {
    "original_author_seeds": "ca99c762162d15f006f0c73b4603ca71fd320a29",
    "aza3_independent_receipt": "740bfe83c881ee61269342070253b3baeaf58aa2",
}


class SiteTransferError(ValueError):
    """No source drift, schema corruption or faulty denominator as a zero effect."""


def git_blob_sha(raw: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw).hexdigest()


def pinned_source_url(entry: dict) -> str:
    if (entry.get("repo") not in ("Yogita25/Multitrophic_Interactions", "zuizui0223/aza3")
            or len(entry.get("commit", "")) != 40
            or entry.get("path") not in ("seeds.csv", "data/evidence/cirsium_arvense_four_arm_joint_response_independent_js_v1.json")):
        raise SiteTransferError("unrecognized author/source archive pinned path")
    return "https://raw.githubusercontent.com/" + entry["repo"] + "/" + entry["commit"] + "/" + entry["path"]


def verified_input(entry: dict, *, raw: bytes | None = None, fetch_pinned: bool = False) -> bytes:
    if raw is None:
        if not fetch_pinned:
            raise SiteTransferError("real source not provided; explicitly use --fetch-pinned")
        try:
            request = Request(pinned_source_url(entry), headers={"User-Agent": "SDEA-v12-original-source-audit/1.0"})
            with urlopen(request, timeout=60) as response:
                raw = response.read()
        except Exception as exc:
            raise SiteTransferError("SOURCE_UNAVAILABLE_OR_DRIFT: " + str(exc)) from exc
    if not isinstance(raw, bytes) or git_blob_sha(raw) != entry["blob_sha"]:
        raise SiteTransferError("ORIGINAL_SOURCE_GIT_BLOB_MISMATCH")
    return raw


def parse_original_seeds(raw: bytes, *, strict: bool = True) -> tuple[list[dict], dict]:
    try:
        reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
        required = set(SOURCE_NUMERIC) | {"Site_Code", "Block_Code", "Treatment_Code", "Vial_Code"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise SiteTransferError("original seeds source missing required columns")
        rows = list(reader)
    except UnicodeDecodeError as exc:
        raise SiteTransferError("original CSV not UTF-8") from exc
    if strict and len(rows) != 615:
        raise SiteTransferError("original source must have exactly 615 head rows")
    for row in rows:
        numbers = {}
        for col in SOURCE_NUMERIC:
            try:
                n = float(row[col])
            except (ValueError, TypeError, KeyError) as exc:
                raise SiteTransferError("missing/non-numeric " + col) from exc
            if not math.isfinite(n) or not n.is_integer() or (col != "Missing_Seeds" and n < 0):
                raise SiteTransferError("invalid seed head count in " + col)
            numbers[col] = int(n)
        if numbers["Flowers"] <= 0:
            raise SiteTransferError("zero floral exposure")
        if numbers["Flowers"] != sum(numbers[z] for z in ("Viable_Seeds", "NViable_Seeds", "Damaged_Seeds", "Missing_Seeds")):
            raise SiteTransferError("original accounting identity not satisfied")
        row["_numeric"] = numbers
    eligible = [r for r in rows if r["Site_Code"] in SITES and r["Treatment_Code"] in ARMS]
    if strict and len(eligible) != 612:
        raise SiteTransferError("eligible original head row count drift")
    return eligible, {"raw_head_rows": len(rows), "eligible_head_rows": len(eligible),
                      "ineligible_head_rows": len(rows)-len(eligible),
                      "negative_signed_missing_raw_heads": sum(r["_numeric"]["Missing_Seeds"] < 0 for r in rows)}


def aggregate_complete_blocks(eligible: list[dict], *, strict: bool = True) -> tuple[list[dict], dict]:
    grouped = defaultdict(lambda: {"head_count":0, "flowers":0, "viable":0, "insects":0})
    for row in eligible:
        key = (row["Site_Code"], row["Block_Code"], row["Treatment_Code"])
        if not key[1] or key[0] not in SITES or key[2] not in ARMS:
            raise SiteTransferError("invalid source site/block/arm in included rows")
        v = grouped[key]
        v["head_count"] += 1
        v["flowers"] += row["_numeric"]["Flowers"]
        v["viable"] += row["_numeric"]["Viable_Seeds"]
        v["insects"] += row["_numeric"]["Insects"]
    blocks = defaultdict(dict)
    for (site, block, arm), v in grouped.items():
        if v["flowers"] <= 0 or v["head_count"] <= 0:
            raise SiteTransferError("empty inferred plant-proxy denominator")
        blocks[(site,block)][arm] = {
            "viable_fraction": v["viable"] / v["flowers"],
            "emerged_insects_per_head": v["insects"] / v["head_count"],
        }
    complete = [dict(site=site, block=block, arms=arms)
                for (site,block),arms in sorted(blocks.items())
                if all(arm in arms for arm in ARMS)]
    site_n = {site:sum(row["site"] == site for row in complete) for site in SITES}
    if strict:
        if (len(grouped),len(blocks),len(complete),site_n) != (
                164,45,36,{"campus":12,"construction":9,"ditch":8,"rapseed":7}):
            raise SiteTransferError("original inferred-proxy/block eligibility drift")
        arm_head_n = Counter(row["Treatment_Code"] for row in eligible)
        if arm_head_n != Counter({"Control":164,"Pollinator Exclusion":142,"Predator Exclusion":154,"Predator Model":152}):
            raise SiteTransferError("original head treatment rows drift")
    return complete, {"inferred_site_block_arm_plant_proxies":len(grouped),
                      "source_site_blocks":len(blocks),"four_arm_complete_blocks":len(complete),
                      "complete_site_blocks":site_n}


def _avg(x: list[float]) -> float:
    if not x:
        raise SiteTransferError("cannot predict with an empty training/test set")
    return sum(x)/len(x)


def paired_contrast(block: dict, *, intervention: str, metric: str) -> float:
    if intervention not in ARMS[1:] or metric not in MEASURES:
        raise SiteTransferError("unsupported comparator or outcome")
    return block["arms"][intervention][metric] - block["arms"]["Control"][metric]


def _direction(value: float) -> str:
    if value < -1e-12: return "negative"
    if value > 1e-12: return "positive"
    return "equal"


def heldout_site_evaluation(blocks: list[dict], *, intervention: str, metric: str) -> dict:
    """Fold-specific prediction uses three OTHER sites only; no target leakage.

    Site-balanced primary scores use four held-out SITE MEANS as the four
    units, not 36 head/plant proxies treated as four independent studies.
    Per-block scores remain a descriptive secondary outcome.
    """
    if len({b["site"] for b in blocks}) != len(SITES):
        raise SiteTransferError("the 4-site external holdout requires all source sites")
    rows = []
    for site in SITES:
        train = [paired_contrast(b,intervention=intervention,metric=metric)
                 for b in blocks if b["site"] != site]
        test = [paired_contrast(b,intervention=intervention,metric=metric)
                for b in blocks if b["site"] == site]
        pred = _avg(train)
        actual = _avg(test)
        rows.append({"heldout_site":site, "heldout_blocks":len(test),"training_blocks":len(train),
                     "observed_site_mean":actual,"predicted_from_other_sites":pred,
                     "sign_matched":_direction(actual)==_direction(pred),
                     "zero_site_squared_error":actual*actual,
                     "heldout_site_squared_error":(actual-pred)**2,
                     "zero_block_mse":_avg([x*x for x in test]),
                     "heldout_block_mse":_avg([(x-pred)**2 for x in test])})
    site_null=_avg([r["zero_site_squared_error"] for r in rows])
    site_pred=_avg([r["heldout_site_squared_error"] for r in rows])
    if site_null == 0:
        gain = None
    else:
        gain = 1-site_pred/site_null
    return {
        "intervention_minus_control":intervention,"measure":metric,
        "holdout_unit":"one whole original study site",
        "site_mean_mse_zero_reference":site_null,
        "site_mean_mse_train_other_sites":site_pred,
        "site_mean_fractional_squared_error_reduction_vs_zero":gain,
        "site_benefit_count":sum(r["heldout_site_squared_error"] < r["zero_site_squared_error"] for r in rows),
        "site_sign_match_count":sum(r["sign_matched"] for r in rows),
        "block_mse_zero_secondary":sum(r["zero_block_mse"]*r["heldout_blocks"] for r in rows)/sum(r["heldout_blocks"] for r in rows),
        "block_mse_train_other_sites_secondary":sum(r["heldout_block_mse"]*r["heldout_blocks"] for r in rows)/sum(r["heldout_blocks"] for r in rows),
        "folds":rows,
    }


def joint_sign_counts(blocks: list[dict], intervention: str) -> dict[str,int]:
    pairs=Counter()
    for b in blocks:
        v=_direction(paired_contrast(b,intervention=intervention,metric="viable_fraction"))
        i=_direction(paired_contrast(b,intervention=intervention,metric="emerged_insects_per_head"))
        label={"negative":"down", "positive":"up", "equal":"equal"}
        pairs["viable_"+label[v]+"__emerged_"+label[i]] += 1
    if sum(pairs.values())!=len(blocks):
        raise SiteTransferError("joint sign block partition failed")
    return dict(sorted(pairs.items()))


def validate_prior_independent_receipt(report: dict, results: dict) -> dict:
    """Cross-check with aza3's old source report, not a new biological cohort."""
    if report.get("version") != "cirsium_arvense_four_arm_joint_response_independent_js_v1":
        raise SiteTransferError("aza3 published independent receipt has unknown version")
    if (report.get("source_counts",{}).get("heads_all") != 615 or
        report["source_counts"]["complete_four_arm_blocks"] != 36 or
        report.get("source",{}).get("seeds_git_blob_sha") != EXPECTED_BLOBS["original_author_seeds"]):
        raise SiteTransferError("aza3 original source counts or blob disagree")
    labels={"Pollinator Exclusion":"bag_minus_control", "Predator Exclusion":"pred_exclusion_minus_control"}
    for intervention, key in labels.items():
        upstream=report["comparisons"][key]
        for metric in MEASURES:
            own=results[intervention][metric]
            expected=upstream["measures"][metric]["site_means"]
            for fold in own["folds"]:
                if not math.isclose(fold["observed_site_mean"],expected[fold["heldout_site"]],abs_tol=1e-11):
                    raise SiteTransferError("new raw-data output disagrees with previous receipt for "+key+"/"+metric)
        if upstream["joint_response_counts"] != results[intervention]["joint_signs"]:
            raise SiteTransferError("new raw-data joint sign count disagrees with prior receipt")
    return {"prior_receipt_is_not_new_independent_data":True,"source_crosscheck_passed":True,
            "prior_receipt_version":report["version"],
            "prior_original_source_blob":report["source"]["seeds_git_blob_sha"]}


def run_original(raw: bytes, prior: bytes, manifest: dict) -> dict:
    src=manifest["sources"]
    if set(src) != set(EXPECTED_BLOBS) or any(src[k]["blob_sha"] != EXPECTED_BLOBS[k] for k in src):
        raise SiteTransferError("source pin manifest does not match approved originals")
    verified_input(src["original_author_seeds"],raw=raw)
    verified_input(src["aza3_independent_receipt"],raw=prior)
    eligible,basic=parse_original_seeds(raw)
    blocks,unit=aggregate_complete_blocks(eligible)
    results={}
    for arm in ("Pollinator Exclusion","Predator Exclusion"):
        results[arm]={metric:heldout_site_evaluation(blocks,intervention=arm,metric=metric)
                      for metric in MEASURES}
        results[arm]["joint_signs"]=joint_sign_counts(blocks,arm)
    try:
        receipt=json.loads(prior.decode("utf-8-sig"))
    except (ValueError,UnicodeError) as exc:
        raise SiteTransferError("aza3 prior receipt is not JSON") from exc
    integrity=validate_prior_independent_receipt(receipt,results)
    bag=results["Pollinator Exclusion"]["viable_fraction"]
    insect=results["Pollinator Exclusion"]["emerged_insects_per_head"]
    neg=results["Predator Exclusion"]["viable_fraction"]
    gate="CONDITIONAL_ACCESS_COMPOSITE_DIRECTION_PORTABLE_ONLY" if (
        bag["site_benefit_count"] == 4 and bag["site_sign_match_count"] == 4 and
        all(f["observed_site_mean"] < 0 for f in bag["folds"]) and
        insect["site_sign_match_count"] == 4 and
        all(f["observed_site_mean"] < 0 for f in insect["folds"]) and
        neg["site_benefit_count"] < 4
    ) else "NOT_CONFIRMED_SOURCE_STUDY_PATTERN"
    return {
        "protocol":"SDEA v12 original 2022 four-arm Cirsium arvense four-site-out audit (2026-10-29)",
        "evidence_tier":"REAL_ORIGINAL_ARCHIVED_CROSS_SITE_PREDICTION; ONE STUDY; NOT INDEPENDENT_FIELD_REPLICATION",
        "source":{k:{"repo":v["repo"],"commit":v["commit"],"path":v["path"],"blob_sha":v["blob_sha"],"verified":True} for k,v in src.items()},
        "audit_counts":{**basic,**unit,"distinct_study_sites":4},
        "comparisons":results,"crosscheck":integrity,
        "decision":gate,
        "strict_limits":[
            "Bagging is a compound mesh-exclusion intervention, not randomized pollen contact alone",
            "Insect emergence counts are not individual seed predators or tested pollinator visits",
            "Plant proxies are site-block-arm, not proven exact maternal plant identities",
            "36 matched blocks nested in FOUR sites, not 36 independent populations or studies",
            "Held-out sites belong to the same preprint/year; no new taxon and no prospective validation",
            "No observed phylogeny, molecular floral toolkit, re-expression history, or descendant establishment",
            "The zero-reference comparator is a forecasting baseline, not a biological null of zero treatment effect",
            "Joint sign counts are descriptive: do not infer cross-partner causation or evolutionary route redundancy"
        ],
    }
