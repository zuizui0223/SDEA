# SDEA v13 — Independent botanical original-study test of stage-specific predictive portability

**Date 2026-10-29. Retrospective, exploratory, source-audited.** This branch retains every v1–v12 numerical readout and scientific interpretation without modification. The v9–v12 original observations are two already-published *Cirsium* studies (one USDA *C. undulatum* study, one *C. arvense* author study); v13 opens a **different published original botanical study** in *Ipomopsis*. This is **not** an independently observed replication of the v12 bagging intervention, nor proof of latent ancestral floral regain.

## 1. Research question and exact negative cross-study gate

Does an observed flower morphology coordinate contain *predictive* information beyond source phylo-ecological background when predicting two distinct reproductive-stage outcomes across previously unseen years?

- **Original source:** Campbell et al. (2022), *American Naturalist*, DOI `10.1086/716740`, original registered data `10.7280/D1KM49`, Zenodo record `5083399`, original author CSV `Selection_by_pollinators_and_predators.csv`, published MD5 `12ed3dbe6d269c5f808f9c6b230d468c`.
- **Earlier source audit:** aza3 `data/evidence/ipomopsis_two_stage_trait_selection_independent_original_js_v1.json` at commit `d54aba2dcac3ecbfd4c4fb6bd5b0b88562f3595d`, Git blob `df57a0f13ac30ab2253969917eb5617b1cf64552`. The receipt's 437 rows, 12 columns, 367 paired width/stage records, site/year/type counts, source missingness, and 41 `seeds > seedsinit` records are source identity checks. That preexisting audit **does not** count as a second independent source study; it already reported exploratory within Site×Year×Type associations and leave-one-year-out correlations, not the new held-out-year *proper predictive-score* model.
- **v12 non-equivalence:** *C. arvense* 2022 has a compound organza bag vs Control intervention and viable seeds / original flowers, and insects emerging / head. The *Ipomopsis* file has *observed* corolla width and `seedsinit` plus `flyescape`. There is **no common randomized treatment**, no common reproductive denominator and no head/plant identity crosswalk. Consequently no direct cross-study bagging effect prediction or pooled selection coefficient is identified. **Do not make one up.**

## 2. Sources, eligibility and inference unit

Public source raw bytes are downloaded **only inside a reproducible GitHub Actions job**, using immutable/original author registration and published MD5 plus an independently SHA-pinned aza3 summary. If source fetch fails, the CSV's MD5 changes, the aza3 Git blob changes, or the expected original schema/denominators drift, fail closed: **`ORIGINAL_SOURCE_UNAVAILABLE_OR_DRIFT`**, never zero ecological effect.

The original published file has **437 rows**, 12 columns `Site,Year,Type,width,color,nectar,sepal,totpinene,indoleN,seeds,seedsinit,flyescape`. Source years 2010–2018, original `Site` codes `agg`, `hyb`, `ten`, and `Type` A/H/T. Per the original author's dataset meaning, `hyb` and `ten` are labels from the same upper experimental field: harmonize them to **one** physical-location indicator `upper`, versus `agg`; still retain the **original** site codes in the locked count audit. Do not represent the three labels as three spatially independent sites.

All nine whole-year folds use the **same complete-case rows for both stages**, requiring nonmissing measured `width`, `seedsinit`, and `flyescape`: expected 367 original rows. No missing response is imputed as zero; each row's original `Type`, `Site`, `Year` is retained. The source does **not** expose stable `plantID`, and repeated measures, genetic individuality, pedigree and between-year identity cannot be verified. `437` is a *source record count*, not 437 independent plants or evolutionary origins. The source also has 41 rows where `seeds > seedsinit`; do not claim the 2 stage variables satisfy an exact per-seed population conservation equation or infer individual eaten seeds.

## 3. Endpoint-specific nested predictive models

The **entire original source year** is withheld for testing, nine folds (2010–2018). For each fold train on the other eight source years; never inspect test outcomes to choose model, parameters or transformation. Two response functions are evaluated independently on the same 367 eligible original rows:

1. **Formation stage:** `Y_1 = log1p(seedsinit)` (initiated seeds, not an individual pollination-mediated causal effect). The shared baseline design is intercept + physical location indicator + source `Type` H/T indicators (A reference). An augmented model adds **training-standardized** measured corolla width.
2. **Fruit-fly escape stage:** `Y_2 = flyescape` (author's fraction of **fruits** escaping fruit-fly predation, not seeds and not viable achene survival). Its baseline includes the same physical location and Type plus **training-standardized** `log1p(seedsinit)` as a preceding measured stage; the augmented model adds measured **training-standardized** corolla width. Conditional prediction here is **not** evidence of causal direct effects or formal mediation by seed initiation. Clamp predicted fractions to `[0,1]` for both arms only.

Use identical closed-form **ridge least squares penalty = 1.0**, on all non-intercept coefficients, for baseline and augmented models in both stages. The penalty and feature engineering are frozen in source. Fit centering/scaling on training rows **only**. No year indicator is fitted (test year is unseen). No hyperparameter selection by the held-out year is allowed. The original-source file bytes and all nine fold labels determine every number.

For each of nine years and both endpoints report: training N, held-out N, baseline squared-error MSE, augmented MSE, and `gain = MSE_baseline − MSE_augmented` (positive improves). Aggregate (a) equal-year mean MSE gain and (b) held-out-row-weighted mean MSE gain. A **strong source-internal readiness gate**, not a significance test, requires both gains > 0 **and** positive year-specific gain in at least seven of nine years *separately at a given endpoint*. This 7/9 criterion is an exploratory stress gate, not a prospective preregistration or an empirical universal transfer threshold.

Do **not** sum/average MSEs for `log1p(seedsinit)` and `flyescape`: scales, meaning and biological denominators differ. Any mixed, zero or negative result is retained. Do not recode absence of measured genomic state or visitor identity as failure to evolve.

## 4. Independence, source provenance and negative controls

- The 2010–2018 samples originate in *one* original published common-garden/source study. Nine year folds are not nine independent species, field experiments or evolutionary origins. Even though *Ipomopsis* is separate from the two *Cirsium* sources, neither source here has the same intervention and response semantics necessary for transport of v12's bagging coefficient.
- A source count / MD5 / Git blob mismatch is a **source gate error**, not a biological null. Synthetic data tests pass through the same parser only with explicit `strict=False` and cannot be substituted for the original source result.
- Both baseline and augmented training regressions always use **identical rows**, fixed background columns and penalty. Altering held-out features must not alter training estimates; no outcome-based row pruning. Complete-case selection is an observational availability constraint, and its bias cannot be removed with these records alone.
- The `flyescape` numerator is fruits, not flowerheads or seeds. Avoid multiplying `seedsinit × flyescape` as a true seed accounting identity. Do not label increased `flyescape` observed plant fitness if seed formation is absent or if observational follow-up is incomplete.
- Generate no per-observation output or original raw-file copy in SDEA: version only the source manifest, executable standard-library code, unit tests, and an **aggregate** file with nine year folds and endpoint-specific metrics.
- If the corresponding existing source receipt reports signs/correlations, those are **prior same-data results** (audit benchmarks only). The v13 question concerns held-out-year predictive MSE conditional on initial-stage observations, a different estimand. Do not claim this as an independent original-data replication of aza3's analyses.
- All genomic ancestral toolkit, functional-pigment pathway, molecular origin and multigeneration offspring establishment layers remain **unmeasured**. The high-level SDEA genetic re-evolution theory remains conditional/synthetic, regardless of the result of this stage-specific predictive test.

## 5. Reproduce

Python 3.11+, standard-library only. In CI:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_external_stage_transfer.py --fetch-pinned --output /tmp/sdea_v13.json
python scripts/verify_fixture.py results/v13_ipomopsis_stage_transfer.json /tmp/sdea_v13.json
```

The final fixture comparison is added **after** an initial successful source-authenticated run, keeping the first source-observed result visible and preventing an outcome-sensitive silent revision. In offline environments, use `--original` and `--aza3-receipt` to point to **separately supplied identical original pinned bytes**; output remains source-checked. If external source is unavailable, do not fabricate an original outcome or publish a frozen zero.

**Claim ceiling:** a new independent botanical published-original-data held-out-year predictive *measurement test*; not a causal cross-species effect comparison, not a new observed natural adaptive evolutionary recovery, and not an independently registered outside-sample forecast.
