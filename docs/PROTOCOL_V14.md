# SDEA v14 — does v13 corolla-width forecasting survive **within-year, within-type, within-field** controls?

**Status:** 2026-10-29 retrospective post-v13 exploratory challenge. This is *the same Campbell et al. 2022 Ipomopsis original public study* as v13, **not another biological replication**. Prior v1–v13 results, including v13's positive 7/9-year seed-initiation stage and negative 5/9-year fruit-fly escape stage, are protected and must not be rewritten or pooled into a causal selection effect.

## Question and the failure mode

The v13 baseline used source Type and physical field but not *within-year × field × Type* variation, and the observed 367 complete records span years with changing type/field composition. A year-held-out overall `width` predictive gain may therefore be explained by a width association between ecological/contextual cells, not individual-level within-cell information.

This version evaluates the specific conditional claim:

> Does measured corolla width add heldout-year prediction when the width slope is learned **only from within-original-Year × source-field × source-Type differences**, with all test-year feature centering performed without reading the heldout response?

This is an observational predictive-adequacy and transportability challenge, not a phenotypic-selection coefficient, a field intervention, reproduction mediated by a verified pollinator, an independent genetic ancestry replication or an ability to restore latent floral traits.

## Source and cohort hard gates

Original Campbell et al. 2022 *American Naturalist*, DOI `10.1086/716740`, DOI `10.7280/D1KM49`, original Zenodo 5083399 file `Selection_by_pollinators_and_predators.csv` with **published MD5 `12ed3dbe6d269c5f808f9c6b230d468c`**. Prior aza3 same-source source-audit receipt at immutable commit `d54aba2dcac3ecbfd4c4fb6bd5b0b88562f3595d`, expected Git blob `df57a0f13ac30ab2253969917eb5617b1cf64552`. All inputs fetched read-only by a network-authorized GitHub Actions job. If the author file is unavailable, a profile/MD5 changes, the prior receipt Git SHA changes, or the previously frozen v13 JSON is not reproduced on the same original bytes, the **v14 job fails closed**, not a synthetic success or biological null.

The same 367 rows with observed width/`seedsinit`/`flyescape` from 437 published source records are used for BOTH stages. Nine heldout whole-year folds 2010–2018. Original `Site` codes `agg`, `hyb`, `ten`, source `Type` A/H/T. `hyb` and `ten` are **not distinct physical sites**: collapse these to the one `upper` physical field for the primary context condition; retain three source site labels only as an *additional stricter annotation* sensitivity. Original `plantID` is unavailable and actual individual independence/within-plant repeated observations cannot be verified.

## One-row unit, and two prospectively testable retrospective sensitivity modes

Mode 1 **physical_location**: grouping key `(Year, physical_field∈{agg,upper}, Type)`.

Mode 2 **original_site_code**: grouping key `(Year, original Site∈{agg,hyb,ten}, Type)`; **not** interpreted as three independent physical sites.

All nine whole years are held out once. For each fold, the **baseline** is exactly the old v13 train-only Type + two-physical-field model: stage 1 predicts `log1p(seedsinit)`; stage 2 predicts `flyescape` conditional on observed `log1p(seedsinit)`. Fixed v13 ridge λ = 1 and existing v13 source feature/standardization definitions are unchanged.

**New augmentation differs from v13 in the training slope, not in the baseline:** In training years only, within each `(Year, location, Type)` group, subtract the group means of width **and response**. Scale width using the **training-year** global SD and estimate one ridge slope β of centered response against centered width (penalty 1) pooled across groups. If a training group has one record, its centered width is zero and cannot identify β. This strips out group mean width differences before estimating β.

For stage 2, additionally remove within-training-group variation explained by **previously observed initiated seeds** from the width column and from the fruit-escape response column with ridge-penalized one-column projections (fixed penalty 1). The stage-2 added width coefficient thus tests **approximate partial** rather than causal/structural independence from measured initiated seeds; it does **not** imply pollen-deposition mediation.

For TEST year, do **not** use that year's response in any mean, scaling, fit or hyperparameter choice. In stage 1, the new regressor is `(width − mean(width | heldout Year, field, Type)) / SD(width | training years)`, which uses **only the heldout cohort's observed widths**. In stage 2, group-center width and **known prior-stage observed initiated seeds** using the same heldout cohort's **features only**, residualize width using training-derived projection, then apply the train-only slope to the old baseline prediction. For response `flyescape`, clamp both baselines and augmented predictions to `[0,1]` as previously.

**Crucial operational limit:** the feature-only TEST-year cohort mean is *transductive*, meaning the investigator has all unlabeled width (and previous-stage seed) observations from the heldout year's group. This is not a forecast for a **single isolated future individual without knowledge of its contemporaneous cohort**, and does not qualify for prospective deployment without an independent pre-measured reference population. The test-year outcome is never read by these transformations.

## Decision and endpoints

For each source year and separately for two response stages, score exactly the same heldout source rows under unchanged baseline and with a within-context width correction. `gain = MSE_baseline − MSE_augmented`. Report all nine fold estimates, test group counts, singletons, fraction in groups of size ≥2, train-only slopes, original-source identity, and independent reproductive-stage response scales. Summaries: (i) original-record-weighted mean MSE gain, (ii) equal-source-year mean gain, (iii) number of folds with positive gain. An **exploratory stringent readiness gate** requires both gains > 0 and positive per-year gains in at least 7/9 years for each mode and endpoint separately. Since v13 outcomes motivated this v14 design, passing is **not independently preregistered confirmation**; failing blocks promotion of a robust within-context effect.

### Mandatory falsifiers

1. Original author file MD5, 437 rows, 367 complete rows, all source missingness, site/type/year counts, and source receipt Git blob SHA must match prior source locks.
2. Entire 9-year heldout folds must contain no training-year row; group response means must use training years only.
3. TEST-year feature centering must not access test responses (`seedsinit` for stage1; `flyescape` for stage2). Stage 2 may use `seedsinit` only as **already observed previous-stage input**.
4. Each TEST-year context-centered width (and stage2 prior input) must sum to zero within groups. No group of size 1 contributes a nonzero deviation.
5. Synthetic **between-context-only** relationships must give zero/approximately zero fitted within-context slope and no invented width gain. Synthetic within-context signal must give positive gain.
6. A test response changed after the prediction is produced must not affect its prediction. Source year group selection cannot use heldout outcomes.
7. All modes and endpoints are included, including negative signs; no posthoc selective outcome dropping or labeling `hyb` and `ten` as two geographical replications.
8. Neither `flyescape` nor initiated seeds should be re-labeled as final viable achenes or genotype inheritance. The source has 41 `seeds > seedsinit` records: an identity such as `seedsinit × flyescape = final viable seeds` is prohibited.
9. Study-specific width prediction is **not the same intervention and response** as v12 *C. arvense* organza bagging; do not pool their MSE or causal effects. Whole years are not independent species/phylogenetic branches.
10. Previous v13 fixture is independently recomputed **before** v14 proceeds. The new v14 fixture is committed only after source-pinned CI has produced a successful original-source report, not inferred from synthetic tests.

## Execution

Python 3.11+ standard library only; no new third-party dependency:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_within_context.py --fetch-pinned --output /tmp/sdea_v14.json
```

The output must retain separate `seed_initiation_log1p` and `fly_escape_fraction` analyses and both location annotation modes. The original 437-row source file is not committed to SDEA. Failure to retrieve the source is `ORIGINAL_SOURCE_UNAVAILABLE`, **not** a zero measured genetic memory effect.
