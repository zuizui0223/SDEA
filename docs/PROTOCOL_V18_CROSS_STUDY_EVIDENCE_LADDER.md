# SDEA v18 — frozen-receipt cross-study information ladder

**14 November 2026 | retrospective cross-source audit, not a newly sampled plant study, meta-analysis, or independent evolutionary replication.**

This optional, read-only follow-up to SDEA v14–v17 asks a deliberately narrower, falsifiable question: **which observed information layer adds predictive value when the original site's or shelter's entire response is held out, and which apparently positive contribution fails its source-internal identity/context control?** The three original studies are not made exchangeable. All existing model parameters, GenBank/Zenodo source objects, original outcome choices, original observations, versioned receipts, and earlier decisions remain untouched. This branch is stacked on the SDEA v17 Draft branch (unmerged) and must not be confused with main.

## Inputs and original grains

| Original independent study | Model evidence | Prediction endpoint | Original heldout block | Status |
| --- | --- | --- | --- | --- |
| *Ulex parviflorus*, Castellanos et al. (2023) | v15 flower/elevation and v16 relatedness-ID placebo | same-season individual fruit/scar fraction | 6 entire sites, **202 matched individual plants** | original author files source-MD5 locked |
| *Sinapis arvensis*, Höfer et al. (2022) | v17 rainfall/flower supply, shape and measured visitor count per *one plant observation effort* | log1p estimated same-season seeds/plant | 5 whole shelters, **52 plants** | original author files source-MD5 locked |
| *Ipomopsis*, Campbell et al. (2022) | v14 original v13 width forecast challenged by within-year/context features | log1p initiated seeds | 9 whole years, **367 rows without validated stable plant ID** | same original study as v13, **not** a new independent source |

Reused **aggregate** JSON receipts only: `results/v14_within_context.json`, `v15_ulex_partial.json`, `v16_ulex_genomic_placebo.json`, `v17_sinapis_process.json`. The source individual rows are not downloaded, recopied or refitted here. The 1,000 genomic and 200 visitor placebo remappings are *source-internal sensitivity references*, not 1,200 independent plant experiments.

## Exact retrospective computation

For each original site/shelter leave-one-block-out MSE pair, define `gain_b = MSE(baseline)_b − MSE(candidate)_b`. Report the equal-block average, original-record-weighted average, `mean(gain_b)/mean(MSE_baseline_b)` as the **ratio of equal-block means**, the mean of per-block relative gains (not the same estimator), positive-block count, and the minimum/maximum average after omitting one existing block's gain.

**Omitting one existing fold here is an influence calculation on already-produced folds**, NOT a new cross-validation, new prospective heldout fit, cluster bootstrap, effect confidence interval or p-value. Training samples overlap between original leave-one-site/shelter-out folds. We never combine different study responses into one shared MSE or one mean percentage.

### Verified new arithmetic from the original v14–v17 frozen aggregates

| Within-study conditional comparison | Equal-block MSE gain | Improvement as ratio of equal-block baseline means | Positive blocks | Delete-one-block minimum gain |
| --- | ---: | ---: | ---: | ---: |
| Ulex elevation → + flower area/weight | **+0.00102078** | **+11.69%** | 5/6 | +0.00060571 |
| Ulex flower → + individual genomic relatedness smoother | **+0.00001723** | **+0.22%** | 4/6 | +0.00000186 |
| Sinapis rainfall/flowers → + flower shape | **−0.02303447** | **−9.77%** | 1/5 | −0.03278167 |
| Sinapis rainfall/flowers/shape → + observed plant visitation | **+0.06312187** | **+24.39%** | 4/5 | +0.04803031 |
| Sinapis rainfall/flowers → + observed plant visitation | **+0.05800206** | **+24.60%** | 5/5 | +0.04667131 |

The four Sinapis models share the same plant IDs and holdout folds; the Ulex comparisons share the same site folds. **The Ulex fruit-fraction and Sinapis log-seed MSE scales are incompatible**, and their treatment/predictor contrasts are different. The numbers must not be ranked as biological strengths, pooled as a meta-analytic effect or tested as independent studies.

Source-internal adverse evidence preserved:
* **Ulex** true individual-K gain +0.000017227 is **below** the within-site genotype-ID placebo 97.5th percentile +0.000048235; 312/1,000 identity reassignment gains were at least as large. No robust genomic *individual-identity* gain beyond site structure is admitted, despite the small positive original MSE.
* **Sinapis** visitation M1→M2 gain +0.0631219 exceeds the within-shelter × rainfall-treatment plant-visit-label placebo 97.5th percentile +0.0115707; 0/200 reassignments matched it. This is **observational contemporaneous predictive information**, not experimental pollen deposition/visitation causality.
* **Ipomopsis** previous v13 original source-wide width gain +0.1674166 row-weighted became **−0.0202821** in the stricter within-year/site/type v14 model (6/9 years with a positive within-context gain). This challenges transfer and causal interpretation, not the existence of flower-width biology.

## Scientific interpretation and what is not identified

The **same visible plant trait is not an adequate universal summary** of the ecological process or genealogical state. Observed process can add incremental forecast information, and apparent genetic incremental performance can fail an identity-level control. However these separate original studies **do not identify** a single-plant joint `genotype × interaction × phenotype × future fitness` relationship; neither pollen delivery, parent-offspring transmission nor future re-expression of an ancestral flower colour is measured along the whole chain. In Ulex there is no pollinator-visit record joined by original individual plant ID; in Sinapis individual molecular genotypes are absent. Ipomopsis cannot supply that missing coupled observation.

**Therefore:** this is evidence toward a measurement/observability framework and an empirical test-design target, not confirmation of cross-study genetic evolutionary accessibility. We do not use ANTO v26 Petunieae gene-to-pigment species correlations as if they were the same individuals as these field studies, and do not numerically pool EGWE/izu-core synthetic effects with empirical observations.

## Falsification / fail-closed checks

`src/sdea/evidence_ladder.py` verifies the exact 202/52/367 source units, unique site and shelter IDs, all records in holdout denominators, finite nonnegative MSE, consistency with independently frozen source-level gains, genomic and visitation placebo reference limits, prior v13→v14 context reversal, and strict lack of the joint biological source key. Invalid or drifted source receipts stop instead of silently promoting an imputed zero or a positive score. `tests/test_evidence_ladder.py` includes receipt mutation controls. The code generates only aggregate summaries.

Run:

```bash
PYTHONPATH=src python -m unittest tests/test_evidence_ladder.py -v
python scripts/run_evidence_ladder.py
```

**Submission/merge state:** documentation and retrospective derived arithmetic only; existing draft remains unmerged. No new preregistered primary biological hypothesis, field evidence or meta-analytic estimate. A true statistical synthesis requires independent studies measuring a comparable outcome and baseline contrast, comparable source-sampling grains and an admissible uncertainty model; none is manufactured here.
