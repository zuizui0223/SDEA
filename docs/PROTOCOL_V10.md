# SDEA v10 — can real same-mother route dependence transfer between years?

**2026-10-28, retrospective and exploratory.** This v10 is a *negative-transfer and measurement-validity challenge* to v9. It is neither a prospective natural-history preregistration nor an additional species, and its selection of endpoints is informed by previously inspected year × treatment counts. The old SDEA v1–v9 models/results retain their distinct estimands and locks.

## Data and identification boundary

The only new predictor evaluation uses the **original** West & Louda *Cirsium undulatum* mother-plant source (USDA CC0 DOI `10.15482/USDA.ADC/1522649`), mirrored in empirical-owner repository `zuizui0223/aza3`. Source pinned to the exact v9 Git commit and Git blob SHA. The v9 parser accepts 106 source mother records, excludes the source-inconsistent `2007|327` (`ap.hd+o.hds=1` versus `plt.gsd=20`), and yields 105 coherent maternal plant-year records (38 in 2007, 67 in 2008). It already established the same-mother four-state counts for eight year × treatment cells. SDEA v10 re-fetches and re-verifies the original **raw bytes**, it does *not* copy maternal records into this repository.

At a given mother-year, define `A=1` if the apical head formed >0 source-defined good seeds, `L=1` if the summed nonapical heads formed >0, and `Y=(A,L)∈{00,10,01,11}`. This is a **contemporaneous anatomical seed-production** outcome, not heritability, genetic/developmental modules, organ-removal counterfactual, plant survival, floral trait reappearance or ancestral return. Treating `year` as an environmental and protocol shift is descriptive; no date-consistent same-genotype individual time series is documented here.

## Precisely matched rival predictors

Within each **training year and treatment group**, count `(n00,n10,n01,n11)`. Add a **fixed 0.5 pseudo-count to each cell**, giving

`p_joint(y)=(n_y+0.5)/(n+2)`.

The alternate `p_independent` uses the **EXACT SAME** training posterior marginal values `pA=p_joint(10)+p_joint(11)` and `pL=p_joint(01)+p_joint(11)`, but enforces independent success outcomes:

`p_independent=( (1-pA)(1-pL), pA(1-pL), (1-pA)pL, pA pL )`.

There is no extra marginal-fit parameter in the joint model. The only different state coordinate is the *within-mother association*. The cell smoothing prevents infinite held-out log loss when a category was absent from the training year. It was chosen for finite-score robustness, not tuned to the held-out outcomes. Under an independent relation with equal cell proportions, the two predictors must coincide.

Primary scoring rules for each untouched held-out mother-year:

- **Multinomial log loss** `-log p(y)` (natural log, lower is better).
- **Unhalved four-category Brier** `sum_k [p(k)-1{y=k}]²` (lower is better).

Primary target: **whole-year transfer in both directions**, train 2007/test all 2008 mothers, then train 2008/test all 2007 mothers, matching only the named treatment groups. No target-year response or plant-level features are used in training; whole-year out-of-sample scores are accumulated across all mothers (n67/n38), and per-treatment scores retained to expose transport heterogeneity. The same USDA study, field years, and treatment designs underpin both folds. These are *not* two independent radiations, species or ecological regimes sampled from a random universe.

**Decision gate**: `JOINT_TRANSFER_GAIN_IN_BOTH_YEARS_AND_BOTH_SCORES` only if **joint beats independent on BOTH scores in BOTH direction folds**; otherwise `NO_REPLICATED_JOINT_FORECAST_GAIN`. Any pooled advantage with sign-changing folds must not be described as an independently validated model. In-sample risk and same-year, same-treatment leave-one-mother-out are **secondary** and cannot rescue a failed year-transfer primary.

## Original association versus predictive transport

For a year-treatment stratum, let `pA=P(A=1)`, `pL=P(L=1)`, `p11=P(A=1,L=1)`. Then **exactly**:

`P(A or L)−[pA+pL−pA*pL] = −[p11−pA*pL] = −Cov(A,L)`.

This is a joint-probability identity, not evidence of a causal competition or compensation mechanism. Cross-year changes in covariance mean that the necessity of representing joint state *for contemporaneous arithmetic* does **not** prove that an estimated association is a stable future transition law.

All source counts, both folds' scores, fitted four-category probabilities, independent and joint marginal identity checks, same-year leave-one-mother-out supplementary results and source hash are exported in `results/v10_real_joint_transfer.json`, then compared in CI after an independent original-source rerun. No multiple testing claim, 95% biological CI, or historical-genetic causality is made from these 105 records.

## Negative controls and fail-closed rules

1. Both predictors **exactly share posterior apical/lateral marginals** under every training count table; perturbing a margin only in one model fails.
2. Perfectly balanced equal four-state training counts make independent and joint predictions identical; artificially disjoint/associated states should differ without changing marginal prevalence.
3. Training/year labels are unchanged across arms; never fit on the held-out year, never treat flowerheads as independent maternal replicates.
4. Missing mandatory original source records, unexpected corruption count, altered Git blob, nonfinite/invalid counts, missing treatment-year cell, altered denominators or unverified natural-source downloads **raise** instead of reporting an ecological null.
5. The independent predicted union error at the *observed* same-year marginals must equal minus empirical `Cov(A,L)` to numerical precision in every original cell.
6. A **negative or mixed** held-out comparison is reported unchanged; never adjust the prior, outcome, horizon, or pooling method to manufacture a universally positive joint-state claim.
7. The source is a single study with two distinct years and a small number of year × treatment strata. Results are not a prevalence estimate of ecological or genetic memory, nor independent support for ANTO ancestral trait reconstruction, izu-core Model 3 or a universal EGWE sufficiency law.

## Reproduction

Python 3.11+, standard library only, SDEA Draft PR:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_empirical_transfer.py --fetch-pinned --output /tmp/sdea_v10.json
python scripts/verify_fixture.py results/v10_real_joint_transfer.json /tmp/sdea_v10.json
```

Source pin: `data/v9_source_lock.json` / `sources.mother_csv`. The original data are never copied into the SDEA git tree. If original source is unavailable or the source SHA changes, stop with `SOURCE_UNAVAILABLE` / `source Git blob SHA mismatch`, **not** `NO_REPLICATED_JOINT_FORECAST_GAIN`.

## Implication for general theory

The **validated narrow proposition** is that the full within-mother joint function state carries information not contained in two observed marginal success rates. The **tested portability proposition** is whether the learned dependence structure generalizes to a second year under the same named treatment; a failed transfer decisively blocks stating that this archive demonstrates a universal, stable joint ecological operator. Testing heritable evolutionary accessibility instead requires verified same-plant floral function and genotypes, historical allele origins, measured contemporary partner environment, and multiple descendant cohorts with fitness endpoints, held out by *lineage or population*, not individual heads.
