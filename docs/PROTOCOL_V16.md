# SDEA v16 — is individual genetic identity informative beyond population relatedness?

**2026-10-29. Retrospective, exploratory, independent negative-control challenge to the already-inspected SDEA v15 original Ulex dataset.** This is **not** a new plant field study, a prospective registered randomization, a new gene-function test, or proof of ancestral floral trait regain. The v15 original result stays frozen without alteration.

## 1. Biological question and why the v15 result needs a stronger audit

v15 matched **202 genotyped individuals across six original Ulex parviflorus sites**, with flower-standard petal area/weight, original flower-scar-normalized twig fruits, and SNP-derived relatedness. Its site-heldout M2 positive-relatedness residual smoother added **+0.00001722708** equal-site MSE improvement over the M1 floral-traits model; only **4/6 sites improved**, so the original transfer gate failed. It did not measure heritability or trait recovery.

The new specific hypothesis is: **is this tiny conditional prediction improvement attributable to correctly assigning genomic relatedness to the *individual* mother plants, or could a predictor based on coarse site-wide genetic structure / arbitrary within-site assignments generate the same pattern?**

The source contains no individual plant identifier in the separate **569 pollinator census windows**. This is a partial **G–T–F** present-time prediction, *not* a full G–I–T–F model. There is no reconstructed ancestral pathway, no controlled genotype treatment, and no offspring follow-up.

## 2. Immutable published-original source and locked prior result

Source: Castellanos et al. (2023), *Proceedings of the Royal Society B*, DOI `10.1098/rspb.2023.0141`; data DOI `10.5061/dryad.9zw3r22jp`, Zenodo record **7761289**. See the existing [v15 manifest](../data/v15_source_lock.json) for five original author registered MD5 values.

The runner downloads only the five original files read-only, requires their exact source MD5, reconstructs v15 with `sdea.ulex_partial.analyze`, and **compares the resulting entire report against the frozen [v15 original-source fixture](../results/v15_ulex_partial.json) before entering v16**. If even one numerical value/ID count differs, the v16 experiment stops with a source or old-result integrity error, not a simulated null. The original original-source raw observation files and the 225×225 relatedness matrix are **not copied into SDEA**.

A valid cohort must retain the 202 joined exact `(site, ind)` plants, these original six site labels: `AresM`, `Cheste`, `Chiva`, `Montserrat`, `Simat`, `Sot.Chera`, 225 genotyped source labels and the original source counts and missingness checks. Original `#VALUE!` in `nofruit` stays counted as **missing, not zero**.

## 3. First reference: reproduce v15's true individual genomic assignment

Use exactly v15's whole-original-site holdout. One `ind` appears in only one site; each fold trains the same predeclared M1 ridge model from the other five sites, using only training flower area, dry mass, elevation, and fruit proportions. The baseline for every heldout plant is M1's clipped fruit-fraction prediction. Genomic positive similarity is the existing v15 `max(0, K_ij)` for a heldout source plant `i` and a training-site plant `j`. Denote fixed training residuals `r_j = y_j − M1(x_j)`. The original model is

\[
\widehat y_{i,\mathrm{realG}}=\mathrm{clip}_{[0,1]}\left[ M1(x_i)+ \frac{\sum_{j\in\mathrm{train}}\max(0,K_{ij})r_j}{5+\sum_{j\in\mathrm{train}}\max(0,K_{ij})} \right].
\]

**Every one of the six M1 and real-G MSEs must reproduce exactly (within 1e-11) the already-frozen v15 output.** This prevents changing M1, the reference cohort, the relatedness functional form, or the regression penalty to rescue a result.

## 4. Adverse comparison A — site-averaged relatedness

For a *single heldout site*, average the source-positive `K_ij` **across all unlabeled genotyped test-site individuals i** to obtain one mean nonnegative K weight for each training individual `j`. Use that single vector to compute the training-residual correction for **every heldout individual in the same site**. This removes within-site test-person genomic personalization while preserving test-site ancestry information and the original training residuals. The test-site genotypes can be known (transductive cohort state), but **their reproductive outcomes must not be used**.

Compare site-wise fruit-fraction MSE (1) trait-only M1, (2) site-averaged K, and (3) original individual-K M2. The **identity-specific incremental gain** is `MSE_siteAvgK − MSE_originalIndividualK`. Report all six signed site differences, equal-site and individual-weighted means, and the number of sites with positive identity-specific gain. Note that site-averaged K still retains genomic *site structure*, so this is a *coarse genetic structure comparator*, not no-genome baseline.

## 5. Adverse comparison B — site-preserving genotype identity placebo

Predeclared purely computational negative control: **1,000 site-stratified random bijections**, Python `random.Random(20261029)`. Each replicate randomly reassigns the **original genotyped labels** among eligible phenotype/fruit individuals **inside each original site only**. Every population keeps its exact original set and number of SNP-derived IDs. No source-genomic matrix entry is modified. Source site, fruit outcome, environmental covariate, flora traits, observation grain, whole-site holdout and fixed M1 training residuals **remain unchanged**. The permuted model still predicts test outcomes from K and *training* residuals only.

For each random bijection, calculate the exact v15 M2 minus M1 prediction gain across all six heldout sites. The 1,000-run reference distribution is summarized by min, 2.5th percentile, median, 97.5th percentile, max, and fraction of permuted equal-site gains **at least as large as** the observed true-label gain. Report separate per-site summaries and equal-site / original-individual-weighted aggregation. The full uncollapsed null histories and original plant IDs are **not published**.

**This is NOT a causal genotype treatment permutation p-value.** Exchangeability of genotype labels among plants within one source site is **unproven**; microenvironment, family/genetic assortment and cohort sampling can confound it. The v16 analysis was designed *after* seeing the tiny v15 gain; its rank is an **exploratory numerical placebo stress indicator**, not a pre-registered confirmation.

## 6. Fixed three-part decision rule

Mark `IDENTITY_SPECIFIC_INCREMENT_CLEARS_EXPLORATORY_ADVERSE_CONTROLS` **only** when all hold simultaneously:

1. True-label **equal-site** genomic MSE improvement strictly exceeds the 97.5th percentile of the 1,000 site-preserving label-placebo gains.
2. True-label individual K beats the site-averaged K comparator on both equal-site and individual-weighted MSE.
3. Individual K beats site-averaged K at **at least 5 of 6 source sites**.

Otherwise return `NO_ROBUST_INDIVIDUAL_GENOMIC_IDENTITY_INFORMATION_BEYOND_SITE_STRUCTURE`. **Neither decision** constitutes natural selection, heritability or molecular ancestral floral re-expression. Store full continuous estimates and all negative site comparisons; do not alter the original v15 eligibility, model, penalty, shrinkage, or outcome after inspecting results.

## 7. Negative controls and scientific caveats

- Synthetic no-K-positive-mass cases must return M1 predictions for all three genomic alternatives.
- Identity mapping and site-stratified bijection must retain exact total genotyped labels and same-site membership, one label assigned to one source plant; no arbitrary test-site or across-site shuffling.
- Independent recomputation of the original M2 by its older algorithm must match before any placebo inference. Passing unit tests on synthetic data never substitutes for source authentication.
- Perturbing heldout fruit outcomes after predictions are formed cannot alter predictions. M1 fits and training residuals cannot use fruits from the omitted site.
- The original relatedness matrix can have **negative relatedness values**. Only the original prescribed nonnegative *prediction weights* truncate negatives; raw matrix information is neither rewritten as absence of genotype nor estimated as heritability.
- The plant's fruit fraction comes from source `fruits/scars`, not germinated offspring, viable seed number or lifetime fitness. The 569 pollinator census windows have no `ind` join, and cannot be assigned to genotyped individuals.
- If the real genomic increment is below the within-site placebo or equivalent to site-average K, report an honest **negative identity-specific forecast** even if an earlier published Ulex study independently estimated floral heritability using different methods. Do not claim to refute the original genetic analysis.
- With only six populations in a single study, ordinary row-bootstrapped p-values would overstate independence. This exact placebo is descriptive, and the sites are not independent radiations or evolutionary reversals.

## 8. Reproduction

Python 3.11 standard library only, with original Zenodo source access explicitly requested:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_kinship_identity.py --fetch-pinned --reps 1000 --seed 20261029 --output /tmp/sdea_v16.json
```

The existing GitHub workflow first replays **all earlier v1–v15 fixtures**. The v16 runner independently retrieves the source, rechecks all five original MD5s and the entire frozen v15 output, then computes the new site-stratified null. Only after a successful original-source CI run is a **second** immutable aggregated v16 receipt committed and added to CI for exact numerical re-verification. All earlier source data owners, versions and results are unchanged.
