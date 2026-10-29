# SDEA v12 — Real multi-site cross-context ecological transfer, separately from genetic regain

**2026-10-29. Retrospective, bounded, original-author data reanalysis.** This is **not** a new empirical campaign, a prospective independent-site preregistration, a genetic floral re-expression test or a natural metapopulation. Previous SDEA v1–v8 numerical constructs and v9–v11 genuine *Cirsium undulatum* maternal route scores remain frozen and scientifically separate.

## Question

After the within-study **2007 ↔ 2008** transport failure in v10/v11, can an entirely different original plant study support *any* cross-context prediction? This time test whether a composite manipulation of floral access has an outcome direction/magnitude that transports between **four sites** within a **single 2022 Cirsium arvense field study**. Also carry a **predator-exclusion control contrast** that is not expected to have a uniform seed-response direction based on previous aza3 source work.

The distinction is important: **predictable treatment-associated seed function is not direct proof of evolutionary accessibility, ancestral trait recovery or functional-module inheritance**. Two source studies of different species and outcomes are not a harmonized meta-analysis.

## Source locks and unit definitions

Original source: Karpate et al. Research Square preprint DOI `10.21203/rs.3.rs-6634172/v1`, repository **Yogita25/Multitrophic_Interactions**, pinned commit `8a151c62ed069f775a16b9bc5fbdce209cbb8691`, path `seeds.csv`, expected Git blob `ca99c762162d15f006f0c73b4603ca71fd320a29`. Archive DOI `10.5281/zenodo.15364152`. This was a 2022 field experiment, not SDEA sampling. The source code retrieves the original UTF-8 bytes by **exact immutable URL**, verifies the SHA using the Git blob hash header, checks required headers, per-head count integer/nonnegativity (allow negative signed `Missing_Seeds`), flowers denominators, and the exact identity `Flowers = Viable_Seeds + NViable_Seeds + Damaged_Seeds + Missing_Seeds`.

Source hierarchy:

- **615** original head records, **612** eligible for the four named sites and four named treatment arms.
- **164** `site × source block × treatment` groups (inferred plant proxies, **not verified unique actual maternal plant IDs**).
- **45** source site-blocks in the eligible dataset; **36** site-blocks have all **four** arms, with counts campus **12**, construction **9**, ditch **8**, rapseed **7**.
- Within an arm/proxy, compute `viable_fraction = sum(Viable_Seeds)/sum(Flowers)` and `emerged_insects_per_head = sum(Insects)/n_sampled_heads`. Never average individual heads before dividing and never count heads or proxies as independent site holdout replicates.
- Within every eligible full block, compute differences `bag (Pollinator Exclusion) − Control` and `Predator Exclusion − Control` on **both** outcomes. Other treatments can be audited for accounting but are not used to select a more favorable main finding.

An older **aza3** source receipt is also pinned (commit `d54aba2dcac3ecbfd4c4fb6bd5b0b88562f3595d`, blob `740bfe83c881ee61269342070253b3baeaf58aa2`) for **numeric crosschecking**, **not** a second independent biological dataset. The previous receipt's source cohort, every site mean of both main measures for both two contrasts, and full 3×3 joint signs are verified; source drift or any mismatch raises a blocking error. Raw original head-level data are not vendored into SDEA.

## External-site prediction design

For each intervention × outcome and heldout site `s`:

1. Fit one constant mean block-paired contrast on **only the other three study sites** (n=36 − n_s complete blocks), 
   \[ \widehat\Delta_{-s}=\frac{1}{N-n_s}\sum_{b:site(b)\ne s}(Y_{b,intervention}-Y_{b,control}). \]
2. Score on the **held-out site mean** `\bar\Delta_s`, using `(\bar\Delta_s - \hat\Delta_{-s})^2` against a **zero-effect reference** `\bar\Delta_s^2`. The primary MSE comparison weights all **four study sites equally**, regardless of block count. The zero reference is a forecasting comparator, **not** a causal no-treatment null.
3. Also report per-site **sign agreement**, each site's squared error, and **secondary** held-out block-level MSE, where block errors are still nested in site and do not generate 36 independent site replicates.
4. The predeclared *within-study* descriptive gate requires (i) bag vs control on viable fraction to have **negative** held-out means in all four sites and improve **all four** site mean squared errors relative to zero, (ii) emergent insect response to have a negative site direction in all four sites, and (iii) predator exclusion viable fraction to fail an all-four-sites improvement gate. This selection of contrasts and gates is **retrospective** and informed by previous aza3 findings; no confirmatory out-of-domain p-value is claimed.

Report all outcomes regardless of sign. Do not pool viable fraction and insects per head into a fake common standardized fitness coefficient; insect emergence is not a direct record of insect pollen delivery or seed predators.

## Original-source result and negative control (site-balanced evaluation)

| Intervention contrast | Outcome | Zero-reference MSE | Other-three-site-predicted MSE | Signs matched | Sites with improved squared error |
|---|---|---:|---:|---:|---:|
| Bag − control | viable seeds / Flowers | 0.2791591645 | **0.0086715082** | **4/4** | **4/4** |
| Bag − control | emerged insects / head | 3.4693509540 | **1.9821044202** | **4/4** | **2/4** |
| Predator exclusion − control | viable seeds / Flowers | 0.0025904215 | **0.0044666989** | **2/4** | **2/4** |
| Predator exclusion − control | emerged insects / head | 1.3839071073 | **0.7861209469** | **4/4** | **2/4** |

At the matched-block grain, bag−control decreases viable seed fractions in **36/36 blocks**. Emerged insects go down in 25/36, remain equal in 8/36 and increase in 3/36; **25/36** have both negative. The bag's seed response is **degenerate in sign**, so the 25 simultaneous-negative cases cannot identify a nontrivial statistical correlation, mediation or antagonism. The predator exclusion contrasts have heterogeneous seed-response signs, so it is a useful **measurement-specific negative portability control** rather than a claim of zero predator effect.

These are four held-out *study sites*, not four independent randomized trials, taxa, radiations or evolutionary origins. The site-level MSE reductions show interpolation/transfer **within one preselected observational source structure**, not prediction on an untouched independent study collected later.

## Key biological boundary

The `Pollinator Exclusion` arm is an **organza bag + wire structure** that can change visiting pollinators, head access by antagonists, humidity, light, microclimate and oviposition. It does **not** isolate pollen flow. `Predator Exclusion` involved adhesive barriers and manual removals, while `Predator Model` was a **plastic spider**, not an independently manipulated functional predator. Site-block-treatment is an inferred plant proxy, not guaranteed true random individual matching. The original viable-seed count is not seed germination or lifetime fitness. The raw signed `Missing_Seeds` is a *mass-balance residual*, not independently observed seed predation. Thus neither path-specific fitness, pollinator compensation, true head morphology/spines, genetic LEGO modules nor actual *Cirsium* ancestral reexpression is identified.

## Fail-closed tests and reproducibility

- Original Git blob byte SHA; no replacements with results-only upstream receipts, and no source metadata URL changes.
- Raw count identity and source denominators, all four site counts and full four-treatment block inclusion.
- Two outcome denominators differ (`Flowers` vs sampled heads); signed negative `Missing_Seeds` is retained, not clipped; impossible negative `Viable_Seeds` rejected.
- Source receipt's site means and joint sign count agree with rederived raw source computations.
- All four heldout sites are excluded **entirely** in their own training; mutate heldout test-site values and confirm its fitted prediction remains unchanged.
- Joint quadrant counting has mutually exclusive outcomes and does not imply same-insect/same-seed linkage.
- A missing source or Git SHA mismatch fails with `SOURCE_UNAVAILABLE_OR_DRIFT`, not a null ecological result.
- No choosing an alternative endpoint after an unfavorable LO(S)O fold; report the predator seed and insect negatives.

The 2022 data establish a real **contemporary outcome-specific and context-limited predictive pattern**; the much stronger SDEA general theory of ancestral, heritable floral route reactivation remains untested by this dataset.

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_site_transfer.py --fetch-pinned --output /tmp/sdea_v12.json
python scripts/verify_fixture.py results/v12_site_transfer_real.json /tmp/sdea_v12.json
```
