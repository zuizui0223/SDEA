# SDEA v9 — original natural-data stage and source-admission audit

**Date:** 2026-10-28. **Study class:** retrospective, provenance-pinned qualification and one new **independent raw-data calculation**, not a preregistered natural-evolution test. This version does not change the v1–v8 synthetic fixtures or reassign any published authorship. The biological claim must be gated before model comparisons. The repository contains **no redistributable copy** of the original USDA mother-level CSV or other raw data: execution reads original, independently versioned public repository URLs pinned at Git commits and validates the Git **blob SHA-1 over original bytes** before parsing anything.

## Core scientific question and stopping rule

Can currently available natural records identify **the same lineage's latent developmental route, its evolutionary source, ecological activation, repeated inherited expression, and subsequent reproductive/lineage establishment**? Only a source satisfying *all* parts could support the full SDEA general theory as a prospective causal plant result. Lacking an observed component is **not proof that the biological component is absent**. No score is filled using another study's individuals or an earlier SDEA synthetic variable.

SDEA v9 carries five heterogeneous **convenience evidence units**, deliberately *not* independent taxon replicates: original *Cirsium undulatum* maternal seed routes; Ryukyu *Cirsium* leaf RNA registry; EGWE Honshu–Izu pollination partial state; EGWE *Oenothera* contemporary mating state; and ANTO's non-systematically selected molecular flower-colour cases (including two mechanisms within **one** *Petunia* genus). The sixth pinned source is **izu-core's model-only status**; it is *not* an empirical natural evidence unit.

**Full SDEA trait-reexpression qualification is `NOT_IDENTIFIABLE_FROM_THESE_SOURCES`: 0/5 complete evidence units.** Five here is the size of this deliberately chosen heterogeneous **audit set**, not the denominator of any biological prevalence or cross-clade failure rate. This statement is about data **joint availability** only. It does not negate the original partial ecological and genetic results.

## Direct original-data replay: two anatomical maternal reproductive routes

Original: West & Louda 2021, *American Journal of Botany*, DOI `10.1002/ajb2.1744`, USDA CC0 primary dataset DOI `10.15482/USDA.ADC/1522649`, mirrored in the **aza3** empirical owner repository. Read mother-level `year`, `plt`, `trt`, `plt.gsd`, `ap.hd`, `o.hds` from the source-locked **106 records**. The biological unit is `(year, maternal plant ID)`; **do not link across years** as if they were the same repeatedly measured genetic individual. Do not call the per-head rows or seed counts new independent parents.

For each mother, `A=1(ap.hd > 0)` and `L=1(o.hds > 0)`. Require `ap.hd + o.hds = plt.gsd` and integer nonnegative good-seed values. The source row `(2007, 327)` has `0+1 != 20` and is **excluded from the within-mother paired two-route inference only**. Any additional inconsistent row fails closed rather than being silently discarded. Retain **105 coherent mother/year records** distributed over **8** year×treatment groups.

Report the full four-category table **neither / apical_only / lateral_only / both**, together with the directly identified within-mother probabilities:

\[
\begin{aligned}
P_\mathrm{any}&=P(A\lor L)=P(A)+P(L)-P(A\land L),\\
P_\mathrm{ind}&=P(A)+P(L)-P(A)P(L),\\
\Delta_\mathrm{ind}&=P_\mathrm{any}-P_\mathrm{ind}=-\operatorname{Cov}(A,L),\\
C_\mathrm{obs}&=P_\mathrm{any}-\max(P(A),P(L)).
\end{aligned}
\]

These are **identity/comparison of contemporaneously observed seed-positive anatomical routes**. In particular, the independence predictor is an intentionally inadequate assumption, *not* a causal model that ought to hold. The sign of \(\Delta_\mathrm{ind}\) is observable, but neither its sign nor the two-route gain estimates DNA linkage, source-module retention, adaptive canalization, trait loss/re-expression, pollen deposition, biological complementarity under route removal, or genetically heritable evolution. Both organs already exist in the plants; this is **functional distribution within a generation**, not reappearance of an ancestral organ.

**Audited original 2008 contrasts:** `di` (apical damage × insecticide) has 19 mothers, counts `none=0, apical_only=4, lateral_only=8, both=7`, `P(any)=1`, `Cobs=4/19 ≈0.21053` and `delta_ind≈+0.08864`. In `dni` (apical damage × ambient insect exposure), 12 mothers have `none=5, apical_only=0, lateral_only=7, both=0`, `P(any)=7/12`, `Cobs=0`. In `ndni` the joint seed-success covariance sign is opposite from `di`, and the independence approximation **overpredicts** `P(any)` by ≈4.89 percentage points. None of these is a newly randomized comparison of selection in genotype space; existing species/source analyses remain in aza3, and these 105 mothers are the *same source*.

## Explicit source-reported, NOT independently refitted, empirical layers

- **EGWE Honshu–Izu**: in the pinned 572-row/8-site whole-site heldout task, a functional partial ecological predictor had MSE `1.08773529`; adding distance yielded `1.13209246` (≈4.08% worse). The original archived ecological state is not a fully observed multi-generation genotype state. **SDEA recomputes only the arithmetic of the published locked score.** It does not rerun the original fit or claim new independent support.
- **EGWE Oenothera**: 60 seed-family rows nested in 23 maternal plants; a treatment-only state model has leave-one-mother-out MSE `0.11618623`, adding contemporaneous spatial isolation `0.09186971` (≈20.93% better). This is a **mating proxy/parentage** endpoint, not floral trait re-expression or a measured descendant colour/survival outcome. **SDEA only verifies its locked summary.** The two EGWE results have *different* endpoints; they cannot be pooled as a unified sign on isolation.
- **aza3 Ryukyu Cirsium**: 25 public young-*leaf* RNA-seq libraries, **zero** flower-stage RNA libraries in this registered source and **zero** newly admitted 298 planned Phase-A individuals. No same-individual validated functional colour toolkit, repeated expression or direct future fitness is in the pinned readiness receipt. Absence of a measured floral module is not molecular absence of the module.
- **ANTO**: six selected literature entries; two functionally supported pigmentation regain/re-route examples within **Petunia**, counted as **one genus** not two independent angiosperm replications. These source papers are not ANTO experiments and their trait-specific viable-seed fitness is not qualified. Other rows are phylogenetically inferred, pattern-switch, loss-only, or Cirsium-source metadata, not a random sample of recovery events.
- **izu-core**: original Model 3 remains a registered synthetic intervention programme and is excluded from the natural evidence-unit denominator. SDEA does not import its coefficients and call them plant estimates.

## Strict reproducibility and claim levels

All six input Git commits and blob SHAs are in `data/v9_source_lock.json`. Run with Python 3.11+ standard library only:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_empirical_admission.py --fetch-pinned --output /tmp/sdea_v9.json
```

An optional `--source-root` points to **separately obtained original source files**, at `ROOT/<repo>/<SHA>/<path>`; `--fetch-pinned` is the explicit read-only network option. The default **refuses to pretend local data exist**. A changed or unavailable source fails closed as `SOURCE_UNAVAILABLE` / `source Git blob SHA mismatch` and does not become a no-effect ecological observation. CI prints the precise source-derived JSON to make a frozen result possible without redistributing source raw bytes.

**Interpretation levels:** (1) exact original observed anatomical-reproductive route contingency in real maternal plants; (2) EGWE/ANTO previous source-reported partial ecological/functional/ancestral evidence; (3) SDEA synthetic mathematical possibility (v1–v8); (4) **not yet observed** common same-lineage, multi-generation origin/toolkit/function/fitness test. Level (3) does not upgrade (1)/(2) into level (4). Independent taxa with prospectively frozen holdouts and actual floral genetic pathway assays will be needed for that.

## Falsifiers and negative gates

1. Source contents do not match the **exact pinned raw Git object** => `SOURCE_UNAVAILABLE_OR_DRIFT`, no result imputation or fallback from a study's summary to individual-level raw data.
2. Any additional maternal `apical+lateral != total` mismatch, duplicate `(year,plt)`, non-integer/negative/missing good-seed count, unexpected treatment or year => reject primary same-mother interpretation, not silently drop rows.
3. Coarse independence approximations must be reported alongside **actual joint states** even if delta changes sign or is zero; prohibit a universal positive redundancy inference.
4. If measured functional pigment modules, rooted historical re-expression and descendant fitness are not jointly present *in one lineage*, mark **NOT_IDENTIFIABLE** even if one source paper has colour changes or another source has seed production.
5. No pooling of EGWE individual/fruit outcomes with unrelated Petunia regulatory studies or original USDA mothers; no ecological meta-effect or cross-clade replicate count is estimated.
6. It is a **source/measurement gate** that can change with genuinely new eligible data and correct original IDs; the 0/5 is not a biological zero of regain probability and may be re-audited prospectively with a new version.
