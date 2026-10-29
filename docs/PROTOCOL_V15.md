# SDEA v15 — a genuinely independent original-source plant-level G–T–F partial test

**2026-10-29; retrospective, exploratory and fully source-fingerprinted.** This version does not modify or supersede v1–v14. The user-defined scientific aim is to connect individual genotypic similarity (G), flower traits (T) and a realized reproductive measurement (F) in an original natural plant source, subject to strict external-forecasting tests. It does not try to force the Ipomopsis v13/v14 coefficient or *Cirsium* v12 bagging response onto a non-equivalent new plant system.

## 1. External original source and immutable authentication

Castellanos et al. (2023), *Proceedings of the Royal Society B*, DOI [10.1098/rspb.2023.0141](https://doi.org/10.1098/rspb.2023.0141), deposited original source DOI [10.5061/dryad.9zw3r22jp](https://doi.org/10.5061/dryad.9zw3r22jp), Zenodo record [7761289](https://zenodo.org/records/7761289) for *Ulex parviflorus*. The authors collected individual-level floral measurements, same-plant twig fruit counts, SNP-derived pairwise relatedness, and **separate patch-level pollinator censuses** across six field sites. This is a different original field study, plant family, biological response and experimental design from Ipomopsis and *Cirsium*. It is not a replicated causal SDEA functional-module restoration experiment.

Source file MD5 checks are fixed **before** source biological rows are opened:

| Source file | Original deposited MD5 | Data role |
|---|---|---|
| `fruits.and.mean.floral.traits.csv` | `8ca72b7d9721f6d6279e793a1526c43e` | Plant ID, location, elevation, mean flower area/weight, flower scars/fruits |
| `floral.traits.csv` | `14216bb46af39bfd0b959c1aa092485b` | Flower replicate plant-ID audit |
| `relatedness.matrix.csv` | `ecaa79680ffa25cb45ebdb73c11fe407` | 225 source-reported individual SNP-relatedness entries |
| `pollinator.census.csv` | `78494af547de11c1800af50a74c5642f` | Patch visit effort and insects; **no individual ID** |
| `README.md` | `9a6d62c907a93fe7b67268922f756995` | Source column semantics, especially fruit/failed-flower accounting |

Fetch only through the original **Zenodo record API**; require exact record ID, file inventory and source MD5 before entering the numerical analysis. No raw plant, genotype or pollinator rows may be committed to SDEA. The prior [EGWE source-grain receipt](https://github.com/zuizui0223/egwe/blob/main/artifacts/empirical/ulex_zenodo_source_gate_20261008.json) looked at headers only, **not** original outcome values or an SDEA forecast. It is not independent biological evidence; it provides an external schema benchmark.

## 2. Original grain and strict source gates

- Fruit record key must be `(site, ind)`, unique within original fruit table. The table fields are exactly `(site,ind,elev,weight,area,scars,fruits,nofruit)` as observed by EGWE source schema audit.
- `scars` is the number of source-recorded original flower scars, the denominator of the authors' fruit-set proxy. Require nonnegative integer `fruits` and `nofruit`, positive scars, the identity `nofruit = scars − fruits`, and `0 ≤ fruits ≤ scars`. All invalid original-count inconsistencies cause **BLOCK**, not selective omission (except declared missing-value rows). Fruit set equals `fruits/scars` and is **not lifetime total seed production**, descendant viability or confirmed pollination.
- Mean floral trait coordinates are original `area` (standard-petal area, mm²) and `weight` (flower dry weight, mg), not Ipomopsis corolla width and not Cirsium head angle/spine. `elev` is a covariate, not a randomized ecological change.
- Individual genotypes must come through the `ind` row/column labels of the 225×225 relatedness matrix. Negative pairwise relatedness entries are valid and **must never be silently rewritten as DNA absence**; for one deliberately bounded smoothing operator only their **positive part** is used as a prediction weight. Matrix symmetry, distinct labels and >=90% available numeric pair cells are necessary.
- At least 150 original fruit/phenotype records must link to original SNP-relatedness IDs; at least six source populations and at least eight joined plants each. Reusing an individual ID in different source sites blocks the join; **do not guess fuzzy aliases or count unmatched plants as genetically absent**.
- `floral.traits.csv` has repeated flowers within `ind`; do not treat flower replicates as independent plants. Report trait flower↔fruit key coverage as an audit. The patch-level `pollinator.census.csv` **does not have `ind`**, so do **not** link its observation windows or flower visits to individual plants or infer unmeasured pollen deposition, despite the source paper's natural-history context.
- Source counts and overlap can be discovered prospectively within this retrospective source-audited v15 run. Failure of an eligible source join returns **SOURCE_ACCESS_OR_GRAIN_STOP**; there is no null genetic effect if the prerequisites are missing.

## 3. Frozen six-site LOSO noncausal prediction comparisons

**Predict** contemporary individual twig fruit fraction, tested on one **entire original site** held out at a time (six folds), with each matched plant scored once. All training transformations and model coefficients use the other five sites only. For a source plant `i`, `y_i = fruits_i/scars_i`.

- **M0 / ecological baseline:** intercept and training-standardized `elev`, ridge penalty 5. The intercept is unpenalized; the elevation coefficient is penalized. Squared-error on the observed fruit proportion, equal weight per source individual.
- **M1 / added contemporary floral structure:** M0 + training-standardized `area`, `weight`, and both centered quadratics. Same ridge penalty 5, same complete genotyped individuals and same whole-site holdout. This permits a bounded curvilinear trait response but **does not replicate the original paper's full stabilizing-selection analysis**.
- **M2 / new genotype-similarity layer:** M1 + `sum_j max(0,K_ij) [y_j−M1_fit(x_j)] / [5 + sum_j max(0,K_ij)]` over **training-site individuals** `j` only. This deliberately simple nonnegative-relatedness residual smoother has a fixed shrinkage constant 5; no test fruit result, future individual ID, population-specific random intercept, hyperparameter tuning or source-label mutation is allowed. Genetic relatedness is transductively measured for test plants; without their genotypes, M2 is not operational. **This is NOT SNP-derived heritability, functionally causal genotype, molecular ancestry or genetic evolutionary ability.** The residual smoother is heuristic and exploratory; any positive result must be independently challenged with true historical/pedigree and genotype confounding controls.
- Clamp all predictions to `[0,1]`. Compute six **site-level mean squared errors** and the mean of the six site gains for M1–M0 and M2–M1 separately. Also report individual-weighted gain and the number of sites with improvement. A strict “all six improved” descriptive marker is printed but is **not** a hypothesis p-value, a registered cross-clade transport theorem or a license to optimize model parameters after outcomes.
- The six sites are six populations in **one study**, not six independent publications. Genotyping and fruit measurement are contemporaneous; this analysis is **not** a multigeneration future-state prediction.

## 4. Scientific adverse controls and interpretation ceiling

1. Altering a source byte or a required MD5 → stop, never impute null. Changed original ID, duplicate fruit key, reused individual across sites, missing genotype/phenotype alignment, nonsymmetric or mostly missing relatedness → stop.
2. Training and heldout sites may not overlap; heldout fruit-set values cannot influence fit, source cohort standardization, kinship residuals, or tuning.
3. Artificial variation of the test fruit response must leave predictions unchanged. Genetic-similarity shrinkage tending to infinity recovers the M1 prediction; zero/negative similarity weights do not create a trait or genomic effect by construction.
4. A genotype-based residual prediction is not proof of heritability; a positive genotype–fruit prediction may reflect regional genetic structure, spatial/environmental confounding, or sample selection.
5. The patch-level pollinator census is quarantined from individual `ind` joins. Its observation/flower effort is not an individual pollinator–genotype interaction assay. No `site`-only join of census visits to focal plant genotypes is permitted as a same-plant variable.
6. Do not compare MSE magnitudes with v13/14 `log1p(seedsinit)` or fruit-fly escape, and do not directly transport the v12 bagging treatment effect: they are **different phenotype measures, reproductive endpoints, and designs**.
7. Original paper prior findings (area stabilizing selection and weight heritability) remain prior scientific work; v15's incremental whole-site forecasting question is a different, explicitly **posthoc** estimand and neither discovers their genetic mechanism nor overrides their analysis.
8. Reaching G–T–F partial admission does not imply functional ancestral-pigment pathway memory, effective pollinator visits or genetically inherited phenotype regain. Full ANTO×EGWE×izu-core field validation remains **not identifiable**.

## 5. Reproduction and next decision

Python 3.11+, standard library only:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_ulex_partial.py --fetch-pinned --output /tmp/sdea_v15.json
```

First source-authenticated CI run prints full **aggregate-only** JSON. If data are accessible and the original join gate succeeds, freeze its exact typed output as `results/v15_ulex_partial.json`, add an independent source reread/replay comparison, then rerun the full SDEA suite. If a gate fails, freeze a **bounded source-access/identifiability-stop receipt** rather than changing the input or fabricating positive results. The original per-plant CSV, relatedness matrix and pollinator census stay only in the external Zenodo archive.
