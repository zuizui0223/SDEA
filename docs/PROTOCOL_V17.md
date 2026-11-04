# SDEA v17 — is direct plant-level visitor state predictive beyond flowering supply across sheltered rainfall contexts?

**Frozen analysis design: 2026-11-05. Original independent plant data, retrospective exploratory reanalysis; NOT new field samples or a genetic trait-regain experiment.** Earlier SDEA v1–v16 numeric results and all negative controls are unchanged. This is deliberately a different source with actual same-plant floral, visitor and reproduction keys, rather than re-tuning the Ulex genomic-similarity model.

## 1. Original source and novelty boundary

Höfer, Lindner, Ayasse and Kuppler (2022), *Reduced seed set under water deficit is driven mainly by reduced flower numbers and not by changes in flower visitations and pollination*. Original open data DOI [10.5061/dryad.66t1g1k5f](https://doi.org/10.5061/dryad.66t1g1k5f), [Zenodo original record 7335086](https://zenodo.org/records/7335086). This experiment manipulated water supply (Control, PrecipRed, Drought) for natural field-raised *Sinapis arvensis*, with original plants nested in five rain shelters and treatment plots. **The source authors already reported lower flower and seed supply under drought and limited mediation through altered floral attractiveness**; SDEA does NOT claim that biological result as novel. Nor does the dataset contain measured genotype alleles, descendants, restored ancestral pigments or sufficient actual pollination deposits.

Source has original fields:

- `Plantlist.csv`: `plantID`, `motherID`, `treatment`, `shelter`, `plotID` — **one plant**; experimental maternal family labels are not molecular genotypes.
- `Floral_traits.csv`: same `plantID`, `motherID`, rainfall treatment, `shelterID`, plot, daily/total source flowering, floral morphology, original author estimated total `seedsperplant` and other floral/reproductive traits — **one plant**, contemporaneous seasonal measurements.
- `Int_data.csv`: `plantID`, `visitorgroup`, same original family, rainfall treatment, shelter and plot, `NrInteraction`, `observetime` — **multiple visitor-group rows for each plant**. Critically, `observetime` is *a plant's full observation effort repeated once per visitor-group row*; summing it across groups would spuriously lower visit intensity. Sum the `NrInteraction` counts across groups but divide by `observetime` **once**, then express visits per plant-observation hour. This is total observed visits, not observed successful pollination or proven pollen source.
- `README.md`: original field definitions and exposure semantics.

Read-only download uses the original Zenodo API with exact record ID and four fixed original MD5 checks from [the source manifest](../data/v17_source_lock.json). The source is expected to contain **53 original plant-list rows, 53 floral rows, 427 visitor-group rows**. The source and experimental `plantID`, `motherID`, `treatment`, `shelter`, `plotID` must match exactly across the three observation tables. Unknown numeric tokens, invalid nonnegative seed counts or zero observation time, duplicate identifiers, unexpected shelter/treatment and source digest mismatch are **data-quality stops** and must NOT be converted to an ecological zero. Declare missing continuous floral/source fields separately and exclude only corresponding complete-case plants. Require ≥35 complete original plants and ≥5 matched plants from each of the five shelters, or block the biological model. **No original individual source CSV rows are committed to SDEA**.

## 2. Single strictly fixed predictive experiment

The primary response is `log1p(seedsperplant)` from the original source. `seedsperplant` is an author's estimated same-season seed quantity, not measured daughter recruitment, viable seeds proven by germination, genetic trait transmission or plant lifetime fitness. We use **five leave-entire-rain-shelter-out folds**, with one heldout prediction per eligible individual and equal-shelter as well as original-plant-weighted proper **squared prediction error**. Since there are only five experimental shelters in one source study, this does not amount to independent multi-experiment replication.

All models share identical eligible plants, precipitation labels, fixed `ridge=4` coefficient penalty except intercept, and training-only standardization. No withheld reproductive response is used to fit, select, standardize or permute predictors.

- **M0: plant flowering supply and rain treatment:** rainfall treatment categorical labels + `log1p(sumflower)`. `sumflower` is a **same-season** flower count; controlling for it is a **conditional contemporaneous forecast**, not a causal direct rainfall effect (it is downstream of rainfall manipulation).
- **M1: floral morphology:** M0 plus original `meandisplayflower` and `meanpetalwidth`. No genetic or environmental causal identification.
- **M2: direct observed interaction layer:** M1 plus `log1p(sum(original NrInteraction across visitor groups) / (single original observetime in hours))`. The encounter layer is measured on the SAME plant. No claims that *each* recorded taxon transferred pollen.
- **M3: visit-only alternative:** M0 plus same observed visit intensity. This is a secondary, label-agnostic process adequacy diagnostic.

Primary process increment = `MSE(M1)−MSE(M2)`; morphology increment = `MSE(M0)−MSE(M1)`. Report **each of five shelter-specific signed gains**, the equal-shelter gain and individual-weighted gain; *never* omit a failing shelter or replace the ecological endpoint with a visually appealing success score.

## 3. Source-internal visitor-identity adverse control

After source schema and model are frozen, run **200 deterministic, one-to-one shuffles** of the observed plant visit-intensity coordinate **within each source `(shelter,treatment)` stratum** (RNG seed 20261105). Maintain entire original plant/floral/reproductive rows, rainfall treatments, actual plots, original observation effort/guild count lists, and site holdout membership unchanged. Recompute the same fivefold M2 forecast after each label permutation. Report reference equal-shelter process-gain median, 2.5/97.5 percentiles, and fraction >= real correspondence. This is a **noncausal descriptive placebo**, **NOT a valid genetic/random-assignment p-value**, because visit rates are not randomized among plants conditional on these limited labels and source maternal background, phenology and microhabitat remain possible confounders.

**Strict *exploratory* readiness gate** for this chosen predictive family: original M2 vs M1 must improve the mean in both weighting schemes, at least four of five shelters, and exceed the 97.5th percentile of this frozen placebo reference. Every failure is preserved. A passed gate establishes at most **source-conditional same-plant process predictiveness under a coarse block holdout**, never causal pollinator mediation, genomic evolution, universal cross-plant process sufficiency or a newly identified ancestral route.

## 4. Mandatory negative tests / scientific boundaries

1. Original registered 4-file MD5/Zenodo record/file-name identity, source counts, header and 5-shelter overlap gates.
2. Plant ID plus maternal family, rainfall treatment, original shelter and plot must agree exactly across tables; incomplete interaction row records are **missing**, never `0` visits.
3. `observetime` must agree across multiple insect-group rows per source plant; sum visits but NOT observation effort. Treat visitor taxon identity as observed groups, not separate plants or independent effect estimates.
4. Regression M0–M3 shares original complete-case plants, fixed ridge and train-only transformations; whole heldout shelters cannot leak into the training equations. Test seed outcomes may be perturbed after fitting but must not modify predicted values.
5. When original visit intensity has no variation, the M2 process increment should be zero to numerical tolerance. When a known synthetic visit effect exists, M2 predictions change and a label-shuffle sensitivity is non-inert.
6. Unmeasured true pollen deposition, individual genotype and descendant floral expression remain **NOT IDENTIFIABLE**. Mother-ID is maternal family, not a molecular or independent evolutionary ancestral trait module. Total seed supply and flower counts are contemporaneous rather than a prospective multigeneration series.
7. Original field study's rainfall treatment grouping is the authors' experimental design; the **M2 increase** in predictive accuracy cannot be promoted to a rainfall-treatment intervention on bee behavior, species-specific pollen delivery or a genotype×visitor fitness interaction.
8. No numerical pooling with `Cirsium` mesh-bag interventions, `Ipomopsis` observed width seed-stage regression, `Ulex` individual SNP relatedness. All have distinct treatments, phenotypes, biology and denominators. A sequence of separate source-specific MSEs is not a meta-analysis.
9. A selected source response is *source author estimated seed output*, not success of heritable trait reacquisition. All v1–v16 synthetic and prior independent-study negative results remain explicitly intact.

## 5. Reproduce

Standard-library Python 3.11+:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_sinapis_process.py --fetch-pinned --reps 200 --seed 20261105 --output /tmp/sdea_v17.json
```

The original source report is printed in a single machine-readable prefixed JSON record in successful CI, excluding original plant ID and raw visits. Only after actual published source download/replay succeeds is its result frozen as a versioned fixture and byte-checked on the next run. Access/source-schema failure is an epistemic STOP, not a null ecological result. The SDEA Draft PR remains open pending independent field, genotypic and descendant-layer validation.

## Original-source outcome after fixed v17 design (2026-11-05)

Original read-only four-file author MD5 authenticated by [GitHub Actions run #41923722762](https://github.com/zuizui0223/SDEA/actions/runs/41923722762): 53 original plant IDs and 427 nested visitor-group records, one missing original seedsperplant, **52 eligible exactly linked plants**, 5 rain shelters (11,9,8,12,12). There is no individual genomic genotype or descendant trait axis; motherID is a family marker, not an assay.

**Heldout log seed output (not flower or seed survival):** Adding measured floral shape to rainfall+total season flower supply M0 worsens mean predictive MSE **−0.0230344721** (1/5 shelters improves); direct total observed plant visitor intensity added after morphology M1 improves equal-shelter MSE **+0.0631218749**, individual-weighted **+0.0623222860**, and **4/5** shelters improve. Adding interaction intensity directly after M0 rather than after morphology improves **all 5 original shelter fold means** (equal-shelter **+0.0580020583**), but this model is separately measured. Shelter 4 is the negative M1→M2 fold **−0.0248047028**; never omit.

**Within-shelter and original rainfall-treatment 200-permutation observational placebo:** Reassigning plant-specific visitation labels within each site/treatment shows median predictive increment **−0.0094199242**, 97.5th percentile **+0.0115706959**; **0/200** equals or exceeds true source assignment **+0.0631218749**. The fixed exploratory readiness gate passes. The placebo ranking is NOT a randomized ecological causal p-value, and the source original author already showed rainfall-treatment flower/seed supply changes; this new analysis cannot attribute the residual predictive contribution to successful pollen delivery alone.

All exact original-source aggregate fold results, source validity, plus both the positive and negative interpretation boundaries, are frozen in [v17 original-data receipt](../results/v17_sinapis_process.json). The prior SDEA v1–v16 fixtures stay unchanged. Full ancestral molecular genotype, actual individual source pollen movement and next-generation offspring fate are still NOT IDENTIFIABLE in this original archive.
