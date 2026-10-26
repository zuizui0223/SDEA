# SDEA — State-dependent Evolutionary Accessibility

**An auditable, exact synthetic benchmark for integrating ANTO × EGWE × izu-core.**

> **Research status (2026-10-26):** v1 is a *constructed one-step model and inference stress test*, not a fit to natural systems, a replication of the other repositories, or evidence that floral traits have re-evolved. Numerical coefficients and the causal trade-off are explicit assumptions. The three parent programmes and their preregistered results remain independent.

## One scientific question

**When two lineages show the same current visible trait, what hidden present state must be measured to forecast their future expression, genetic accessibility and persistence?** History can inform this forecast when the present state is partially observed; once the future-relevant full state and forcing law are known, it is redundant *under the stipulated Markov model*.

## v1 benchmark: four falsifiable contrasts

1. **Hidden joint toolkit (ANTO):** two fully white populations have identical marginal and pairwise retention of three modules, but different triple retention. Their next-step expression risk differs. Predictive history value disappears on conditioning on the complete current toolkit; an incomplete/noisy assay may not remove it.
2. **Cross-layer alignment (EGWE):** two spatial systems have identical coarse q and toolkit totals but differ in *where* the complete toolkit meets effective ecological interaction. Their next-step joint endpoint differs. When patch q is equal, the difference must be exactly zero.
3. **Persistence vs genetic options (izu-core):** assurance can improve the one-step survival gate while a *declared assumed cost* decreases target allele support and the joint future-access endpoint. With that cost set to zero, assurance cannot reduce target support in this kernel.
4. **Closed vs reopened support:** with zero target allele frequency and no mutation or immigration, finite allele draws cannot spontaneously restore target support; allowing explicitly defined input can reopen it. This is not an ancestral-floral trait regain probability.

A model failure to pass a prespecified negative control invalidates the claimed mechanism in this benchmark. None of these outcomes by itself demonstrates a biological universal law.

## Exact kernel (simplified, not Model 3)

For a current latent module mask \(m\in\{0,\ldots,7\}\), local interaction \(q\in[0,1]\), reproductive assurance \(a\), diploid offspring count \(N\), and target-allele frequency \(p\), the four quantities are:

- First expression \(E=\mathbf{1}(m=7)\alpha q\).
- One-step survival/reproduction gate \(S=q+(1-q)a\).
- Post-selection target frequency \(p_s=pw/(1-p+pw)\), where \(w=\exp[\beta(q-0.5)-ca]\). One-way mutation and migrant replacement are applied if explicitly enabled.
- Probability that at least one of \(2N\) independent allele-copy draws retains the target: \(G=1-(1-p_\mathrm{eff})^{2N}\).

The **one-step joint endpoint** is \(J=E S G\), conditional on three independent gates. It is **not** phenotypic fixation, multi-generation establishment, Mendelian Model 3, or population-level extinction probability. \(a\) here is an *abstract assurance coordinate*, not observed self-compatibility. \(m\) is an abstract functional toolkit, not a verified gene set. Synthetic environment \(q\) is not partner richness.

## Reproduce

Python 3.11+, standard library only:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_benchmarks.py --output results/v1_exact_benchmarks.json
```

Results are in [`results/v1_exact_benchmarks.json`](results/v1_exact_benchmarks.json). Full model decisions, falsifiers and evidence tiers are in [`docs/PROTOCOL_V1.md`](docs/PROTOCOL_V1.md). Source ownership and pinned documents are in [`docs/SOURCE_REGISTRY.md`](docs/SOURCE_REGISTRY.md).

## Boundary / next step

SDEA owns **cross-model definitions, matched-state counterexamples, validation contracts and eventual held-out cross-system tests**. It does not silently copy, recalibrate or reinterpret the parent repositories' simulations. Actual Cirsium molecular/experimental work belongs in [aza3](https://github.com/zuizui0223/aza3); natural fragmentation synthesis belongs in [egwee](https://github.com/zuizui0223/egwee). Next research gate: a common empirical observation design that independently measures functional module compatibility, ecological partner state, finite allelic support and descendant expression across time. Until then, *universal, causal and field-confirmed* claims are **NOT established**.

## v2 — exact multi-generation lineage-state benchmark (2026-10-26)

v2 extends the **one-step v1 fixture without altering it**. A finite Markov recursion now carries joint functional modules, 0..6 target allele copies, two-state ecological interaction, a separate reproductive occupancy gate, and previous expression history through eight generations. The model distinguishes **first expression**, **three consecutively survived expression updates**, **a surviving lineage that ever completed a streak**, and **genetic/functional support still present at horizon**. These are not equivalent to natural establishment or fixation.

For the constructed complete-toolkit case, first expression by H=8 is 0.818665; a 3-consecutive survived run ever is 0.117417; still occupied at H after a 3-run is 0.041242. Matched 2-patch states yield 0.117417 versus 0.008750 for the cumulative 3-run event when only the resource-to-environment alignment is swapped. Different recovery endpoints have **distinct conditional assurance-cost crossover frontiers**, not one universal recovery threshold.

See [the complete v2 methods, predictions, adverse controls and inference ceilings](docs/PROTOCOL_V2.md), [the exact v2 machine readout](results/v2_multigeneration_exact.json), and [`src/sdea/multigeneration.py`](src/sdea/multigeneration.py). Reproduce with:

```bash
python scripts/run_multigeneration.py --output results/v2_multigeneration_exact.json
PYTHONPATH=src python -m unittest discover -s tests -v
```

The kernel is a declared synthetic **independent allele-copy/binomial** approximation, *not* izu-core's Mendelian Model 3; optional expression-triggered toolkit damage is *not* ANTO's empirically measured molecular cost. No measured rates or live plant trait regain are inferred. v2 should be reviewed within the existing Draft PR before any merge.


## v3 — explicit diploid mating-structure bridge (2026-10-26)

**New falsifiable mechanism:** equal present-day phenotype, census and target-allele count are not sufficient to predict next-generation target-locus loss when selfing couples maternal gametes. With N=4 and two target copies, two heterozygotes and one homozygote have equal loss risk 10.0113% under neutral random outcrossing. Under the stipulated q=.75, a=1, δ=.2 mating kernel, their **conditional-on-reproduction** next-generation loss risks differ (11.2304% versus 13.9990%). These are synthetic, not observed genetic or flower-colour outcomes.

The exact one-generation identity is `P(aa)=(1-i)[(1-ū)^2+η Var(u)]`, where η is the viable selfed-offspring fraction and Var(u) is weighted between-parent variance in gamete A-probability. The *conditional* probability of losing all A copies among N recruited diploids is `P(aa)^N`. Consequently the v1–v2 **cost-free assurance leaves allele support unchanged** property is not portable to a mating policy that also changes selfing: assurance may increase the viable-reproduction gate yet reduce conditional genetic support even when direct allele cost c=0. Environmental genotype selection can reverse the risk ordering; none of these signs is universal.

A standard-library exact finite-multinomial generator tracks Mendelian gametes, selfing versus outcrossing, viable seed conditioning, genotype-specific parent weights, phenotype activation, 3-consecutive survived expression, mutation/immigrant inputs and feedback in a *single-locus, constant-N-conditional-on-occupancy* model. It is **not** the full izu-core Model 3, a natural-colour reversal, or a replicated empirical claim. Read [v3 protocol](docs/PROTOCOL_V3.md) and [frozen v3 synthetic result](results/v3_mendelian_exact.json).

Reproduce with `PYTHONPATH=src python -m unittest discover -s tests -v`, `python scripts/run_mendelian.py --output /tmp/sdea_v3.json`, then `python scripts/verify_fixture.py results/v3_mendelian_exact.json /tmp/sdea_v3.json`. v1/v2 results and source-specific constraints are retained without rewriting prior outcomes.


## v4 — pollen versus seed gene flow and dual-route compatibility (2026-10-26)

v4 is a **new and independently labeled conditional synthetic model**. A fixed source patch supplies pollen or donor-produced diploid seed genotypes to a finite receiving patch with two independently assorting Mendelian loci representing candidate ancestral-compatible **A** and alternative-compatible **B** routes. The donor is a stipulated reservoir, not an evolving metapopulation; recipient recruitment is conditional on viable local reproduction and fixed N=2, with no natural phenotype or fitness measurements.

**The key identification counterexample:** Set both recipient diploids to lack route A and both source diploids to carry homozygous A. Compare a one-generation pollen pulse at `m_p=0.4` with a seed pulse at `m_s=0.3250393701`, matched for the chance that each recruit *carries* A. Both generate exactly **28.7696%** generation-one visible-A expression and **54.4428%** generation-one A support; nevertheless A support at generation four after input stops is **26.9981%** (pollen) versus **42.3286%** (seed). Immediate visible A alone cannot identify the movement policy (balanced optimal accuracy **50%**), but a complete genotype census can (**77.2214%**) under this exact hypothesis pair. The genetic-copy inputs are intentionally *not* matched: the experiment reveals why two mechanisms can yield equal photos but different inherited resources.

A second contrast starts with locally retained route A and a B-only donor: pollen can construct dual-capable **A+B offspring in the first generation** (54.4428% chance at the tested pollen value), whereas seed-only import supplies route B but creates **no** initially dual-capable offspring. Later local crosses can create joint capacity under both. Route A/B labels are abstract genetic competence, not measured ancestral pigment pathways, true re-evolution or ecological establishment.

Read [the complete v4 protocol, counterexamples and falsifiers](docs/PROTOCOL_V4.md), [the finite Mendelian source/recipient model](src/sdea/landscape.py), [the frozen machine readout](results/v4_landscape_exact.json), and [24 independent v4 tests](tests/test_landscape.py). Reproduce with:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_landscape.py --output /tmp/sdea_v4.json
python scripts/verify_fixture.py results/v4_landscape_exact.json /tmp/sdea_v4.json
```

**Claim ceiling:** no coupled source dynamics, genetic fixation, fitness-mediated survival, Cirsium regain or calibrated gene-flow rate; phenotype expression has a stipulated environment-dependent law. All prior v1/v2/v3 source contracts and fixtures remain unchanged.
