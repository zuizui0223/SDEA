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


## v5 — evolving two-patch metapopulation (2026-10-26)

Unlike v4's static donor, v5 lets **both** finite diploid patches evolve, send/receive pollen and seeds in either direction, lose occupancy and reoccupy through incoming donor-produced seeds. Both candidate genetic routes A and B may be combined within an offspring by independently assorting Mendelian loci; their visible expression is a separate *assumed* q-dependent gate. Source genotype and population availability are no longer fixed, and no external immigrant reservoir exists.

**Exact one-step operator controls:** An initially empty patch cannot be colonized by pollen-only flow (no maternal plant), but can be colonized by viable donor seeds. With two A-type adults in the other patch and the stated synthetic coefficients, the one-step colonization chance at seed-rate 0.35 is **0.458994**. In an A-only recipient and B-only donor contrast, B-seed immigration increases first-step recipient occupancy **0.681481→0.827679** while reducing A genetic support **0.681481→0.538165**; paternal B pollen can instead form an A+B individual in the very first generation (**0.143508** chance of at least one) while donor B seeds alone cannot do so in that same step.

**Exploratory H=7 Monte Carlo:** nine synthetic settings × **9,000 independent entire two-patch landscape replicates** (81,000 histories). A vacant recipient patch is initially colonized at least once in **67.00%** of seed-only trials, versus **0%** under pollen-only, yet only **6.37%** of the seed-only trajectories still occupy that patch at H. Applying an additional 0.32 catastrophe probability to the evolving donor patch reduces initial recipient colonization from 67.00% to **59.86%** and recipient A support at H from 6.37% to **2.44%**. A separate A/B source configuration finds occupied both patches at H in **10.27%** with two-way seed+pollen flow versus **4.10%** with seed alone. Finite-run Monte Carlo sampling intervals are embedded in the JSON; **they do not measure real biological uncertainty**.

**First occupation is not recolonization:** v5 explicitly distinguishes `colonized_initial_vacancy` from `recolonized_after_extinction`. Foreign pollen cannot create a mother where no adult survives, and a gene lost from **both** patches cannot reappear without mutation or an external source (not modeled). Pollen-specific paternal transmission and seed-produced diploid migrants are separately tagged upon entry only: this is **not a full ancestry pedigree**.

[Complete v5 assumptions, falsifiers, exact identities and limitations](docs/PROTOCOL_V5.md) · [v5 finite-state meiosis/flow + Monte Carlo source](src/sdea/metapopulation.py) · [81,000 synthetic-path results and MC-only intervals](results/v5_metapopulation_mc.json).

Reproduce via Python 3.11+, standard library only:
```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_metapopulation.py --replicates 9000 --output /tmp/sdea_v5.json
python scripts/verify_fixture.py results/v5_metapopulation_mc.json /tmp/sdea_v5.json
```

These findings are *kernel-specific counterexamples*, **not natural population recovery, measured A/B floral re-expression, or the same dynamics as the source repos**. Keep the PR Draft for source and ecological validity review. v5 excludes explicit donor seed-production budgets, large landscapes, genotype-dependent selection and real pollinator movement; v1–v4 remain distinct versioned results.


## v6 — finite source seed budgets across three coupled patches (2026-10-26)

**What changes:** Each occupied maternal source has a finite B seed-production **attempt** budget per generation, shared by local recruitment and any of the three destination patches. Paternal pollen location is drawn from a separate 3×3 movement operator. Mendelian offspring occupy at most K=2 adults per destination after optional genotype-weighted **over-capacity recruitment**; catastrophe, first colonization, later recolonization and visible A/B candidate expression are tracked separately. All 3 patches can dynamically become sources or go extinct. This is a new **synthetic** mechanistic counterexample, not a wild-population calibration, molecular ancestor-regain experiment, or source-repository reexecution.

**Exact new result (finite shared source):** with 3 independent attempts, source viable-seed probability S=.755625 and allocations .2/.4/.4, each initially vacant receiver colonizes with probability **.66029688** but both colonize in the same generation with probability **.38245797**, rather than the **.43599197** predicted by multiplying the marginals. The discrepancy is **−.05353400**, because no single viable seed can be delivered twice. At B=1 joint colonization is structurally impossible; at B=8 it is .88823424 vs naive independence .89079203. The exact inclusion–exclusion proof and source-retention formula are in [v6 protocol](docs/PROTOCOL_V6.md).

**Exploratory H=8 ensemble:** 8 distinct parameter cells × **3,000 independent whole-three-patch simulations** (24,000 landscapes). One A-only source with B=3 and all-local seeds survives at H in **20.03%**, vs **4.53%** with strong export (.2 retention/.4 to each receiver), while initially vacant patch1 is colonized at least once in **0% vs 84.17%**. Raising source attempt budget **B=3 to B=8 under the same strong allocation** changes H8 source occupation **4.53%→97.70%** and receiver colonization **84.17%→100%**. These are conditional fecundity/allocation contrasts, not claims that dispersal is universally harmful or that 8 is a natural critical budget.

**Contrary/limited findings:** In an initial A-only / empty / B-only chain, patch1 ever simultaneously hosts one A+B-capable individual in **16.67%** under seed-only chains, **19.13%** under seed+pollen chains, and **26.00%** with fully connected seed+pollen, but adding the stipulated genotype-settlement selection gives an *inconclusive* −0.70 percentage-point effect on patch0 A support (95% MC interval approximately −1.78 to +0.38 pp). These are hypothesized route capacities, **not proven functional molecular re-expression**.

All prior v1–v5 fixtures remain untouched. See [v6 numerical result](results/v6_finite_budget_mc.json), [generative code](src/sdea/finite_budget.py) and [27 new adversarial tests](tests/test_finite_budget.py).

Reproduce (Python 3.11+, standard library only):

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_finite_budget.py --replicates 3000 --output /tmp/sdea_v6.json
python scripts/verify_fixture.py results/v6_finite_budget_mc.json /tmp/sdea_v6.json
```

**Claims not established:** actual source seed budgets, pollinator movement rates, flower-color A/B identity, trait fixation, long-term ecological establishment, natural island effect direction, and independent model-class robustness outside this invented generator. Preserve Draft review pending these gates.


## v7 — matched six-patch network–genotype alignment (2026-10-27)

**The original v6 research gap is now tested directly:** all network policies begin with the same six fully occupied patches, the same A-only versus B-only two-locus diploid genotypes, the same q=.75, the same K=2 and per-occupied-source finite seed-attempt budget B=3, and (within active policies) the same two-neighbor degree, 25% local retention and identical incoming/outgoing flow marginals. The connected clustered and alternating 6-cycles are **graph-isomorphic**; only which genotype-bearing patches are neighbors changes. Two disjoint A-only/B-only triangles form a strict no-cross-route negative control. Source budgets are matched initially and per **occupied** source thereafter; total realized output is an **endogenous demographic mediator** and is not fixed over later generations.

**Analytical generation-one theorem:** The cross-route outcross zygote production probability for a source site is `(1−a(1−q)) q × opposing-route paternal-flow mass`. Independent B attempts yield `P(any A+B seed formed)=1−product_i(1−h_i)^B`. For B=3, this probability is **96.8125%** for clustered, **99.9996%** for alternating and **0** for separated components; expected hybrid seeds are **2.9953 vs 8.9859** (clustered/alternating). These are **viable zygotes produced**, not recruited, observed flowers, genotype fixation, or ecological establishment. The first-generation event is independent of where formed seeds are sent (pollen and seed are distinct operators).

**Monte Carlo H=7 result:** 12 network × operator scenarios, **2,500 whole-six-patch landscape histories per scenario = 30,000**. For **both active seed and pollen flow**, final A+B-capable adult presence is **60.28% clustered vs 62.80% alternating**: paired difference **+2.52 percentage points** with simulation-only 95% interval **[−0.10,+5.14]**, so a directional long-horizon difference is *not resolved*. In seed-only it is 26.88% vs 26.00% (−0.88 pp); in pollen-only 82.92% vs 82.92% (0 pp), both unresolved within the paired intervals. With the two isolated triangles, within-individual A+B capacity is structurally impossible despite separate persistence of each route elsewhere. A late-added 2×2 crossed seed–pollen topology check is **exploratory**; no supported final-access interaction emerged.

**Observation saturation:** varying B in a separate exact mathematical sensitivity shows that a binary 'at least one hybrid seed' indicator saturates while the expected hybrid-seed **number** still differs threefold at fixed B=3. The exact balanced-prior Brier improvement of knowing the connected topology falls from 0.0226945 at B=1 to 0.00025394 at B=3. This is a constructed predictive-state counterexample, not a measured natural-history signal.

Code: [v7 six-patch kernel](src/sdea/network_topology.py) · [v7 protocol, exact proof and falsifiers](docs/PROTOCOL_V7.md) · [v7 frozen machine record](results/v7_network_topology_mc.json) · [27 new adverse/normalization tests](tests/test_network_topology.py). All v1–v6 fixtures are preserved. Reproduce via:
```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_network_topology.py --replicates 2500 --horizon 7 --output /tmp/sdea_v7.json
python scripts/verify_fixture.py results/v7_network_topology_mc.json /tmp/sdea_v7.json
```
**Evidence ceiling:** synthetic fixed-q, no mutation, no source budget reallocation by adult genotype, no empirically validated developmental or genetic pathway, no real seed/pollen observations, and no source-repository independent replication. The observed H7 null/inconclusive results are retained.
