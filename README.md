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


## v7b — directed seed reachability after identical first censuses (2026-10-27)

**Complementary to v7, not a replacement.** The six-patch v7 compares *degree-balanced, isomorphic* network placement with matched incoming/outgoing totals; its exact first-hybrid-zygote contrast is positive but many seven-generation retention comparisons are inconclusive. **v7b** instead isolates a deliberately different property: the **presence/absence of a directed seed path** in a three-patch system. Its cycle versus backflow graphs have identical **per-source seed budget** B=3, local retention 60%, export 40%, initial A-only genotype distribution, q=.75 in all patches/years, selfing, germline model and local pollen routing. They **do not** preserve the same incoming column sums or all graph degrees; these are part of the deliberately varied routing intervention.

**Sharp exact witness:** At generation zero, only patch0 has two A-only diploids; patches1 and2 are vacant. The cycle 0→1→2→0 and backflow 0↔1,2→0 have identical source0 and (empty) source2 routing rows, so their **full generation-1 genotype/population state laws are identical**, with P(patch1 colonized)=**0.6602968801** and P(patch2 colonized)=0. At generation2 the forward cycle has **P(patch2 occupied)=0.1904815097**, whereas backflow has **exactly zero** because no initial source can reach patch2 by seeds. The analytic binomial formula and an independent categorical enumeration agree. Pollen-only mobility cannot create a new mother in the empty unreachable patch.

**Predictive-state consequence:** With the full post-generation-1 genotype state known but future seed routing hidden, the exact binary t2 patch2-occupancy Bayes Brier risk is **0.0701901781**. Including the operator reduces it to **0.0451396014**; the exact gain is **0.0250505768**. This is **missing future forcing/transition law**, not non-Markov ancestral memory. A same-operator negative control has Brier gain 0.

**Exploratory seven-generation validation:** Eight separately analyzed scenario cells × **4,000 independently seeded whole-three-patch histories** = 32,000 synthetic landscapes. From the one-A-only source, patch2 is ever colonized in **45.15%** under cycle and **0%** under backflow. In a *different* initial family, with A-only patch0 and B-only patch2, the event of a joint A+B-capable individual ever appearing in patch1 is **43.20% cycle** and **46.90% backflow** (independent Monte Carlo difference **+3.70 points**, simulation-only 95% interval +1.52 to +5.88 points). The balanced network is an exploratory third treatment: it changes direct incoming flow to patch1 and must not be described as matching the first-census control. The symmetry control (two identical A sources, cycle vs reversed cycle) gave −0.25 pp for patch1 ever-colonized with a simulation interval crossing zero.

**Falsifiers and provenance:** [precise v7b protocol](docs/PROTOCOL_V7B.md), [independent restricted generator](src/sdea/network_reachability.py), [27 new tests](tests/test_network_reachability.py) and [v7b deterministic numerical record](results/v7b_network_reachability_mc.json). Runs in the *same Draft PR* but does not alter v7 or older fixtures:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_network_reachability.py --replicates 4000 --output /tmp/sdea_v7b.json
python scripts/verify_fixture.py results/v7b_network_reachability_mc.json /tmp/sdea_v7b.json
```

**Inference limit:** mathematical and finite synthetic probability only, not empirical gene-flow, actual ancestral trait restoration, universal network superiority, independent cross-taxon confirmation or source-repository replication. Realized total viable seed production varies endogenously after different colonization trajectories even though each occupied source has a fixed budget.


## v8 — structural genetic access versus selection-dependent establishment (2026-10-28)

**New versioned synthetic test** retaining all prior v1–v7b protocols and fixtures. The v8 model imports the existing SDEA v7 six-patch Network, degree-balanced clustered/alternating/components graphs, two-locus Mendelian founders, K=2 recruited adults per patch, B=3 seed attempts per occupied source and fixed q=.75. Only the post-zygotic relative recruitment weight w of viable A+B seeds is changed: 0 lethal, 1/3 competitively disfavored, 1 neutral, and 3 competitively favored. These are constructed *selection policies*, not molecular or field estimates.

**Exact stage separation:** For a fixed seed pool [A+B,A,B] with K=2, P(A+B recruited) = 1−2/[(w+1)(w+2)], giving **0%, 35.7143%, 66.6667%, 90.0000%** across the four weights. Initial Mendelian viable A+B seed production is selection-invariant: **0%** for disconnected A/B triangles, **96.8125%** for clustered and **99.9996%** for alternating connected graphs. A genotype that is not formed cannot be made by selection at recruitment.

**Executed H=7 factorial:** 12 policies × 2,000 independent whole-six-patch synthetic histories (24,000), with paired PRNG keys and separate production/recruitment streams. H7 joint adult presence: clustered **0 / 54.40 / 62.00 / 68.05%**; alternating **0 / 54.70 / 64.30 / 69.70%** under lethal/disfavored/neutral/favored policies. Favored−neutral H7 difference **+6.05 pp** clustered (paired simulation 95% +4.78 to +7.32), **+5.40 pp** alternating (+4.18 to +6.62). Nonlethal selection changes genetic retention while the any-patch-occupied H7 indicator remains **unchanged** within a network under this genotype-independent reproduction kernel. Hard lethal exclusion also lowers occupancy.

**Important negative:** Alternating versus clustered H7 joint adult retention remains statistically unresolved under either neutral (difference +2.30 pp, paired simulation 95% −0.58 to +5.18) or favored (+1.65 pp, −1.15 to +4.45). Large exact differences at generation-one zygote formation do not guarantee a long-horizon directional advantage.

[v8 protocol and falsifiers](docs/PROTOCOL_V8.md) · [reproduction/selection source](src/sdea/selection_filter.py) · [full frozen model result](results/v8_selection_filter_mc.json) · [new adverse tests](tests/test_selection_filter.py).

Run Python 3.11+ standard library: 
- PYTHONPATH=src python -m unittest discover -s tests -v
- python scripts/run_selection_filter.py --replicates 2000 --horizon 7 --seed 1182026 --output /tmp/sdea_v8.json
- python scripts/verify_fixture.py results/v8_selection_filter_mc.json /tmp/sdea_v8.json

**Evidence ceiling:** no independent source-repo simulation reproduction, measured plant natural selection, actual trait re-evolution, pollen genetics calibration, experimentally established A/B pathways or long-run lineage establishment. Per-occupied-source budgets match ex ante; total realized descendant seed output is endogenous. Keep Draft for mechanistic and cross-species validity review.


## v9 — pinned real-data evidence admission (2026-10-28)

**First direct original-source replay, not another synthetic evolutionary result.** A new standard-library source-pinned auditor retrieves the original *Cirsium undulatum* CC0 USDA mother-plant CSV from [aza3](https://github.com/zuizui0223/aza3) at an immutable commit. It verifies original Git blob bytes before reading 106 maternal records, reports the sole incoherent `2007|327` record, then separately calculates **105 coherent same-mother paired apical versus lateral good-seed functions** in eight year×treatment strata. In the 2008 damaged/insecticide cell (n=19) the observed counts are **apical-only 4, lateral-only 8, both 7, neither 0**. In the 2008 damaged/ambient cell (n=12), **7 lateral-only, 5 neither**. These are **contemporaneous anatomical reproductive routes**, not ancestral flower trait re-expression, genetic redundancy, experimental organ-removal benefits or inherited trait recovery.

A second *source-reported-only* receipt gate checks [EGWE](https://github.com/zuizui0223/egwe) Honshu–Izu 8-site ecological source scores and 23-mother *Oenothera* mating scores, the [ANTO](https://github.com/zuizui0223/anto) nonrandom Petunia/cross-clade mechanism registry, [aza3](https://github.com/zuizui0223/aza3) Ryukyu **young-leaf only** RNA source readiness, and [izu-core](https://github.com/zuizui0223/izu-core) as an explicitly model-only research source. The joint prospective ancestral floral-genotype/history/expression/descendant fitness contract is **not identifiable in these sources**. `0/5` is a descriptive tally across five deliberately heterogeneous evidence units, **not** a wild-plant event prevalence.

Original biological source records remain in their original repositories; only derived **reproducibility code, SHA locks and bounded results** enter SDEA. See [v9 protocol, source ownership and falsifiers](docs/PROTOCOL_V9.md), [six exact immutable source blob pins](data/v9_source_lock.json), [auditor](src/sdea/empirical_admission.py), [offline gate tests](tests/test_empirical_admission.py), [reproduction script](scripts/run_empirical_admission.py).

Run with network **explicitly enabled**:
```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_empirical_admission.py --fetch-pinned --output /tmp/sdea_v9.json
```

If source access or checksum validation fails, the result is `SOURCE_UNAVAILABLE_OR_DRIFT`, **never** a null result or an invented zero. EGWE and ANTO scores are source reports; SDEA does **not** refit those original datasets. This stage is empirical *admission/partial-stage confirmation only*, **not** a validation of the entire SDEA evolutionary theory.

**Externally rerun original data admitted:** [frozen v9 source-checked result](results/v9_empirical_admission.json), independently regenerated at the pinned Git blobs by [GitHub Actions run #40610725397](https://github.com/zuizui0223/SDEA/actions/runs/40610725397). It contains only bounded derived summaries (not redistributed plant-level raw CSVs). Read `full_contract_decision=NOT_IDENTIFIABLE_FROM_THESE_SOURCES` as a **measurement limitation**, not a biological absence-of-reexpression result.


## v10 — external-year transfer of same-mother reproductive-route association (2026-10-28)

**Retrospective, exploratory REAL ORIGINAL-DATA prediction test; not a new simulation, not independently sampled species.** SDEA v10 reuses v9's SHA-locked USDA CC0 source from [aza3](https://github.com/zuizui0223/aza3), with **106** original `Cirsium undulatum` mother-year rows, excluding solely the **original inconsistent 2007|327** mother, leaving **105** valid mother-years in 8 year×treatment cells. Raw plant-level rows remain in the empirical-owner repository and are never copied into SDEA.

**Strict information contrast:** Predict 4 mutually exclusive mother-year good-seed states (neither, apical-only, lateral-only, both). Fit a symmetric per-cell pseudocount **0.5** on the **training year only**. The joint categorical model and its independence/product-marginals control have **exactly identical fitted apical and lateral success marginals**, differing solely in retained within-mother covariance. Score untouched held-out mothers by natural-log 4-category loss and **unhalved** categorical Brier. Both directions use the same treatments: **2007→2008 (67 mothers)**, **2008→2007 (38 mothers)**. Within-year leave-one-mother-out is transparently secondary.

**Crucial negative transfer:** joint-minus-independent per-mother log loss / Brier = **+0.057425 / +0.009949** for 2007→2008 and **−0.117906 / +0.000913** for 2008→2007. The prespecified joint-transfer gate fails because joint information **does not improve both metrics in both held-out directions**. Pooling opposite signs does not rescue it. In no-damage/ambient (`ndni`), within-mother apical–lateral good-seed covariance changes sign between source years (**2007 −0.04938, 2008 +0.04889**), a direct descriptive illustration of unstable ecological association. Two years from the *same published study*, not replicated lineage-level experiments, so no universal year effect or random-sample inference is justified.

Read [v10 protocol](docs/PROTOCOL_V10.md), [source-locked pure-standard-library predictive algorithm](src/sdea/empirical_transfer.py), [23 adverse/negative-control unit tests](tests/test_empirical_transfer.py). The code reads original bytes through v9 Git SHA enforcement; the previous v1–v9 numerical fixtures remain unchanged. The v10 run/report is only a test of **contemporaneous anatomical maternal reproductive routes**, not heritable joint molecular modules, causal removal/replacement, ancestral phenotype recovery or long-run ecological establishment.

Reproduce with Python 3.11+:
```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_empirical_transfer.py --fetch-pinned --output /tmp/sdea_v10.json
```

**Do not promote joint-state marginal descriptive adequacy to cross-year predictability.** Unmeasured pollinator/seed-predator and year/protocol differences may explain failure; they are not experimentally decomposed. The full evolutionary accessibility chain remains empirically unidentified.

The original pinned USDA data were re-read independently by GitHub Actions: [frozen v10 real-data cross-year forecasting receipt](results/v10_real_joint_transfer.json). This records the **negative transfer** as a primary output. Numerical source provenance: [original source CI execution #40614715787](https://github.com/zuizui0223/SDEA/actions/runs/40614715787). No original maternal row-level data are redistributed.


## v11 — decomposition of year-shifted maternal reproductive routes (2026-10-28)

**Posthoc natural-source transport diagnostic, not a new original biological cohort or another synthetic evolutionary result.** An original-data runner SHA-verifies the same West & Louda USDA CC0 maternal CSV mirrored by empirical owner [aza3](https://github.com/zuizui0223/aza3) and freezes the original 106/105 mother-year cohort split (38 mothers in 2007, 67 in 2008, eight year×treatment cells, source anomaly 2007|327 excluded only from joint status). The raw mother-level CSV is never copied into SDEA.

**Exact year-shift decomposition:** For the four-category same-mother good-seed distribution, denote `Q_year` the product of the two observed anatomical route marginals and `c_year` their binary covariance. The equality `P2008-P2007 = (Q2008-Q2007)+(c2008-c2007)*(+1,-1,-1,+1)` holds exactly. In no-apical-damage/ambient-herbivory treatment (`ndni`), both apical and lateral success probabilities fall from **7/9 to 2/15** and covariance changes **−0.049383→+0.048889**. Full four-state total variation across the years = **0.8000**, marginal-product shift TV = **0.7017**, association-vector TV norm = **0.1965**. These norm magnitudes cannot be added and they do **not** identify causal contribution percentages.

**Train-year-only shrinkage counterfactual:** Fit a single fixed four-cell association weight `lambda` from the predeclared grid (0,.1,.25,.5,.75,1) by training-year leave-one-mother-out **log loss**. All candidate predictions share exactly the original v10 smoothed apical/lateral Bernoulli marginals; only within-mother dependence differs. 2007 training selects **lambda=.25**; in untouched 2008 test data (67 mother-years), the held-out loss change relative to independence is **+0.008410 log loss / +0.001466 categorical Brier** (both worsen). 2008 training selects **lambda=.75**; on untouched 2007 (38 mothers), differences are **−0.103753 / −0.001242** (both improve). The two-direction improvement gate is **not met**, so joint association transfer remains unreplicated. These are retrospective choices made *after viewing v10*, not independent prospective model-selection evidence.

This is a **genuine mixed/negative real-data result** about **contemporaneous anatomical seed-production routes**, not ancestral pigment re-expression, genetic module redundancy, causally established year mechanism, mating system evolution, heritable offspring, or observed natural ecological establishment. Fréchet bounds express what the two individual route margins alone can and cannot identify; transportability is a separate question.

[v11 source-pinned protocol](docs/PROTOCOL_V11.md) · [numerical and exact-identity source](src/sdea/empirical_drift.py) · [30 new adversarial tests](tests/test_empirical_drift.py).

Reproduce (requires original source fetch authorization):
```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_empirical_drift.py --fetch-pinned --output /tmp/sdea_v11.json
```

All v1–v10 sources and fixtures remain unchanged. A source checksum mismatch or absent original data **fails closed** rather than creating a fake null.

**Executed source and frozen negative record:** [v11 source-verified drift and train-only shrinkage fixture](results/v11_real_transfer_drift.json), generated from the original pinned USDA Git blob by [GitHub Actions original-source run #40618487013](https://github.com/zuizui0223/SDEA/actions/runs/40618487013), no source rows copied into SDEA. The fixture stores actual exact signed distribution-shift vectors and the complete candidate lambda sensitivity; do not overread them as causal ecology/selection or cross-clade gene restoration.


## v12 — four-site real ecological response transport (2026-10-29)

**New, genuinely independent original biological source from v9–v11, but *not* independent field observations conducted by SDEA:** 2022 *Cirsium arvense* heads from Karpate et al. [original author source commit](https://github.com/Yogita25/Multitrophic_Interactions/blob/8a151c62ed069f775a16b9bc5fbdce209cbb8691/seeds.csv), exact Git blob SHA `ca99c762162d15f006f0c73b4603ca71fd320a29`. SHA-gated 615 original head records → 612 eligible → 164 **site×source-block×treatment plant proxies** → 45 source site blocks → 36 complete **four-arm matched site blocks** (campus 12, construction 9, ditch 8, rapseed 7). The earlier independently implemented [aza3 calculation](https://github.com/zuizui0223/aza3/blob/d54aba2dcac3ecbfd4c4fb6bd5b0b88562f3595d/data/evidence/cirsium_arvense_four_arm_joint_response_independent_js_v1.json) is pinned and checked **only as an arithmetic cross-validation of the same original 2022 study**; it is NOT a second species/site-sampling cohort.

**Entire-site leave-one-site-out (LOSO) external-context forecast:** For each one of the four heldout sites, predict the within-block **treatment−Control** effect from the other three sites, then compare **4 equally weighted heldout SITE MEANS** against a zero-effect forecasting baseline. Primary viable-seed-per-flower proportion for compound Pollinator Exclusion *bagging*: site-balanced MSE `0.279159→0.008672`, **4/4 negative heldout site means**, **4/4 sites improve on the zero predictor**, and original complete blocks **36/36 negative**. A separate insect-emergence-per-head response: MSE `3.469351→1.982104`, **4/4 directional sign agreement**, yet only **2/4** sites improve squared-error magnitude (strong site heterogeneity). Same-block seed-and-insect decreases occur in 25/36; there are 8 zero-insect-change and 3 insect-increase blocks. **Because seed-change sign is negative in all 36, this joint-sign count does not by itself identify a nontrivial response-covariance or mediation effect.**

**Crucial adverse control:** Under Predator Exclusion−Control, the viable-fraction mean varies in sign across sites; site-balanced MSE **worsens** `0.002590→0.004467`, sign agreement **2/4**. Predator-exclusion insect emergence direction agrees 4/4, but its magnitude improves in only 2/4. Hence predictive transport is **outcome- and intervention-specific**, not a universal floral-pathway law. No difference in p-value is translated into a genetic or causal phenotype return.

**Scope:** A site×block×treatment row is an inferred plant proxy, not a verified mother plant ID or independently randomized treatment; an organza pollinator-exclusion bag also changes other insect access and microclimate. Insects emerged cannot be assigned one-to-one to seed predation. Viable seed classification is not germination or lifetime fitness. Four sites are **one** study in one field year, not four independent biological replications, and this analysis is retrospective with previous aza3 original-data result inspected first. No heritable genetic route, flower-color reversal or natural long-term establishment is established.

[Full v12 scientific protocol](docs/PROTOCOL_V12.md), [read-only pinned source lock](data/v12_source_lock.json), [source-pinned four-site original analysis code](src/sdea/site_transfer.py), [28 negative/adversarial tests](tests/test_site_transfer.py). The previous v1–v11 results are unchanged.

Reproduce (standard library, network only when explicitly enabled):
```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_site_transfer.py --fetch-pinned --output /tmp/sdea_v12.json
```

Frozen source report will be included as an additional fixture after the original-source runner has succeeded in independent GitHub CI. The entire cross-layer ancestral flower-history/trait-restoration/fitness sequence remains **NOT IDENTIFIED**.

**Source-replayed evidence receipt:** [v12 frozen genuine original-head held-out-site results](results/v12_site_transfer_real.json). Reproduced from the author-pinned original 615-head source in [CI source run #40807391636](https://github.com/zuizui0223/SDEA/actions/runs/40807391636). The frozen source readout contains only aggregate site fold metrics, no sampled individual head or plant raw records. Missing/changed originals invalidate CI rather than being imputed as ecological zero.


## v13 — independent Ipomopsis stage-specific year-heldout forecast (2026-10-29)

**New independently sourced original published botanical study — but not an independently randomized replication of the Cirsium v12 bagging contrast.** Source: Campbell et al. 2022 American Naturalist, [DOI 10.1086/716740](https://doi.org/10.1086/716740); [original open dataset 10.7280/D1KM49](https://zenodo.org/records/5083399), published file `Selection_by_pollinators_and_predators.csv`, **author MD5 `12ed3dbe6d269c5f808f9c6b230d468c`**. A [previous original-file audit from aza3](https://github.com/zuizui0223/aza3/blob/d54aba2dcac3ecbfd4c4fb6bd5b0b88562f3595d/data/evidence/ipomopsis_two_stage_trait_selection_independent_original_js_v1.json) is separately checked by Git blob SHA; this is the **same Ipomopsis study** as the prior audit and not a second new experiment. The published cohort has 437 original rows from 2010–2018 with no stable plant individual ID; 367 rows have paired width, initiated-seed and fruit-fly escape measures.

**New falsifiable predictive estimand, distinct from prior correlation audits:** Hold out **one whole original source year** (9 year folds). In each fold, test whether adding *measured corolla width* to a fixed type+physical-site baseline improves heldout original-row MSE for (a) `log1p(seedsinit)` and separately for (b) `flyescape` (fraction of fruits escaping fruit-fly attack), whose baseline *also conditions on observed initiated-seed log count*. Fixed ridge penalty 1, all scaling learned on training years only. Track MSE improvement by year, equal-year and original-row-weighted, with an exploratory strict 7/9-year readiness gate per endpoint. **Do not pool the two endpoint scores**, infer direct causal effects from conditional stage regressions, or treat years/rows as independent species or ecological origins.

**Cross-source STOP:** v12's treatment is **physical organza Pollinator Exclusion vs Control**, with viable seeds/original Flowers and head insect emergence. Ipomopsis has no matched bagging intervention or equivalent denominator. Therefore **the numeric Cirsium v12 bag effect CANNOT be directly transported, pooled or independently validated by this study**. A scientifically valid outcome here is only a second botanical original-study *stage-specific forecasting* assessment.

[v13 source/interpretation protocol](docs/PROTOCOL_V13.md) · [source immutable receipt and original MD5](data/v13_source_lock.json) · [new read-only source and held-out modeling code](src/sdea/external_stage_transfer.py) · [29 new fail-closed tests](tests/test_external_stage_transfer.py). Source identity mismatch causes a CI failure; no synthetic zero is substituted. The source raw original row file is **not copied into SDEA**. v1–v12 model and original-source fixtures stay untouched.

Reproduce on network-authorized CI:
```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_external_stage_transfer.py --fetch-pinned --output /tmp/sdea_v13.json
```

The v13 machine readout becomes a versioned frozen fixture only after the original registered MD5 and earlier aza3 receipt SHA are verified on an independent GitHub Actions run. No genetic ancestral memory, homologous morphology or trait re-evolution is empirically tested by this accession.

**v13 source-verified result:** [frozen genuine Ipomopsis 9-year heldout-stage readout](results/v13_ipomopsis_stage_transfer.json), originally regenerated from registered MD5 plus aza3 Git SHA by [GitHub Actions source run #40828380864](https://github.com/zuizui0223/SDEA/actions/runs/40828380864). In 367 complete source rows, adding width improved heldout-year initiated-seed prediction in **7/9 years** (row-weighted MSE gain **+0.167417**, equal-year gain **+0.134617**) but fruit-fly escape fraction prediction in only **5/9 years** (equal-year gain **−0.0000672**, despite a tiny row-weighted gain). Under this explicitly exploratory rule, the predictive width contribution is **stage-specific**. The result is noncausal; it is not a replicated Cirsium bagging treatment effect or molecular ancestry test.

 
## v14 — within-year field/type conditional forecasting (2026-10-29)

**A retrospective strict check of the already analysed Campbell et al. 2022 Ipomopsis original source, NOT a new field cohort or independent experiment.** Original author Zenodo CSV Selection_by_pollinators_and_predators.csv remains published-MD5-locked (12ed3dbe6d269c5f808f9c6b230d468c), and the prior aza3 study audit remains git-blob SHA-locked. All 437 original rows, 367 complete measured width/seed-initiation/fruit-escape rows, nine original years and original missingness are source-checked. The previously frozen SDEA v13 result must first reproduce exactly from the same input bytes. No source row-level data are committed.

**Precise newly tested hypothesis:** The positive v13 heldout-year width prediction for initiated seeds may merely reflect between-Year, between-Type or between-field composition. Estimate a width slope from TRAINING data centered WITHIN each original Year × physical field × original Type. Compare prediction with a fixed v13 Type+physical-field baseline on an entire heldout year. The heldout feature itself is width centered within its Year × field × Type **using only the test cohort's unlabeled width inputs**, not its outcome. A separate stricter source Site-code × Year × Type version is sensitivity only: hyb and ten are code labels from the same upper physical field, not independent localities.

The two original reproductive-stage endpoints are fit separately: initiation is log1p(seedsinit); fruit-fly escape is a FRUIT fraction, whose baseline already conditions on measured initiated seeds. The stage-two conditional width coordinate is partialled against previous-stage initiation differences *inside training context groups* before estimating the width slope. Source-year-heldout prediction always uses only the training response; stage-two test initiated seeds are allowed as already observed previous-stage input. All source years, negative results, context coverage and original endpoint scales are retained. Fixed ridge penalty 1, no tuning on the heldout-year outcomes.

**Important operational limitation:** group means estimated from all **unlabeled** heldout-year widths constitute **transductive cohort-conditional prediction**, not a model for a solitary new plant without a premeasured comparison cohort. The exploratory 7-of-9-year readiness gate is NOT preregistered natural evidence; this design is a post-v13 challenge, not independent confirmation. No raw phenotype/genome ancestry, trait reappearance, actual seed viability or natural selection is identified. Numerical Cirsium v12 organza-bag intervention effects are incomparable to observational Ipomopsis width/fruit-escape outcomes.

[Full v14 protocol](docs/PROTOCOL_V14.md) · [MD5/SHA and failure-mode manifest](data/v14_source_lock.json) · [source-dependent and outcome-firewalled model](src/sdea/within_context.py) · [32 new adverse tests](tests/test_within_context.py).

Reproduction (Python 3.11+, public source access required): PYTHONPATH=src python -m unittest discover -s tests -v ; python scripts/run_within_context.py --fetch-pinned --output /tmp/sdea_v14.json

The empirical v14 numerical fixture will be frozen only after successful independent source-pinned CI, without rewriting v1–v13 or translating synthetic model effects into empirical data.
