# SDEA v8 — structural accessibility versus selection-dependent realization

**2026-10-28 JST. Status: new synthetic mechanism test; not natural plant evidence.** This v8 experiment was designed as an additive robustness test of SDEA v7's finite diploid network model, **not** as a replacement for its frozen numerical result or as an independent reexecution of ANTO, EGWE or izu-core. The v7b three-patch directed-reachability experiment addresses a different estimand and remains unmodified. Nothing in the present project establishes ancestral floral pigment restoration, genotype-specific field fitness or a universal fitness parameter.

## Question and three separate biological objects

When genetic crossing across patches **can** create an A+B-capable diploid, must the resulting viable zygote be **recruited**, and must its genetic competence still be present after seven generations?

A viable dual-capacity seed *formed* is neither an adult that survived recruitment nor a viable persistent phenotype. This model distinguishes:

1. **Structural reachability:** does the pollen graph have cross-route edges, so a viable seed carrying at least one A allele and one B allele has positive formation probability under the founders and mating law?
2. **Post-zygotic recruitment:** holding the seed pool fixed, what is the probability at least one A+B seed is among K adult recruits, under the given selection operator?
3. **Finite-horizon retention:** after selection has changed inherited genotypes, can any A+B diploid still be found at horizon H=7? Also report **global coexistence of A and B in different individuals** and **any patch remaining occupied** separately.

All A and B labels mean *abstract biallelic route-associated genetic capacities*, not an assayed ancestor, an actually functional biochemical combination, flowers of specific colours, or an established fitness advantage.

## Matched model and selection-only intervention

Reuse SDEA v7's existing `network_topology.Network`, `routing`, `INIT`, `_cross` and `_first_hybrid_exact` source. Six patches initially hold two homozygous A-only adults in each of patches 0–2 and two B-only adults in each of patches 3–5. All share the same q=0.75, K=2, B=3 seed production attempts **per occupied maternal patch**, autonomous assurance 0.45 and selfed viability 0.80. Two unlinked biallelic loci undergo Mendelian meiosis. No mutation, immigration, genotype-dependent parental fecundity, mating preferences, real trait feedback or measured fitness is added.

For each connected configuration, the seed-destination and pollen-father matrices are the *same* degree-balanced 6-cycle, either `clustered` or `alternating`. Both have the same row/column sums, 25% seed retention, 37.5% routing to each neighbor, 18 initial seed attempts, site count, q and all initial genotypes. The graph labels are isomorphic as unlabelled graphs: the contrast concerns **genotype-to-network alignment**. The two disjoint triangles `components` form a no-cross-route negative control (not connected to the other two, but they retain the same 2-regular degrees and row/column marginals).

**Stage 1 — produce seeds:** Reproduce v7's self versus outcross attempts, viable-seed gate, paternal location, independent Mendelian gametes, and one unreplicated seed destination. A positive or negative *adult recruitment fitness* must **never** change the actual seed-production kernel in the same generation. The code uses separate PRNG streams for seed production and settlement. Identical source states plus identical production seed keys yield byte-identical seed pools across all selection policies at generation 1.

**Stage 2 — recruit adults:** Give each joint-capable A+B seed the relative *settlement weight* w; give every single-route/nonjoint seed weight one. Pick at most K=2 adults without replacement by successive normalized weights. The four fixed policies are:

| Policy | Joint recruit weight `w` | Interpretation in this synthetic operator |
|---|---:|---|
| `lethal` | 0 | hard post-zygotic exclusion: joint seeds produced but cannot become recruits |
| `disfavored` | 1/3 | lower relative priority **under crowding** |
| `neutral` | 1 | genotype-blind finite recruitment |
| `favored` | 3 | higher relative priority **under crowding** |

A strictly positive w does **not** add a mortality filter when enough adult slots are empty; a seed survives and is recruited if admissible and there is room. The w=0 case is qualitatively different: it is a **hard inviability/establishment gate** rather than the limit of a finite competitive advantage among otherwise admissible candidates. Thus it may also lower demographic occupancy. No genotype-specific parental fecundity or absolute physiological fitness is estimated.

**Stage 3 — propagate:** Under identical population and environmental initial states, update six finite populations for H=7 generations. Later source genotypes and occupied patch counts may vary endogenously after selected recruitment; later total realized seed production is **not** held constant and is not a separately identified selection-mediation coefficient. Track production, recruitment and surviving genetic support as separate events. No time-invariant repeated visible phenotype or natural establishment is inferred.

### Exact conditional recruitment identity

For a specified three-seed pool containing one `(1,1)` joint A+B seed, one pure A seed `(2,0)` and one pure B seed `(0,2)`, with K=2 recruit slots, sequential weighted sampling yields

\[
P(\text{joint seed recruited}\mid\text{that pool},w)
 = 1-\frac{2}{(w+1)(w+2)}.
\]

The missed-joint events are the two disjoint paths A→B and B→A. The formula is independently verified by recursive enumeration of disjoint ordered recruit sequences and stochastic draws.

| Joint weight | Exact joint recruitment chance, fixed pool |
|---:|---:|
| 0 | 0% |
| 1/3 | 35.7143% |
| 1 | 66.6667% |
| 3 | 90.0000% |

**This is an exact conditional mathematical example**, not the probability a randomly produced landscape contains this particular three-seed pool. A genotype with absent B cannot be given B by changing recruitment fitness. No selection intensity can manufacture a genotype missing from the germline production support.

The same fixed-pool forecasting task under 50:50 prior on `w=1/3` versus `w=3` gives Bayes-optimal Brier improvement `(0.9−5/14)^2/4 = 0.0736734694` by revealing the **future settlement operator**. This is omitted-forcing information, not intrinsic ancestral evolutionary memory or an empirical predictive-model validation.

### Source-exact first-generation genetic formation under the original v7 generator

For the frozen v7 homozygous founders and any of the above settlement-selection policies, let `c_i` be pollen migration mass from source i to a patch carrying the opposite candidate route. The exact probability of **at least one viable A+B seed produced** is

\[
P(\text{at least one hybrid zygote produced})
=1-\prod_{i=0}^{5}[1-(1-a(1-q))q\,c_i]^B.
\]

At B=3, `components` has zero probability, `clustered` has **0.9681251027**, `alternating` **0.9999960766**. The exact expected number of *produced* viable joint seeds is **0 / 2.9953125 / 8.9859375**, respectively. **These probabilities are invariant to w** because selection is solely downstream of this production process. This is not a model for natural functional/fitness establishment.

## Fully executed synthetic sensitivity: 12 cells × 2,000 whole-landscape replicates

The complete factorial is **3 network alignments × 4 selection policies = 12 scenarios**, 2,000 independently seeded *complete six-patch landscape trajectories* each, all at H=7: **24,000 synthetic histories**. Pair the same PRNG seed index between scenarios, with a **separate seeded reproduction and recruitment stream per generation**. Paired Monte Carlo (MC) uncertainty is computed from the variance of within-seed 0/1 differences; it is **not** a natural sampling confidence interval or evidence of taxon-level replication. The first-generation analytic calculations are exact, and v8 preserves the separately frozen v7 and v7b outputs.

| Same seed+pollen network | Settlement weight | P(any gen-1 joint seed produced) | P(any gen-1 joint adult recruited) | P(joint adult present H7) | P(any patch occupied H7) |
|---|---:|---:|---:|---:|---:|
| Clustered | 0 | 96.95% | 0% | 0% | 31.70% |
| Clustered | 1/3 | 96.95% | 85.75% | 54.40% | 85.00% |
| Clustered | 1 | 96.95% | 91.50% | 62.00% | 85.00% |
| Clustered | 3 | 96.95% | 95.40% | 68.05% | 85.00% |
| Alternating | 0 | 100.00% | 0% | 0% | 8.75% |
| Alternating | 1/3 | 100.00% | 100.00% | 54.70% | 84.15% |
| Alternating | 1 | 100.00% | 100.00% | 64.30% | 84.15% |
| Alternating | 3 | 100.00% | 100.00% | 69.70% | 84.15% |

Rows reflect **model-generated MC proportions**, while the analytic production probabilities above apply exactly. A finite sample can show 100.00% when the true event chance is less than one. All four `components` selection regimes produce **zero joint seeds and zero joint recruits throughout the horizon** because no A/B gene contact exists across the isolated subnetworks.

### Predeclared paired model comparisons

- Within clustered, H7 joint-adult chance: `favored − neutral` **+6.05 percentage points** (paired MC 95% +4.78 to +7.32); `disfavored − neutral` **−7.60 pp** (−9.10 to −6.10).
- Within alternating: `favored − neutral` **+5.40 pp** (+4.18 to +6.62); `disfavored − neutral` **−9.60 pp** (−11.11 to −8.09).
- In *both connected* graphs, the nonlethal selection treatments change the **joint genotype outcome** while the `any patch occupied H7` indicator has exactly the same value for every paired random key in this implementation. This is a structural result of genotype-independent seed-production viability and positive relative recruitment weights; it does not show genetically neutral demography in natural populations. The `lethal` policy removes seedlings and disrupts that occupancy invariance.
- Under `neutral`, alternating − clustered H7 joint-adult difference is **+2.30 pp** (paired MC 95% −0.58 to +5.18); under `favored`, **+1.65 pp** (−1.15 to +4.45). Neither H7 connected-graph directional contrast is resolved by these simulations, **despite the large exact first-generation production contrast**. Under `disfavored` the difference is +0.30 pp (−2.75 to +3.35).
- First-generation seed production contrasts are exactly identical for all four w policies with each initial matched network and seed key. They cannot be used as new independent evidence of selection effects.

**Scientific conclusion bounded to this model:** Network genotype-contact changes the **chance of constructing** a combined trait route; downstream ecological recruitment weights change its **chance of being retained**, including eliminating adult establishment at w=0 without retroactively preventing its formation. The structural network signal at generation one and the net H7 signal are distinct estimands. Selection is not identified as explaining the v7 H7 topology null, nor are all evolutionary-genetic routes fixed over longer horizons.

## Falsifiers, negative controls and inferential firewalls

1. **No retrocausality:** with identical initial states and production seed keys, w cannot change first-generation seed pools or first-generation hybrid seed counts. Test paired equality for each initial seed and each topology; a violation invalidates the claimed separation.
2. **No emergent genetic modules:** genotype-dependent adult fitness may reweight available recruitment candidates, but cannot produce an A+B seed if all lineages carry A-only or B-only in disconnected genetic components without any pollen/seed bridge or mutation. Zero across seven generations for all w.
3. **Hard post-zygotic gate:** w=0 eliminates all adult joint genotypes, even though connected configurations generate joint seeds. It need not preserve occupancy or total population size; no positive w policy should be called lethal.
4. **Overflow-specific positive selection:** at fixed pool size ≤K containing viable joint seeds, positive finite competition weights do not change which seeds are recruited. Their effect requires competition for fewer slots; w=0 remains a separate hard exclusion rule.
5. **Unique seeds and Mendelian copying:** never count the same produced seed as occupying two sites. Total before-recruitment pools equal seeds produced and post-recruit adult total cannot exceed that number. All A/B copies come from legal biallelic Mendelian gametes.
6. **Statistical independence:** one complete six-patch history, not the six patches or seven generations, is the simulation sampling unit. Reusing random keys creates paired *within-program* contrasts but is not biological replication. Report paired intervals and H7 negative/inconclusive contrasts without choosing a new horizon to force a discovery.
7. **Exact genotype support vs fitness:** conditional recruitment of a specific seed in a particular fixed pool is not unconditional multi-generation ecological establishment. A perfect selective advantage never reopens a mutation-free route that has been lost globally.
8. **Model-source segregation:** v8 imports v7's shared network/probability source, but it is **not the full izu-core individual Model 3** and does not reproduce EGWE's original eco-genetic fragmentation model or ANTO's failure-cost branching process. Real Cirsium phenotype, ancestry and seed observations remain exclusively in aza3.
9. **Realism warning:** the fixed q does not represent observed pollinator visitation, selection affects only competitive juvenile settlement, no donor fecundity or environmental feedback responds to genotype, and no independent taxa/lineages are analyzed. Only an independent natural-genetic and demographic panel can support a broad ecological theory.

## Commands and verification

Python 3.11+, standard library only:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_selection_filter.py --replicates 2000 --horizon 7 --seed 1182026 --output /tmp/sdea_v8.json
python scripts/verify_fixture.py results/v8_selection_filter_mc.json /tmp/sdea_v8.json
```

The benchmark uses only separate, predeclared selection × network cells and is not fit to earlier outcomes. Frozen results are in [`results/v8_selection_filter_mc.json`](../results/v8_selection_filter_mc.json), source in [`src/sdea/selection_filter.py`](../src/sdea/selection_filter.py), tests in [`tests/test_selection_filter.py`](../tests/test_selection_filter.py). The prior v7 and v7b fixtures remain as recorded.
