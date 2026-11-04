# SDEA v7 — Matched seed budgets, network placement, and evolutionary accessibility

**Date: 2026-10-27. Status: fully constructed, source-separated synthetic model and 7-generation finite-horizon experiment.** v1–v6 are preserved. No natural pollen/seed movement, candidate pigment-pathway function, fitness, ancestral regain, or long-term field persistence has been inferred. Numerical costs, reproductive dynamics and genotype-to-flower expression are assumptions rather than estimates.

## Scientific question and primary matched comparison

**Does the placement of reproductive links relative to initial genotype spatial structure change the formation and eventual retention of compound genetic capacity when total initial seed-production opportunity and source/destination movement marginals are matched?**

The principal **connected-ring** treatments are *isomorphic as unlabelled graphs*; they differ in **genotype-to-edge alignment**, not generic six-node degree, number of edges, or connectedness. A closed-components treatment adds a structural impossibility control. Do not frame the clustered-versus-alternating test as a proof that one unlabelled graph topology is intrinsically better.

Initial six patches are fully occupied with two A-only homozygous diploids each in patches 0–2 and two B-only homozygous diploids each in patches 3–5. There are two abstract, unlinked, biallelic route-associated loci with diploid genotypes `(A copies, B copies)`; `A=(2,0)`, `B=(0,2)`. Neither locus has verified molecular ancestral status or is mapped to flower color. Every patch has the same dimensionless effective pollination index **q=.75** for all seven generations; all initial phenotypes are *unobserved/absent*. There is no mutation or external allele introduction.

Each occupied patch has K=2 adult slots and **B=3 seed-production attempts per generation**; the initial total reproductive attempt budget is exactly **18** in all conditions. Selfing is attempted with `s=a(1-q)` with assurance `a=.45` and seed viability `1-δ=.80`. Outcross viability is q conditional on an eligible father. Two loci segregate independently through Mendelian gametes. Parents can be chosen for separate attempts, but a successful seed is an individual with **one maternal source and one destination**; a seed cannot be cloned to multiple destinations. A neutral without-replacement adult recruitment lottery retains at most K seeds per patch. No genotype-dependent reproductive fitness or immigrant advantage is built into this model. Local population absence blocks pollen fertilization (no mother), but seeds from occupied patches can colonize it. There is no additional catastrophe or seed bank in the baseline.

**Important qualification:** Future realized whole-landscape seed attempt totals are **not held constant**, because a patch produces seed only if it has an adult. Different network placements can change occupancy and hence realized future budgets. Fixing that realized output post hoc would condition on a potential mediator. The controlled effect is an *ex ante network-policy contrast* with exactly matched reproductive and genomic starting conditions and a common per-occupied-patch budget, not a global constant-output intervention through all seven generations.

## Strictly matched networks

Undirected 2-regular graphs over the same six labelled patches:

| Label | Node cycle(s) | A–B interpatch directed links |
|---|---|---:|
| `components` | `(0,1,2)` and `(3,4,5)` | 0 |
| `clustered` | `(0,1,2,3,4,5)` | 4 |
| `alternating` | `(0,3,1,4,2,5)` | 12 |

When active, each source assigns `.25` of either seed destination or pollen father-location weight to itself and `.375` to **each of its two neighbors**. Every routing matrix is symmetric and doubly stochastic: all six outgoing and incoming marginal weights equal 1. All have **six undirected edges**, degree 2 at every node, identical local-retention share and identical total input/output. Only the **specific contact matrix** differs. An `identity` matrix (all local) is an explicit **no-interpatch-flow control**, not degree matched to active networks.

Operators are crossed in three initial blocks: seed-only (seed graph + local fathers), pollen-only (local seed destination + pollen graph), both (the same graph in both). Each block includes `components`, `clustered`, and `alternating`. A no-flow baseline is added. **Two additional crossed cells** (seed clustered/pollen alternating and seed alternating/pollen clustered) and the 2×2 interaction analysis were appended **after observing preliminary 7-generation endpoint saturation**; report them as exploratory, not as independent prospective confirmation. In total there are **12 model cells**.

## Analytical first-generation genetic assembly theorem (model conditional)

Under the frozen homozygous A-only versus B-only founders, let `c_i` be the probability that a viable outcross paternal source for mother-patch `i` comes from a patch carrying the **opposite** route. For a routing matrix `M`, `c_i=Σ_(j of opposite route) M_ij`. Because every cross-type A×B zygote has genotype `(1,1)`, and each initial seed attempt is independent conditional on the fixed founders,

\[
  h_i=(1-a(1-q))q\,c_i,\qquad
  P(\text{at least one viable A+B seed produced at generation 1})
       =1-\prod_{i=0}^{5}(1-h_i)^B.
\]

`q=.75`, `a=.45` and `B=3` imply `(1-a(1-q))q=.665625`. No movement across routes (`components` or local fathers) yields **zero** immediately produced dual-capacity zygotes. The clustered cycle yields expected hybrid seed count **2.9953125** and an at-least-one probability **0.9681251026678412**. The alternating cycle yields expected count **8.9859375** and probability **0.999996076577506**. These are **seed production** probabilities, *not recruitment, sustained expression or ecological establishment*. They are independent of seed destination routing for this first formation stage. The exact formula was checked separately against seeded independent Monte Carlo, Mendelian crossing and a seed-route-swap adverse control.

### Binary observation saturation and state-information audit (exploratory)

A second **post hoc analytical audit**, performed after seeing the H7 saturation, varies the **same** fixed per-source attempt budget `B` across separate hypothetical experiments. *Within each B*, the two connected ring policies still have exactly matched seed opportunities and genotype states; comparisons **between B values** do change resource availability and must not be interpreted as network-only interventions.

| Attempts per source B | P(hybrid seed at first generation), clustered | Alternating | Brier risk improvement from knowing the network label |
|---:|---:|---:|---:|
| 1 | 0.68293406 | 0.98422794 | 0.02269450 |
| 2 | 0.89946919 | 0.99975124 | 0.00251412 |
| 3 | 0.96812510 | 0.99999608 | 0.00025394 |
| 8 | 0.99989786 | ≈1 | ≈0.000000003 |

With an equal prior over clustered/alternating labels, let `p_C` and `p_A` be their exact binary hybrid-production probabilities. Bayes-optimal binary Brier risk under a coarse label-ignorant predictor is `B_coarse = p̄(1−p̄)`, whereas including the network gives `B_network = [p_C(1−p_C)+p_A(1−p_A)]/2`; the exact gain is

\[
   B_{coarse}-B_{network}=\frac{1}{4} (p_A-p_C)^2.
\]

This is a **structural demonstration of information loss when the presence/absence endpoint saturates**: at fixed B=3 the alternating network yields three times as many *expected* hybrid-seed production events as the clustered network (8.98594 vs 2.99531), while an any-seed indicator differs by only about 3.19 percentage points. There is no empirical prediction training/test dataset, no evolutionary clock, and no proof that photographs or observed flower-color states follow this event law. The difference in Bayes risks measures the constructed forecasting task only.

The `components` result is a **structural zero for the entire future process** if both pollen and seed matrices are confined to separate A-only and B-only triangles, with zero mutation/external immigration: no individual can ever contain both genetic route loci. The two disconnected components can each retain a different genetic route, so **both routes present somewhere across the six patches** is not equivalent to an **individual carrying both**.

## Multi-generation experiment and endpoints

Each of the **12 cells** is run for **H=7** generations, with **2,500** whole-landscape Monte Carlo replicates each: **30,000 six-patch histories**. The same seed key (`981027 + 1_000_003 × replicate index`) is used across all cells to construct descriptive paired contrasts. The unit for Monte Carlo precision is the **complete six-patch history**, not a seed, site, individual, or generation. As altered policies consume the pseudo-random stream differently, common seed keys do not guarantee optimal coupling; the empirical variance of within-key differences is used in the reported paired normal-approximation MC intervals. Full code and readout use `random.Random.random()` for portable deterministic draws across Python minor versions.

Report distinct events: (i) at least one viable dual-capacity zygote **produced** at generation one; (ii) a dual-capacity diploid **recruited** in that generation; (iii) dual-capacity diploid **ever recruited** by H; (iv) any such diploid still **present at H**; (v) both loci present somewhere across the landscape at H; (vi) any patch occupied at H; (vii) cumulative seed attempt/output/export sums. Also record same-patch genotype co-occurrence and phenotype-activation *only as descriptive synthetic filters*; there is no selection imposed from these labels.

**Primary comparisons:** alternating minus clustered within matched seed-only, pollen-only, and both-active treatments; and clustered minus disconnected-components as structural controls. The latter changes whole-network connectedness. Positive or negative direction beyond the one-generation exact expression is **not guaranteed**; no multiple-testing-corrected confirmatory long-horizon model-selection claim is made.

## Frozen executed results (portable PRNG, Python standard library)

| Treatment | Ever recruited A+B by H7 | A+B adult present H7 | Both routes somewhere H7 | At least one patch occupied H7 |
|---|---:|---:|---:|---:|
| Seed-only clustered | 0.9620 | 0.2688 | 0.2892 | 0.5684 |
| Seed-only alternating | 0.9692 | 0.2600 | 0.2780 | 0.5256 |
| Pollen-only clustered | 0.9988 | 0.8292 | 0.8728 | 0.9648 |
| Pollen-only alternating | 1.0000 | 0.8292 | 0.8596 | 0.9660 |
| Both clustered | 0.9992 | 0.6028 | 0.6296 | 0.8228 |
| Both alternating | 1.0000 | 0.6280 | 0.6556 | 0.8380 |
| Both disconnected components | **0** | **0** | 0.3556 | 0.8364 |

Numbers are conditional *synthetic* cell probabilities. See machine results for exact stored values and Wilson intervals; rounded table is navigational, not a replacement for the JSON. The rounded values have been checked against the frozen v7 result fixture.

The **H7 final A+B adult outcome** (alternating − clustered), with 95% *paired simulation-precision* intervals, is:

- Seed-only: **−0.88 percentage points** `[−3.26, +1.50]`.
- Pollen-only: **0.00 percentage points** `[−1.95, +1.95]`.
- Both operators: **+2.52 percentage points** `[−0.10, +5.14]`.

Thus direct first-generation joint-zygote production can differ sharply **without a resolved direction in seven-generation joint-genotype retention** between the two *connected* cycles. This is neither proof of biological equivalence nor evidence that gene flow invariably improves future genetic options. The long-horizon null/inconclusive contrasts are retained as main findings.

The added exploratory seed × pollen topology factorial defines for any chosen indicator `Y`:

`DID = E[Y | seedAlt,pollenAlt] − E[Y | seedAlt,pollenCluster] − E[Y | seedCluster,pollenAlt] + E[Y | seedCluster,pollenCluster]`.

At first viable hybrid-seed production the DID is **exactly zero** under the first-generation germline production operator, because seed destinations are chosen *after* producing zygotes. The H7 final A+B indicator DID is exploratory and its paired MC interval includes zero (see machine results). Do not interpret this as an identified genetic-drift component or pathway-specific biological fitness effect.

## Falsification controls and evidence restrictions

1. Each topology must match all row/column marginal movement weights, node degree, total edges, retention, initial genotypes, q and source B. If not, the matched comparison is invalid.
2. Every seed production event consumes at most one finite source attempt and yields at most one seed in one destination; seed-production accounting equals total pre-recruitment seed pools.
3. A pollen-only network must not colonize a patch lacking a maternal individual; a seed network may. No external input and no locus B initially present anywhere implies no spontaneous B or AB.
4. No cross-route edges and no external introduction imply no AB capacity can arise; adding inter-route edges only makes first-gen formation *possible*, not automatically established at H.
5. Changing seed routing while keeping pollen routing identical must not change first-generation **produced** dual-capacity seed probability (even though recruit placements can change). A disconnected A/B source system with no cross-route movement must remain genotype-mixed-negative for every seed.
6. The exact closed-form hybrid-production probability must agree with finite stochastic reproduction; check both extreme k=0 and nontrivial cases. E.g., `B=0` forces zero reproduction, and all-empty systems remain empty.
7. Do not label the binomial Monte Carlo precision intervals as natural-system statistical confidence intervals. Approximate individual tests are uncorrected exploratory comparisons. A null or sign-indeterminate H7 contrast remains reported without changing horizons or cell definitions to recover a positive result.
8. The present-state sufficiency hypothesis remains conditional: complete genotype arrangement, maternal and paternal routing matrices, resource budget, occupancy and future forcing are needed for the declared Markov transition. The work does **not** prove all such coordinates can be measured in wild plants.

## Reproduction

Python 3.11+, standard library only:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_network_topology.py --replicates 2500 --horizon 7 --output /tmp/v7.json
python scripts/verify_fixture.py results/v7_network_topology_mc.json /tmp/v7.json
```

The new protocol, source and fixture are versioned in SDEA only. ANTO retains route accessibility/re-expression distinctions, EGWE state-sufficiency/ecological operators, and izu-core its own independently defined Model 3. The present generator is a **new cross-repository robustness thought experiment**, not a port or source-level independent replication of the three. For external testing, independently measured realized pollen paternity, source-limited seed output, parent genotypes, environmental interaction and multigeneration descendant phenotypes (including the ancestry/re-expression alternatives) remain necessary. *aza3* retains all natural Cirsium specimen/experimental ownership.
