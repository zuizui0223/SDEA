# SDEA v6 — A shared finite seed budget changes joint evolutionary access across three patches

**Status: new post-v5 synthetic mechanism and exploratory 2026-10-26 Monte Carlo.** This protocol is a transparent definition of the model and its falsifiers, **not** a prospective preregistration, independent ecological experiment, source-repo replication, molecular floral-regain measurement, or fitted landscape model. No previous v1–v5 fixture, registered outcome, or biological evidence rank is changed. Original thistle observations remain in `aza3`.

## Primary general-theory question

Does matching **each** destination's probability of receiving a candidate trait imply that a landscape has the same probability of **simultaneous** access, when all destinations draw on the **same finite reproductive source**? Can seed export simultaneously increase colonization and reduce the donor's retention of its original genetic route? Does the pattern persist when seed output per donor is relaxed, pollen can move paternal alleles, and genotype-dependent recruitment is switched on?

Three separate claims are being challenged:

- **ANTO interface:** A/B route-associated alleles can be assembled within one diploid, but are only abstract _capacities_; neither their co-occurrence nor their visible expression establishes ancestral re-expression or viable descendants.
- **EGWE interface:** A patch-wise marginal arrival or trait-presence summary can fail to be a sufficient statistic for a **joint** future outcome. Here the missing source state is a finite **shared seed budget and allocation law**, not an asserted irreducible historical memory. Conditioning on complete source state and future process trivially closes a Markov kernel.
- **izu-core interface:** Viable reproduction, seed allocation, maternal/paternal separation, finite Mendelian genotype transmission, retained local alleles, current occupation and future reachability need not move together. These synthetic processes **are not the frozen Model 3 simulator**.

## Explicit model and event order

Three patches (`i=0,1,2`), each occupied by `0..K` diploid adults (`K=2` in reported scenarios). Each individual has two unlinked biallelic loci represented by `(k_A,k_B)`, each copy number in `{0,1,2}`. No mutation, external unsampled reservoir, genotype age or true ancestry is modeled. q is effective interaction state (`.25`/`.75`), not richness, floral color, wind, true pollen delivery or measured island distance. At time 0 all visible phenotypes are stipulated unobserved.

Each generation is synchronous and consists of:

1. **Seed production attempts:** each *occupied source patch* has `B` attempts, **total**, not B for each receiver and not B per adult. Draw a maternal diploid uniformly from surviving source adults. A potential self has attempt probability `s=a(1-q)` and succeeds with probability `1-δ`. Otherwise attempt outcross, which succeeds with probability `q` only if its paternal route points to a patch containing a permissible father (for a local father, **exclude the mother**; no father causes failure). Outcross paternal patch choice is a source-indexed stochastic matrix `P` with rows summing to 1.
2. **Genetics:** Each successful event yields one viable diploid seed by two independently assorting Mendelian gametes from the mother and chosen father. The maternal **genetic source** remains the patch where that seed was produced, even if the paternal gene arrived via pollen from another patch.
3. **Seed allocation:** Send every viable seed **exactly once** to one of the 3 patches using a maternal-source-indexed row of seed dispersal matrix `D`. Source-retained seeds and exported seeds compete for the **same B source attempts**. There is no exogenous seed subsidy, same seed delivered to multiple destinations or separate unlimited donor export pool. No fixed proportion of population size is guaranteed to reproduce.
4. **Density-dependent settlement:** If >K viable seeds arrive at a destination, sample **without replacement** a size-K subset. The selectable weight of genotype `g` in destination q is `exp[beta_A*(q-0.5)*k_A + beta_B*(0.5-q)*k_B]`. At <=K arrivals, all recruit regardless of fitness weights. Neither fitness coefficient has been measured in any natural thistle and baseline beta values are 0.
5. **Extinction shocks:** apply an independent post-recruitment patch catastrophe. If it occurs all recruited adults disappear; there is no seed bank in that patch. A currently empty patch has **no maternal reproduction**; paternal pollen alone cannot recolonize it.
6. **Environment and expression:** draw next low/high q (baseline `P(high|high)=.8`, `P(high|low)=.12`), optionally modified by **previous** adult A/B visible expression in that patch, clamped to `[0,1]`. After recruitment and q update, an adult with ≥1 A copy can start A expression with stipulated probability `.75q`; an adult with ≥1 B copy can start B expression with probability `.75(1-q)`. A within-individual simultaneous activation is resolved with a 50:50 tie. Visible expression is **not** a genotype mutation, genetic fixation, maternal fitness measurement or proof of a real ancestral pigment.

`run_trajectory` tracks first colonization of a patch initially vacant, recolonization **after later extinction**, cumulative first visibility, joint A/B genetic capacity in the same plant, per-generation seed attempts, successful seeds, exports and admitted recruits. Counts count **whole three-patch histories**, not individual seeds or patch-years as independent ecological replicates. `ensemble` uses independent fixed pseudorandom landscape seeds within each synthetic scenario; scenario-to-scenario uncertainty bands use independent-sample binomial Monte Carlo standard errors and are **not inferential CIs for real systems**.

## Exact shared-source theorem (single occupied source, two empty recipients)

Fix the donor reproductive state within the first generation and assume each of its `B` attempts independently yields a viable seed with probability `S`. Conditional on a viable seed, maternal allocation chooses destination 1 with probability `m_1` and destination 2 with probability `m_2`, with `m_1+m_2≤1`. Let `u_j=S*m_j` and let `c_j` be independent post-settlement catastrophe probabilities. Assume destination capacity at least 1 and no other sources. Then:

\[
\begin{aligned}
 P(C_j)&=[1-(1-u_j)^B](1-c_j),\\
 P(C_1\cap C_2)&=[1-(1-u_1)^B-(1-u_2)^B+(1-u_1-u_2)^B](1-c_1)(1-c_2).
\end{aligned}
\]

Consequently, in the no-catastrophe cell, `P(C1∩C2)≤P(C1)P(C2)` and the inequality is strict for `B≥1`, positive probabilities at both sites and nondegenerate conditions (except saturation/degenerate cases). **This is negative dependence induced by the mutually exclusive fate of each seed**, not a stochastic claim about spatial interactions in nature. A second independent budget per destination can have the same `P(C_j)` but would predict the product for their joint occurrence.

Exact donor retention (assuming no incoming source and donor capacity≥1) is:

\[
P(\text{source occupied after step})=[1-(1-S m_0)^B](1-c_0),
\]

where `m_0` is its retained-seed allocation probability. Increasing export allocation at **fixed B** reduces this retention chance while potentially raising receiver colonization. This observation is a model-specific allocation trade-off, not adaptive trait loss, inbreeding depression or an empirical dispersal cost.

**Analytical witness:** One source has two A-only adults, q=.75, assurance=.45 and depression=.20. Two recipient patches are empty. Hence `S=.755625`, `B=3`, allocations `(m0,m1,m2)=(.2,.4,.4)`, no catastrophe:

| Exact event | Probability |
|---|---:|
| Each recipient colonized (marginal) | 0.660296880140625 |
| Both colonized from finite shared budget | 0.382457969156250 |
| Both if recipient events wrongly treated independent | 0.435991969923443 |
| Joint minus independent | −0.053534000767193 |
| Source retains ≥1 locally allocated viable seed | 0.388310211580078 |

At B=1, colonizing **both** distinct receivers in the first step is impossible. At B=2, q=1 and exactly two viable seeds allocated 50:50, the joint event is `.5`, not `.75²=.5625`. The difference becomes smaller in a fixed cell when more seed attempts saturate both receivers (e.g., B=8 at the same per-attempt rates: joint `.88823424` vs independent `.89079203`), but saturation alone is not a universal convergence law at fixed expected seed number.

## Transparent exploratory eight-generation numerical experiment

**Eight new alternative kernels** (not a preregistered confirmatory comparison), each with 3,000 seeded whole-landscape trajectories, H=8. Initial source test: patch 0 has two A/A-route (A-only) diploids, patches 1 and 2 empty. Vary B and the 3×3 seed routing matrix:

- `single_local_B3`: all seeds stay at maternal patch, B=3.
- `single_moderate_B3`: maternal retention .65; each other patch .175, B=3.
- `single_strong_B3`: maternal retention .20; each other patch .40, B=3.
- `single_strong_B8`: the same strong export but B=8. Other ecological and reproductive parameters stay unchanged.

Initial mixed-route test: patch 0 holds two A-only adults, patch 1 is empty, patch 2 holds two B-only adults. B=3 and no genotype-selective recruitment at baseline. Compare:

- `distinct_seed_chain_B3`: nearest-neighbor seed routing, pollen retained locally.
- `distinct_seed_pollen_chain_B3`: same seed routing, pollen also restricted to chain neighbors.
- `distinct_seed_pollen_connected_B3`: all-to-all seed and pollen with .65 maternal-patch retention.
- `distinct_seed_pollen_connected_B3_selection`: identical connected routing but hypothetical `beta_A=beta_B=1.4` affecting **only over-capacity recruitment**.

No natural seed budget, pollinator range, local extinction rate or genotype fitness is inferred. Any scenario comparison estimates a **conditional model difference** only; ecological sampling hierarchy is the whole 3-patch replicate, not the number of flowers/seeds/time steps.

### Model outputs and adverse findings

The first-step joint dependence is an **exact finite-budget result**. The H=8 estimates below are finite Monte Carlo proportions (3,000 whole landscapes per cell) with source-fixed seed schedules; intervals in the machine JSON reflect Monte Carlo uncertainty only:

| Scenario | Source patch0 occupied at H8 | Initially empty patch1 colonized at least once by H8 |
|---|---:|---:|
| local seed retention, B=3 | 20.03% | 0% |
| moderate export, B=3 | 5.30% | 67.17% |
| strong export, B=3 | 4.53% | 84.17% |
| strong export, B=8 | 97.70% | 100.00% |

The unadjusted synthetic strong-versus-local source-occupancy difference is **−15.50 percentage points** (independent-ensemble 95% Monte Carlo interval −17.11 to −13.89 pp), while initial empty-patch colonization rises **+84.17 points**. Increasing **both** seed attempts and immigrant/local seed supply from B3 to B8 under strong routing changes source occupancy by **+93.17 points**; this does **not** isolate a dispersal effect from fecundity. This strong budget sensitivity is precisely the alternative mechanism under review.

For the mixed A/B scenario, fraction of histories that *ever* contain one dual-capable individual in patch1 by H8: **16.67%** (seed chain), **19.13%** (seed+pollen chain), **26.00%** (connected seed+pollen) and **28.13%** (connected plus selective recruitment). The selection contrast for **patch0 A support** is **−0.70 pp** with independent-ensemble 95% Monte Carlo interval about **−1.78 to +0.38 pp**, i.e. **inconclusive in this finite synthetic sample**. Do not claim selection generally erodes or saves A based on this trial. Adding pollen under a chain also changes viable reproduction and source extinction; no single clean controlled direct effect of pollen on joint retention is identified by the descriptive comparison.

## Negative and falsification controls

1. **Seed conservation:** per source and generation, `successful_seeds≤attempts≤B`; every successful seed is allocated to exactly **one** destination, such that `Σproduced=Σdestination_arrivals`; exported+retained equals produced; settled recruits cannot exceed arrivals or local K. No seed cloning among recipients.
2. **Maternal absence:** empty source creates no seed; pollen alone never creates occupancy in empty recipient.
3. **Paternal specificity:** foreign pollen changes offspring genetics only for viable outcrosses with resident mother. It is not the same mechanism as seed immigration.
4. **Closure:** an allele absent from all patches never appears because mutation and external immigration are disabled. A+B within an individual is not visible trait restoration.
5. **Shared-budget analytical calibration:** with B=1, two initially empty receivers cannot both colonize in a generation; with B=2, q=1 and 50:50 allocation, exact joint = 0.5. Exact colony+source retention formulas agree with independently seeded finite draws to predeclared tolerance.
6. **Competition:** genotype selection acts only when viable arrivals exceed capacity; zero effect when all arrivals fit. Conditional positive selection is not a measured genetic coefficient.
7. **Extinction and recolonization:** local shock occurs *after* recruitment, a previous loss flag is distinct from a patch initially empty at time zero, and re-colonization requires incoming seed in the restoration generation.
8. **Independent source reconstruction:** all new exact and Monte Carlo outputs are regenerated from code and checked against one separately frozen JSON fixture. Existing v1–v5 fixtures are never silently rewritten.
9. **No unqualified generality:** thresholds, signs, long-horizon occupation and convergence depend on capacity, B, dispersal topology, selection, q, donor viability and future stochastic forces. No retrospective relabeling of an exploratory cell as confirmatory evidence.

## Reproduction and source boundary

Standard-library Python (>=3.11):

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_finite_budget.py --replicates 3000 --output /tmp/sdea_v6.json
python scripts/verify_fixture.py results/v6_finite_budget_mc.json /tmp/sdea_v6.json
```

Full source: [`src/sdea/finite_budget.py`](../src/sdea/finite_budget.py); new negative tests: [`tests/test_finite_budget.py`](../tests/test_finite_budget.py). Immutable Monte Carlo/analytic record: [`results/v6_finite_budget_mc.json`](../results/v6_finite_budget_mc.json). The experiment is independently written **within SDEA** and does not re-run ANTO, EGWE or izu-core code. The use of an abstract `A` route does not establish phylogenetic ancestry, trait re-expression rate, individual-level seed viability in a real species, ecological adaptation or cross-clade generality.

**Next biological gate:** measure source-constrained actual propagule production and realized dispersal, pollen parentage, same-individual functional route genotypes/expression, independent ecological q, and surviving descendants with whole-population held-out prediction. The real Cirsium molecular/field programme remains under [`aza3`](https://github.com/zuizui0223/aza3); no new sampling, consent or field validation was performed for SDEA v6.
