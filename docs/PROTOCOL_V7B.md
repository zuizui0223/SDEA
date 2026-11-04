# SDEA v7b — Matched reproductive budget, different directed dispersal topology

**Date:** 2026-10-27. **Role:** complementary *three-patch* directed-reachability stress test to the existing *six-patch* SDEA v7 degree/flow-alignment analysis. SDEA v7 is not replaced or reinterpreted. **Evidence level:** a new, explicitly synthetic model-class restriction and exploratory finite-horizon sensitivity. Results were generated from prescribed contrasts and inspected in this session; this is **not** a retrospective biological preregistration or a natural replication of ANTO, EGWE, izu-core or SDEA v6. Previous v1–v7 code, results and inferential boundaries remain unchanged.

## Central question

> When the *initial* reproductive budget, allelic state, mating system, local environmental quality, selfing and fraction of seeds exported are identical, does changing **which patch is connected to which** alter future access to a lost local site or the assembly of A/B candidate genetic routes?

This sharply differs from the v6 B=3 versus B=8 comparison, which changed reproductive supply. The **opportunity** per occupied maternal patch stays fixed at three seeds per generation in every v7 scenario. Realized total seed production across the landscape can subsequently diverge because occupancy evolves; that divergence is an *outcome*, not a newly changed input.

## Restricted generative model

- **Three patches** with capacity K=2 diploid adults apiece; two freely assorting, unlinked biallelic loci. `A=(2,0)` and `B=(0,2)` are *hypothetical route-associated abilities*, not validated ancestral pigment modules or floral colors.
- **Shared seed production per maternal source:** B=3 attempts per occupied patch per generation. Each attempt produces at most one viable zygote and assigns it to exactly one patch. A failed attempt produces no seed. The donor cannot clone one seed into multiple destinations.
- **Environment:** effective pollination value q=0.75 in **all** patches in **all** generations. No temporal forcing change and no phenotype–environment feedback in the primary network comparisons. The q value is not observed pollinator richness or visitation frequency.
- **Reproduction:** assurance a=0.45 and selfing attempt share s=a(1−q)=0.1125. Selfed viable zygote weight s(1−δ) with stipulated δ=0.2; a distinct local paternal parent is required for outcross and an outcrossing seed succeeds with probability q. The primary pollen-location matrix is identity, so immigrant pollen is *absent* from the network-only experiment. Parent gametes segregate Mendelially at each locus, with no mutation, selection or demographic catastrophe.
- **Competition:** At most K=2 surviving recruits per receiving patch, with uniform lottery without replacement when more than K seeds arrive. No seed bank, perennial adults or additional mortality mechanisms. All adults are replaced by recruited offspring each generation. Source and recipient updates use the previous adult state synchronously.
- **Initial state:** primary family has `((A,A), (), ())`, one A-only source and two genuinely vacant patches. A second, **separately analysed** family uses `((A,A), (), (B,B))` to study joint A+B capacity. A third symmetry control uses `((A,A), (), (A,A))`.
- **Time horizon:** one/two-generation exact operator checks plus H=7 exploratory independent landscape Monte Carlo, n=4,000 independently seeded *whole landscapes per treatment*, not patches/propagules counted as replicate systems. Eight scenarios = 32,000 simulated landscapes, not 32,000 independent wild islands.

### Strictly matched seed-routing treatment matrices

Rows identify maternal **source**, columns the **destination**. All tested matrices have diagonal 0.6, row-sum 1 and exactly **0.4 export probability per seed**. Seed production B=3 and all other parameters are identical.

```text
cycle (0 -> 1 -> 2 -> 0)     return (0 <-> 1, 2 -> 0)
[[.6, .4, 0],               [[.6, .4, 0],
 [0, .6, .4],                [.4, .6, 0],
 [.4, 0, .6]]               [.4, 0, .6]]

balanced                        reversed cycle (symmetry control)
[[.6, .2, .2],              [[.6, 0, .4],
 [.2, .6, .2],               [.4, .6, 0],
 [.2, .2, .6]]              [0, .4, .6]]
```

**Important caveat:** row sums, source retention, total export, q and B are fixed. **Column sums/incoming connectivity and directed paths are not always matched**, because precisely these are part of the network structure under intervention. Both cycle and reverse have column sums 1; cycle versus return differ in column sums. Do not call the manipulation matched *immigration* or interpret effects as independent of destination-side incoming flow.

## Exact two-generation result: same full t1 state law, different t2 future

In the primary one-A-source setting, only patch0 has adults at t0, and its routing row is **identical** in cycle and return. The complete joint genotype/population state distribution at the first recruited census is therefore identical under the two networks. The first-generation offspring from A/A mothers and fathers are all A/A. In particular, donor's per-seed viability is

\[
S_2 = s(1-\delta)+(1-s)q = 0.755625,
\qquad s=a(1-q)=0.1125.
\]

The per-attempt probability of successful seed arrival to patch1 is `p_01 = 0.4*S_2 = 0.30225`. Thus

\[
P(C_1(t=1))=1-(1-p_{01})^B=0.660296880140625
\]

under **both** treatments. Because there is no 0→2 connection and patch2 begins empty, it cannot be occupied after t1 under either treatment.

At the second generation, if patch1 holds n recruited diploids, its viable attempt probability is

\[
S_n = \begin{cases}
0 & n=0,\\
s(1-\delta)=0.09&n=1,\\
S_2=0.755625&n=2.
\end{cases}
\]

Let X~Binomial(B,p_01) be the number of seed arrivals to patch1 in the first generation, and N_1=min(X,2) the admitted mothers. For any routing matrix with no 0→2 edge and local pollen, the t2 occupancy probability of initially empty patch2 is

\[
P(C_2(t=2))=\sum_{x=0}^B {B\choose x}p_{01}^x(1-p_{01})^{B-x}
\left[1-(1-m_{12}S_{\min(x,2)})^B\right].
\]

The **cycle** has m_12=.4 and therefore **P=0.1904815097464421**. The **return** has m_12=0, hence **P=0 exactly**. An independent explicit `3^B` categorical allocation enumeration reproduces both values with machine precision. With no external source/mutation and initially empty patches, a seed-directed path from an initially occupied source is necessary for colonization. A pollen-only path is not sufficient without a receiving mother. A positive path is only a *possibility* when seed production/viability remain positive, not a guarantee of long-term persistence.

The exact first-step identity is stronger than simply matching two summary statistics: the **entire full current genotype-state distribution at t1 coincides**, while the *future movement laws* differ. Under a known fixed transition operator this is not a violation of the Markov property; a forecast must condition on that operator.

## Exact observer-design experiment: knowing future routing matters

Suppose the network intervention label is equally likely to be cycle or return and is hidden from the observer at t1. The observer knows the **entire** t1 genotype census, and therefore n1, but not the future routing operator. The binary outcome is whether patch2 is occupied at t2.

For each n1 state, let `p_cycle(n)` and `p_return(n)` be the exact colonization probabilities shown above. The Bayes-optimal expected Brier risk with the operator hidden is `E_n[p_bar(n)(1−p_bar(n))]`, where p_bar is their 50:50 mixture. With the operator revealed it is `E_n[(p_cycle(1−p_cycle)+p_return(1−p_return))/2]`. The difference is the nonnegative identity

\[
\Delta B=\frac14\mathbb E_n\left[(p_{cycle}(n)-p_{return}(n))^2\right].
\]

| Observation given at t1 | Exact expected Brier risk |
|---|---:|
| Complete current genotypes; network law hidden | **0.0701901781122** |
| Complete current genotypes plus known future routing law | **0.0451396013512** |
| Gain from observing law | **0.0250505767610** |

The decrease is ~35.7% of the hidden-law Brier risk. If both candidate operators are the same, the gain is exactly zero. This is a **model-specific experimental-design result about an omitted forcing/operator**, not evidence of irreducible past-history dependence.

## H=7 exploratory finite-history outcomes

Monte Carlo: **4,000 independently seeded whole landscapes per cell** and new seed streams across cells. These values have *simulation* uncertainty only; no natural biological CI. Separate founder families are not pooled. Event definitions are fixed before writing the result table.

| Founder family, seed network | t1 patch1 colonized | t2 patch2 occupied | patch2 ever occupied by H7 | All three occupied at H7 |
|---|---:|---:|---:|---:|
| Single A source, cycle | 65.73% | **19.03%** | **45.15%** | 0.63% |
| Single A source, return | 67.20% | **0%** | **0%** | 0% |
| Single A source, balanced | 38.38% | 31.05% | 67.43% | 0.85% |

The single-source t2 cycle–return Monte Carlo contrast is +19.025 percentage points with approximate independent-simulation 95% interval **+17.81 to +20.24 points**. This is independently consistent with the **exact** 19.048% cycle probability. Note that the balanced network sends only 20% (not 40%) of source0's seeds directly to patch1: the balanced cell does not match the focal patch's first-step *immigration rate* and should not be compared as if it did.

For **two distinct preexisting A-only and B-only sources**, the event "an A+B-capable individual ever occurs in middle patch" was **43.20%** (cycle), **46.90%** (return) and **44.70%** (balanced). Return–cycle difference **+3.70 percentage points** (simulation 95% interval approximately +1.52 to +5.88). Balanced–cycle difference **+1.50 points** with simulation interval approximately **−0.67 to +3.67**: **inconclusive** at this computational precision. All genotype labels are candidate genetic possibilities, not expressed composite floral traits or verified ancestral regain.

**Symmetric adverse control:** Two A-only sources at patches0 and2, with cycle vs its patch0↔patch2 permutation (reverse), are distributionally symmetric for patch1-colonization outcomes. The estimated patch1-ever-colonized difference is −0.25 percentage points (simulation 95% interval approximately −1.17 to +0.67); **the null is not rejected, but not established as a biological equivalence**. Exact graph permutation provides the structural symmetry, rather than relying only on the interval.

### Accounting and denominator warnings

All treatment cells have fixed **per occupied-source** seed-attempt budget B=3. As colonization/extinction changes occupied source number, **total realized seed attempts and viable seed production vary as downstream outcomes**. E.g. one-source H7 cycle averaged 7.43 viable seeds per full history vs 7.45 in return; a balanced network averaged 6.43. These realized totals are not experimentally held equal and must not be used to claim a pure fixed-total-realized-propagule effect.

The model replaces all adults every generation, lacks observed seed banks, long-lived plants, resource competition between different maternal plants, genotype-dependent seed success, environmental feedback, distinct actual floral route expression and natural pollinator histories. Pollination q is fixed by assumption; the nominal pollen matrix is held local to isolate seed topology. A functional network is **not** identical to a geographical distance or island size matrix.

## Adversarial and identifiability gates

1. The cycle and return must have identical first-generation *full joint genetic/population-state distributions*, since only routing from initially empty patch1 changes. A same-seed one-step counterfactual checks this equality for a set of seeds.
2. Exact binomial path calculation and an independent `3^B` enumeration must agree. A future simulation ensemble at H2 must agree with the analytic prediction within a prospective tolerance; a structurally disconnected patch must never appear, including with immigrant *pollen* but no incoming *seeds*.
3. Seed-production attempts per occupied source are bounded; a successful seed cannot be sent to two destinations. All adult counts stay ≤K. Neither A nor B appears without parental input or mutation.
4. Network diagonal and row-wise seed export fraction remain fixed, while incoming sums and directed paths are allowed to change as explicit network treatment coordinates.
5. If the two future dispersal operators coincide, adding the operator label must yield zero exact Brier improvement. q, mating, seed budget, genotype starting state and pollen location must remain matched in the observer calculation.
6. Symmetry control has identical theoretical patch1 fate under a simultaneous network/founder-label permutation. Under simulation, differences should be reported with intervals; intervals spanning zero are **not proof of equivalence**.
7. Do not interpret an underpowered genotype-route contrast, the Brier score or graph reachability as proof of actual floral trait reacquisition, stochastic drift fraction or long-run survival.
8. For natural validation, independent material would require a donor–recipient seed and pollen flow census, *viable* seed budget per source, multilocus parentage/ancestry, real same-individual functional module assays, measured post-colonization genotype retention, offspring performance and whole-landscape holdout.

## Reproducibility

Run with Python 3.11+ and the standard library; no independent external biological data are loaded:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_network_reachability.py --replicates 4000 --output /tmp/sdea_v7b.json
python scripts/verify_fixture.py results/v7b_network_reachability_mc.json /tmp/sdea_v7b.json
```

Frozen numeric readout: [`results/v7b_network_reachability_mc.json`](../results/v7b_network_reachability_mc.json); new standalone source: [`src/sdea/network_reachability.py`](../src/sdea/network_reachability.py); adverse tests: [`tests/test_network_reachability.py`](../tests/test_network_reachability.py). The standalone source is a purposefully constrained **independently implemented** check of v6's seed-budget/Mendelian logic, not the original source EGC or izu-core process itself. SDEA v7 and all older results remain reproducible separately.
