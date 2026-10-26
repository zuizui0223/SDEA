# SDEA v2 — Exact finite-horizon evolutionary accessibility

**2026-10-26; an explicitly new exploratory synthetic model.** The v1 preread contract and its numerical fixture remain unchanged. No new natural observation, measured mutation rate, phenotype reversal, descendant fitness, or independent reproduction of ANTO/EGWE/izu-core is claimed. v2 converts the earlier **one-step** thought experiment into a **multi-generation, fully enumerated finite-state chain**, without importing any parent-model parameter or falsely calling a 3-generation streak fixation.

## Testable questions

1. **Can first-visible expression be common while three consecutive *survived* expression updates and final occupied lineages are rare?** Report the three events on one unconditional denominator.
2. **Can identical present-day visible phenotype and pooled functional resources yield different repeated-expression chances, simply because the joint resources occur where q is high?** Swap mask-to-patch assignments while retaining patch q and mask marginals. Use equal q as the negative control.
3. **Does ancestry/history add forecast information when current joint modules are incompletely observed, and disappear after the full predictive state is observed?** Compare Bayes Brier risks; add uninformative-history and perfect-assay controls.
4. **Can assurance benefit occupancy yet diminish accessibility when a specified direct allele-transmission cost is allowed?** Compare high/low assurance, keeping initial resources, allele copies and forcing law identical. Zero cost is an adverse control.
5. **Are there separate conditional assurance-cost crossover boundaries for transient expression, surviving lineage expression and genetic support?** Evaluate distinct endpoint definitions *before* root-finding, and fail without a sign bracket.

## Finite, conditional state and transition

An occupied synthetic lineage at generation `t` has `(m, k, q, r, seen, ever3)`: `m=0..7` is the joint functional module bit mask; `k=0..6` is number of target allele copies in a constant `2N=6` slots; `q` is one of two effective-interaction values (0.2, 0.8); `r` is current survived-expression streak 0..3; `seen` records any earlier visible expression and `ever3` a previously completed three-generation run. Extinction is absorbing. **All founders have a visibly absent phenotype at time 0 even when `m=7`**; phenotype at subsequent censuses is generated anew.

For each occupied lineage and generation:

1. Expression `E_t~Bernoulli(1[m=7 and k>0] αq)` with `α=.75`.
2. Reproductive occupancy gate `S_t~Bernoulli(q+(1-q)a)` independent of the expression draw conditional on state. Expression counts as a first appearance even if survival subsequently fails in that generation; a streak counts only *successful survived* updates.
3. Conditional on survival, future allele-copy frequency follows `w=exp[β(q-0.5)-ca]`, `p_s=(k/2N)w/[1-(k/2N)+(k/2N)w]`, optional one-way mutation and external migrant replacement, then `k'~Binomial(2N,p_eff)`. **This is NOT Mendelian individual mating, diploid genotype tracking, or the izu-core Model 3 operator.**
4. Each of three modules is retained independently with probability `ρ` if present and reacquired with probability `μ_M` if missing. After a visible expression, optional **hypothetical shared-toolkit collapse** with probability `d` sets `m'=0`. This is not ANTO v15's cost-of-failed-expression mechanism, nor a measured gene-loss event.
5. Ecological q switches with baseline `P(high'|high)=.85`, `P(high'|low)=.15`, plus an optional nonnegative expression-feedback increment (clipped to one). `q` is not richness or visits.

Default: `N=3`, `a=.30`, `β=.50`, `c=1`, `ρ=1`, `μ_M=0`, `d=0`, no allele input, no expression feedback. The focal complete state begins `m=7`, `k=2`, `q=high`. Histories use the exact v1 even/odd parity mask priors, both currently silent; other coordinates are matched. An independent two-patch union is `1-(1-P1)(1-P2)`. There is **no between-patch dispersal** in this analysis; external immigration is a stipulated input, not a coupled landscape transition.

## Explicit non-equivalent endpoints

- `first_expression_ever`: visible once by horizon, including lineages that subsequently go extinct.
- `three_consecutive_survived_ever`: ever completes 3 successive visible AND survived updates, including subsequently extinct lineages.
- `alive_and_ever_three`: still occupied at horizon after a completed 3-update streak at any earlier point.
- `alive_and_current_three`: survives at horizon with 3 successful consecutive updates ending there.
- `alive_and_joint_toolkit_supported`: occupied at horizon with `m=7` and target allele `k>0`, **not** necessarily expressing.
- `alive`: occupied at horizon without any requirement to retain phenotype or target allele.

These are probability-of-events under a finite model. None implies population-level fixation, long-term demographic nonextinction, actual evolutionary restoration or a molecularly measured pathway.

## Completed exact calculations, H=8

| Initial / treatment | First visible ever | Three survived in a row ever | Alive and ever three | Alive at H |
|---|---:|---:|---:|---:|
| Complete toolkit, k=2, q=high | 0.818665 | 0.117417 | 0.041242 | 0.131318 |
| Even parity resource history | 0 | 0 | 0 | 0.131318 |
| Odd parity resource history | 0.204666 | 0.029354 | 0.010311 | 0.131318 |
| Complete toolkit, zero initial allele copies and no input | 0 | 0 | 0 | nonzero |
| Complete toolkit, post-expression collapse probability 1 | 0.818665 | 0 | 0 | nonzero |

Matched 2-patch totals (one complete, one empty toolkit; one high and one low q) generate 3-survived-ever probabilities **0.117417 aligned** vs **0.008750 misaligned**, difference **0.108667**. Equal-q control gives difference zero.

At `c=0`, assurance 0.1→0.9 raises H8 occupancy **0.077828→0.737129** and the probability of an alive lineage that previously reached 3 successive expressions **0.034426→0.238636**. At `c=3`, occupancy changes identically, but the latter becomes **0.025414→0.002717**. The counterfactual transmission cost is constructed, not inferred from observed selfing or inbreeding depression.

### Three distinct exploratory frontiers (same H=8)

For assurance `a_low=.1` versus `a_high=.9` and an initial complete toolkit, numerical bisection of the high-minus-low contrast gives:

| Endpoint | Conditional cost crossover c* |
|---|---:|
| Three consecutive survived expression updates at any time by H | 0.805469 |
| Alive at H with complete modules and at least one target allele | 0.714673 |
| Alive at H after ever completing a 3-survived-expression run | 1.590315 |

Roots are bracketed within `[0,8]`, with residuals below `3e-15` for the declared cells; we did **not** prove global uniqueness or a universal biological ordering. Because the event and the horizon differ, the cost threshold has no definition independent of an explicitly selected endpoint. Changes in module repair, mutation, population size or interaction feedback may alter or eliminate a crossover.

### EGWE-style prediction-sufficiency check

For the H=8 `alive_and_ever_three` event, exact Bayes-optimal Brier risks are: visible alone **0.00512868733**, visible + parity history **0.00510211058**, noisy resource assay alone **0.00501901033**, assay + history **0.00497653075**, true complete module inventory **0.00494265010**, full inventory + history **identical**. The law and other coordinates are held fixed, so residual history would only indicate omitted state/incorrect closure in any external application. Parity-histories with the same prior mask distribution yield no history advantage.

## Falsifiers / controls

- All transition rows sum to one, occupied + extinct mass sums to one at every generation (within stated numerical tolerance).
- Zero toolkit and no module reacquisition: no expression; reacquisition explicitly permits expression. Zero target allele copies and no mutation/migration: no expression; mutation/input can reopen.
- Post-expression damage=1, no repair: first return remains possible, no 3-consecutive survived run possible. This blocks an illegitimate inference from first return to persistence.
- H<3: three-consecutive outcome is identically zero.
- Equal q in both patches: mask-swapping changes no endpoint.
- At neutral transmission parameters, expected descendant allele copies conditional on survival equal current copies; no-input closure is absorbing.
- Cumulative events must be nondecreasing with horizon while horizon survival is nonincreasing in this absorbing-extinction closure.
- Cost frontiers fail rather than fabricate a root for unbracketed conditions or mutation/immigration.

## Reproduction / scope

`PYTHONPATH=src python -m unittest discover -s tests -v` and `python scripts/run_multigeneration.py --output results/v2_multigeneration_exact.json`; fixture checked using `scripts/verify_fixture.py`. Numerical work is exact finite-state probability propagation subject to floating arithmetic, **not** Monte Carlo or a natural data fit. Source theories remain owned by ANTO (higher-order toolkit, retry-dependent costs), EGWE (state sufficiency, interaction/allele alignment) and izu-core (selection, assurance, finite genetic accessibility and survival). v2 tests logical compatibility, not independent source replication.

**Next empirical identification requirements:** a functional floral/nonfloral module assay, confirmed parentage/allele transmission, independently measured partner interaction, time-series descendant flowering and viable-seed performance, both lost-expression and source-turnover alternatives, and whole-lineage held-out validation (aza3 and external independent clades). No SDEA result at this stage identifies trait re-evolution in Cirsium.
