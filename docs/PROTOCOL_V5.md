# SDEA v5 — Two-way evolving metapopulation: finite genetic accessibility and colonization

**Research snapshot:** 2026-10-26 (synthetic, exploratory, after v1–v4). This protocol and its deterministic PRNG readout are **new constructions**, not a prospective ecological preregistration, a reexecution of `anto`, `egwe` or `izu-core`, or an empirically identified Ryukyu floral-recovery route. Versioned v1–v4 endpoints, numerators, fixtures and claim ceilings are preserved. Original thistle trait, ancestry, seed and molecular evidence remains in `aza3`.

## Key question and evidence separation

> When both the local population and the source population can change, go extinct and recolonize, is genetic contact sufficient to restore a *lost local option*, and is a temporarily restored option retained at the horizon?

Do not conflate (i) genetic route availability, (ii) first visible candidate-route expression, (iii) recolonization, (iv) continuous occupancy, (v) genetic fixation, or (vi) documented molecular ancestral re-expression. **Only (i)–(iii), finite-horizon occupancy and population-level matching are approximated here.** Routes A and B are *constructed* biallelic genetic capacities. They are neither empirically established ancestral pigment systems nor fitness-measured alternatives.

## Differences from the v4 static-donor model

Both patches now have separate **changing** populations of 0–2 diploid plants, with two unlinked Mendelian loci `(A copies, B copies) ∈ {0,1,2}²`, independent local ecological-state chains, optional expression→ecological-feedback, independent catastrophes, and bidirectional pollen and seed exchange. Parentage and reproduction draw from **both patches' previous generation**, making all incoming flow contingent on surviving local source parents. No exogenous fixed donor pool, whole-genome lineage history, unlimited species pool, mutation or creation of a new route is allowed.

There are still deliberately substantial limitations: only two patches, no absolute geography or seasonal q, no explicit pollinator individuals, no patch carrying capacity above two, no explicit pollen-limited ovule counts, no genotype-dependent fitness or parental age, no resource budget on **donor exported propagules**, and no independently calibrated reproductive parameters. Donor seed exports are drawn independently per receiving recruit slot, conditional on the donor's viable-seed probability; source seed production is not decremented. Environment q is an index of effective interaction, **not pollinator diversity**.

### Life cycle (simultaneous source/recipient update)

For each occupied patch `i`, its current q is either .25 or .75; viable autonomous selfing has attempt share `s_i=a(1−q_i)` and viability factor `(1−δ)`. Local outcross and immigrant-pollen fertilization have combined viable weight `(1−s_i)q_i` when an admissible father exists. Maternal parents are local; a *distinct* local father is required for local outcrossing. **Pollen migration** replaces only the paternal gamete, with probability `m^P_i` if the other patch has at least one parent. A patch with no mother cannot be fertilized by incoming pollen.

For every one of two recipient recruitment slots, a donor-produced **viable** seed is offered first with probability `u_i=m^S_i S_j^{local}`, conditional on the other patch having seed-producing parents. If no donor seed arrives, the recipient attempts a viable local reproductive event. The donor seed is **formed by Mendelian segregation and the donor's own local selfing/outcrossing**; it is not an adult clone or pollen-only donation. Failed slots remain vacant. An independent patch-local catastrophic shock after recruitment can erase all recruits. Crucially, **seed immigration can colonize an initially empty patch or recolonize one that lost occupancy**, but pollen cannot do either without resident mothers. Nonempty patches can be occupied at a later generation with fewer than two plants. No artificial replacement of empty adults is introduced.

Source and recipient are updated **synchronously** from their previous adult cohorts (not from newly produced migrants); seed and pollen flows may run in both directions. q changes between low/high after recruitment with `P(high'|high)=.80`, `P(high'|low)=.12`, optionally modified by the **previous** observed A/B expression (±.12 in the feedback stress cell). Both alleles assort independently without mutation. At model time zero phenotype is stipulated absent/unobserved; a newly recruited plant with A may initiate A with probability `.75 q`, and a plant with B may initiate B with `.75(1−q)`. If both initiate within one plant, a `.5` tie rule chooses **one** recorded route. Different plants may simultaneously show different routes. The q-activation law is a **hypothesis-defining convenience**, not a biological inference from flowers.

### One-generation exact controls

The code enumerates an **exact, normalized probability distribution per recruit slot** for `empty / local self / local outcross / immigrant paternal pollen / immigrant whole donor seed`, × every one of the nine two-locus diploid offspring genotypes. The longer trajectory is **Monte Carlo**, not exact enumeration across all two-patch lineages. Each PRNG replicate is one entire metapopulation trajectory; its four recruited slots or generations **are not independent ecological replicates**.

With an empty focal patch, two A-homozygous source plants and no catastrophes in the focal patch, let `S_d` be source viable local reproduction per slot. The exact probability of first colonization in one generation is

\[
P(C_1) = (1-\kappa_0) [1-(1-m^S_0 S_d)^N], \quad N=2.
\]

For `q_d=.75`, `a=.45`, `δ=.20`, `S_d=.755625`, `m^S_0=.35` and `κ_0=0`, this is **0.4589937803**. Pollen rate alone cannot colonize an empty patch under this reproductive mechanism (exact probability zero), even if pollen can cross an established population.

A second exact one-step counterexample fixes patch 0 to two `A=2,B=0` diploids, patch 1 to two `A=0,B=2` diploids, patch 0 q=.25, patch 1 q=.75, assurance `.45` and depression `.20`. With no immigration: final generation-1 focal A support **0.681481** and occupancy **0.681481**. With seed input `.35` from B-only donor: A support falls to **0.538165** while occupancy rises to **0.827679**. With pollen input `.45` only, the probability of at least one individual retaining **both** candidate capacities at generation 1 is **0.143508**; with only B-seed import it is exactly zero at this first step. These are not universal pollen/seed rankings and do not imply adaptive superiority.

## Complete finite-population trajectories and reporting contract

Nine **exploratory** scenarios use N=2, H=7 generations and **9,000 independently seeded full-landscape trajectories each**. The integer seeds/parameters are frozen in `benchmark` and saved per scenario in `results/v5_metapopulation_mc.json`. All endpoints are unconditional on survival and include extinct metapopulations in the denominator.

Two initial conditions:
- `empty`: patch0 initially empty, patch1 with two `(A=2,B=0)` plants. Initial ecological indices `(low,high)`.
- `distinct`: patch0 has two `(A=2,B=0)` plants; patch1 has two `(A=0,B=2)` plants. Same initial ecology.

Rates are per destination, so a `.45` pollen rate means each outcross in **either direction** uses the other patch for the father with that probability (when present); `.35` seed rate means an established seed-pressure probability scaled by the other patch's viable local seed probability. The same numbers do not equate migrant genomes, pollen grains, seed counts or ecological effects. Scenarios vary none / pollen-only / seed-only / both, an additional **0.32 catastrophe probability in patch1** to make the evolving source fragile, and a distinct A/B expression→q feedback stress cell. No history outcome is used to fit original parent-model parameters.

### Registered names distinguish colonization and recolonization

- `patch0_colonized_initial_vacancy`: patch0 **started empty** and obtained its **first-ever** recruited adult through immigrant seed.
- `patch0_recolonized_after_extinction`: patch0 was **previously occupied**, then experienced an actual occupancy extinction, then obtained an adult from immigrant seed. This event can occur after initial colonization or for initially occupied patches.
- `patch0_A_supported`: has at least one A allele at generation 7, irrespective of flower expression or demographic viability beyond H.
- `patch0_joint_in_individual`: any single individual carries both A and B route-associated alleles at H, **not biochemical compatibility or simultaneous visible routes**.
- `patch0_joint_ever`: at least one A+B individual appeared in that patch by H, including those lost before H.
- `both_patches_occupied`: each patch has ≥1 recruited adult at H; does not guarantee long-term persistence or retained route alleles.
- `global_A_absent`, `global_B_absent`: absorbing gene-copy absence across both patches (given no mutation/external reservoir).
- `patch0_..._gene_entered_via_pollen/seed`: first-step *gene-copy entrance when the recipient lacked that route allele just before reproduction*; **not a permanent ancestry marker or visible floral origin claim**. Source tags are destroyed when descendants reproduce; later genomic genealogy is not inferred.
- `patch0_any_A_ever` / `...B...`: at least one visible expression initiation of the respective route by H; an already equipped immigrant seed might or might not express.

For each binary outcome: count, probability and Wilson 95% **Monte Carlo sampling interval**. Exploratory cross-scenario differences use a normal 95% interval from two independent 9,000-history ensembles, **without multiple-comparison correction**. These intervals quantify uncertainty from the artificial RNG; they are not ecological effect confidence intervals or source-independent biological replication.

## Numerically executed H=7 results

| Initial / intervention | First colonization of empty patch0 | Recolonization after recorded extinction in patch0 | Both patches occupied at H | A+B individual in patch0 at H |
| --- | ---: | ---: | ---: | ---: |
| Empty / no exchange | 0 | 0 | 0 | 0 |
| Empty / pollen-only | 0 | 0 | 0 | 0 |
| Empty / seed-only | **.67000** | .08600 | .03089 | 0 |
| Empty / seed-only, donor catastrophe .32 | **.59856** | .03189 | .00700 | 0 |
| A patch0 / B patch1, no exchange | n/a | 0 | .00011 | 0 |
| Distinct, pollen-only | n/a | 0 | .00567 | .00489 |
| Distinct, seed-only | n/a | .19567 | .04100 | .00867 |
| Distinct, both exchange | n/a | .21600 | **.10267** | **.03200** |
| Distinct, both + expression-q feedback | n/a | .21944 | .10200 | .03211 |

For empty→seed-only, first colonization by H was **67.0%** but focal occupancy **at H** was only **6.37%**. A successful colonization event **is not** a sustained metapopulation rescue. Catastrophe of the evolving source reduced first colonization by **7.144 percentage points**, with an illustrative 95% MC interval **[5.741,8.548] pp**. More fragile donor conditions reduced patch0 A support at H from **6.367% to 2.444%**, in these small synthetic populations.

In the distinct-route treatment, pollen-only assembled A+B genetic capacity at some point by H in **24.57%** of focal trajectories, yet only **0.489%** still had a dual-capable individual at H. Under both pollen+seed movement, **3.200%** had such an individual at H and **10.267%** retained occupied adults in both patches; effects reflect the *entire declared reproduction, movement and extinction kernel*. Seed movement sometimes **replaces locally A-bearing recruits** with B-only donor offspring: the exact first-step counterexample demonstrates an occupancy benefit **coexisting with an ancestral-allele-support penalty**, even without a direct selection cost on A.

The expression-q feedback cell does not establish a universally favorable recoupling process; high-q focal occupancy changed from **.35044 to .33189** under a chosen A-positive/B-negative feedback law, and the sign depends on the assumed activation/feedback rule. Do **not** call this empirical ecological compensation or universal pollinator feedback.

## Falsification / negative-control gates

1. Exact 9-genotype Mendelian cross matrices and all source-tagged slot distributions sum to one (every declared valid birth/arrival possibility is accounted for).
2. A route absent from **both patches** cannot arise (mutation disabled). Global loss is absorbing; introgression from the other patch is not a de novo genetic return.
3. If both patches are extinct, **neither seed nor pollen exchange** creates organisms without an external reservoir.
4. In an empty recipient, pollen-only never colonizes; viable seed input can. A later true recolonization requires an observed earlier occupancy-loss transition, not merely the initial empty state.
5. Pure pollen from B-only father into A-only mother may create A+B offspring immediately; pure donor-produced B-only seed by itself cannot. At a later generation, mating between A and B plants may assemble the joint capacity under seed migration too.
6. Seed input from B-only source into A-only recipient can increase occupancy **and decrease A support** in one step. This refutes an unconditional beneficial-connectivity assertion for that receptor-genetic endpoint.
7. q=0 blocks pollen-mediated outcrossing under the declared viability law even if donor adults persist; viable selfing may remain. Arbitrary q-feedback does not count as an observed visitor response.
8. Every recorded genetic route entrance has a live donor in the **previous** generation and an explicitly tagged pollen or seed origin. Route labels do not diagnose ancestral gene recovery in real taxa.
9. Repeated runs with the same seeds and config reproduce identical numerical results. CI regenerates the fixture, compares all fields with the established v1 tolerance checker, and retains v1–v4 checks.

## Source interface and next step

- [ANTO](https://github.com/zuizui0223/anto): actual hypothesized functional genetic toolkit and ancestral/alternative re-expression **remain unmeasured** by this 2-locus toy. Route-specific failure cost or real lineage restoration is not implemented here.
- [EGWE](https://github.com/zuizui0223/egwe): local interactions and **within-individual** genetic alignment predict different futures; the q feedback here is imposed rather than inferred from the EGC interaction network.
- [izu-core](https://github.com/zuizui0223/izu-core): finite genetic recruitment, reproductive assurance and survival are separated. This is **not** full Model 3 and no parent implementation/registered Q1–Q4 values are edited.
- [aza3](https://github.com/zuizui0223/aza3): retains the original Cirsium biological data/field experiments; molecular functional modules, genomic parentage, actual seed and pollen rates and reproduced flowering success are still needed.

The most important open questions are **true donor-limited propagule budgets**, more than two patches/individuals, explicit visitor dispersal and effective pollen transfer, tradeoffs between route A and B, repeated expression across *descendants* and independent functional/genomic calibration. Until that work and externally held-out natural taxa exist, the scientific claim is **a set of reproducible structural counterexamples and model-dependent estimates, not a universal evolutionary law**.

### Run

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_metapopulation.py --output /tmp/sdea_v5.json --replicates 9000
python scripts/verify_fixture.py results/v5_metapopulation_mc.json /tmp/sdea_v5.json
```

The code uses Python standard library only. Cross-Python floating-point comparisons use tolerance from the preexisting `verify_fixture.py`, while PRNG streams are pinned. `test_metapopulation.py` contains **27 new deterministic/negative-control tests**. Numerical intervals reflect only finite synthetic trajectory sampling.
