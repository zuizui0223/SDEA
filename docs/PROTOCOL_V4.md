# SDEA v4 — Exact donor–recipient pollen/seed operators and competing expression routes

**Research state: 2026-10-26.** New *constructed* model, separate from v1–v3. This document preserves all previous fixtures and inference decisions. Its numerical parameters were selected for exact synthetic counterexamples, not fitted or validated against original plants. **This model has no proof of ancestral trait re-evolution, no dated ancestry, and no wild-species route-recovery rates.**

## 1. The scientific challenge

Can two movement histories yield **identical currently observable floral reappearance probabilities**, yet leave different future genetic accessibility? Can paternal gene flow combine an existing local capacity with a donor's alternative capacity, when whole-zygote seed movement initially keeps the two capacities in different plants?

The comparison tests the interface among:

- **ANTO:** a complete or alternative latent developmental route is not equivalent to the phenotype observed today; a route introduced by *gene flow* is not intrinsically re-evolution of an ancestral route.
- **EGWE:** spatial and within-individual configuration may be needed beyond pooled marginals or presence-only phenotype to predict future descendants; a complete present-state-plus-forcing Markov closure eliminates extra history *by construction*.
- **izu-core:** pollen-mediated mating, autonomous selfing, finite diploid inheritance and the ability to retain a trait are distinct from one-step expressed phenotype. This is deliberately a **smaller and different** genetic generator, not source Model 3.

## 2. Full generative closure and timing

The two geographic locations are a **maintained donor reservoir** and a **finite recipient**. The donor adult genotype frequencies are fixed throughout all model generations. Only the recipient's genotype counts evolve. There is **one-way flow from donor to recipient**, no source demographic feedback, no local population-size changes, no stochastic plant extinction or recolonization, and no mutation. The generator conditions on viable local seed reproduction and recruits exactly `N=2` diploid plants independently into the recipient each generation. Source seeds are real *modeled zygotes produced by maternal/paternal gametes among donor adults*, not adult clones.

There are two **unlinked, biallelic loci**, each with 0, 1 or 2 copies of a functional route-associated allele. We call them `A` (ancestral-*compatible* route) and `B` (alternative-*compatible* route), without identifying real enzymes, phylogenetic ancestry, pigment hues or flower genes. A diploid has genotype `g=(k_A,k_B)`, with nine possible combinations. The recipient's current complete genotype state is the 9-vector of integer counts, summing to N. At model time zero **all plants have unobserved/absent expression**, regardless of their retained alleles. This is a stipulated initial condition, not an experimentally verified episode of trait loss.

Reproduction proceeds in this order:

1. Select local maternal parents uniformly from recipient adult genotype frequencies. Meiosis is Mendelian at each locus (transmit the functional allele with probability `k/2`); loci assort independently. No mutation or segregation distortion is allowed.
2. A seed attempt is an autonomous self with probability `s=a(1-q)` and an outcross otherwise. A locally selfed viable zygote occurs with weight `s(1-δ)` and an outcrossed viable zygote with weight `(1-s)q`. Condition on viability with `η=s(1-δ)/[s(1-δ)+(1-s)q]`. A selfed zygote draws its two gametes from the same selected local mother; outcross draws a local mother and a paternal pool. `q` is a dimensionless effective interaction index, not observed visitor count or island distance.
3. **Pollen operator:** conditional on outcrossing, replace the paternal adult gamete pool with the donor's at `m_p` (maternal parent always remains local). Pollen migration does **not** directly import complete diploid adults or seeds and cannot fertilize an autonomously selfed zygote.
4. **Seed operator:** conditional on recruited-zygote production, draw a *complete diploid donor-produced seed genotype* with probability `m_s`, otherwise recruit a locally produced zygote from step 2–3. All weights sum to one. `m_s` is a *conditional recruit-slot mixture probability*, not a demographic immigrant establishment rate. Zero-viable-reproduction cases are rejected by the current conditional kernel rather than silently adding recolonization.
5. Draw N independent recruited diploid offspring from the resulting mixed genotype distribution. Sources of this generation's **route allele entry** are tagged as local self, local outcross, donor-pollen paternal or donor-seed whole zygote. Tags are one-generation accounting; after reproduction the model does not retain permanent ancestry labels.
6. Conditional on the recruited adult genotype and its current `q`, route A starts expression with probability `α_A q` if at least one A allele is retained; route B with probability `α_B(1-q)` if at least one B allele is retained. If both initiations occur in one plant, a declared tie weight `τ` determines *which visible route is recorded for that plant*. Different plants may simultaneously express different routes. **The chosen q direction for B is an assumed environmental-response function**, not an inference that alternative pigments universally succeed with fewer pollinators.
7. The first generation in which any recipient plants express is recorded as `none`, `ancestral_only`, `alternative_only` or `simultaneous_A_and_B`. Subsequent generations keep this first-observed category. Later genetic support is tracked independently. Environmental forcing can be fixed or an explicitly specified schedule. In pulse experiments, `m_p=m_s=0` after generation 1.

**Claim limits:** Four-generation first observed expression is neither a three-generation expression streak (v2), nor lineage survival/fixation, nor genetic adaptation measured in a field population. `A` and `B` are *abstract route-associated allele states*, not confirmed molecularly ancestral or newly evolved trait outcomes. A seed's genotype can be pre-equipped without the donor flower being observed expressing it; do not label this genetic import as empirically measured migration of an already-expressed flower.

## 3. Analytical comparison: identical photos, different genetic descendants

Let the recipient initially be `A=0,B=0` in both diploids; the donor be homozygous `A=2,B=0` in both diploids. Fix q=0.6, a=.35, δ=.15, `α_A=.8`, and no B expression. Denote the viable-outcross fraction 

\[ f_o = \frac{(1-a(1-q))q}{(1-a(1-q))q + a(1-q)(1-\delta)}. \]

For **pollen-only**, each recipient recruit receives one donor A allele with probability `c=f_o m_p`; otherwise it has no A. For **seed-only**, each recipient recruit receives a donor-made homozygous A/A genotype with probability `c=m_s`. Choose `m_p=.4` and `m_s=f_o*.4=0.3250393700787402` to **match carrier probability per recruit**. This deliberately does **not** match A gene-copy counts, because that is the hidden mechanism under test.

Under either intervention, at least one carrier at generation one has exact probability

\[P(G_A(1))=1-(1-c)^N\]

and at least one visible ancestral-compatible expression has exact probability

\[P(E_A(1))=1-(1-c\alpha_A q)^N.\]

For N=2 they are **0.5444281481** and **0.2876958989** respectively **under both interventions**. But expected A allele copies at the first recipient census are **0.6500787402** after pollen versus **1.3001574803** after seed, exactly double.

After both flow operators are stopped after this first pulse, the exact conditional recipient recursion differs:

| Source pulse | First expression by generation 1 | A support at generation 1 | A support at generation 4 | First A observed by generation 4 |
|---|---:|---:|---:|---:|
| No movement | 0 | 0 | 0 | 0 |
| Donor pollen, m_p=.4 | .287695899 | .544428148 | **.269980753** | .435743340 |
| Donor seed, m_s=.32503937 (matched carrier) | .287695899 | .544428148 | **.423285993** | .499755633 |
| Donor seed, m_s=.4 (equal nominal parameter, *not* matched carrier) | .347136000 | .640000000 | .507476735 | .590893087 |

At the *same* carrier risk and photo-only distribution in generation one, the probability of retaining route A at generation four differs by **15.33 percentage points**. This is a **model-specific finite-genotype inheritance result**. It is neither a universal superiority of seed immigration nor a measured evolutionary restoration probability.

### History identifiability

Suppose prior probability of the pollen-only and matched-seed-only interventions is 1/2. In generation one, a binary observation of **any A expression** has exactly the same distribution under both: total variation TV=0, so even the Bayes-optimal identification accuracy is **50%** (chance). Observing the complete individual-genotype counts yields TV **0.5444281481**, for an optimal equal-prior identification accuracy **77.2214%**. This is an exact contrast between observation schemes under a specified closed kernel. The flow-route identity is a **historical label**, and the full genotype sample is richer present-state information. It does not prove irreducible history dependence in a fully observed natural process.

## 4. Pollen and seed differ in assembling *joint* route capacity

Hold the local population fully `A/A` and route B absent; hold donor fully `B/B` and route A absent. The pollen immigrant fertilizes an A-bearing **resident mother**, permitting **an A+B-capable hybrid in one offspring**. The seed immigrant carries a donor B-bearing offspring but **does not immediately bring the resident A allele into the same individual**.

At equal nominal `.4` pollen-only versus `.4` seed-only interventions (not an equal-gene-flow comparison):

| Source flow | B route present in recipient at gen 1 | Any individual carries *both* A and B at gen 1 | Gen-1 B expression only |
|---|---:|---:|---:|
| No movement | 0 | 0 | 0 |
| Pollen only | .544428 | **.544428** | .079909 |
| Seed only | .640000 | **0** | .165888 |
| Both pollen and seed | .835994 | .352013 | .224291 |

Thus the composite genotype can be produced immediately by pollen in this highly specified genetic scenario even when seed input is more likely to produce B carriers. By generation four, seed descendants can also form joint A+B genotypes through subsequent local mating. This is *not* a proof of ancestral-pathway rescue or real biochemical co-functionality; even two alleles with both capacities do not establish a floral trait or descendant fitness.

**Provenance:** When B starts absent from the recipient, first-generation route input is partitioned exactly into `none`, `pollen_only`, `seed_only` and `both sources represented in distinct recruits`. The last case has positive probability under simultaneous flow, so pollen- and seed-entry incidences must not be added as if disjoint counts. Source tags expire after a generation; later first-expression route labels alone do not reveal immigration history.

## 5. Environmental route filtering is an assumption to challenge

For a balanced donor carrying either A only or B only, with `.4` pollen and `.4` seed movement:

| Recipient q | First A-only expression | First B-only expression | Simultaneous A-and-B expression |
|---:|---:|---:|---:|
| .2 | .068727 | **.352294** | .018428 |
| .5 | .218283 | .218283 | .037813 |
| .8 | **.408403** | .079752 | .027434 |

This is the *mathematical consequence of choosing* A activation proportional to q and B activation proportional to (1-q), **not** novel evidence that real flowers shift from ancestral to alternative traits with environmental change. Replacing this response law with realistic, independently estimated activation functions is a prerequisite to ecological prediction.

With initial low q=.2 and a one-generation pulse, a subsequent externally imposed change to q=.8 in generations 3–4 raises the four-generation A-first cumulative probability from **.140968** (low-q throughout) to **.196479**. Future forcing therefore matters even if the starting genotype and visible phenotype match.

## 6. Adverse controls and falsifiers

- **Mass conservation:** Meiosis, mixed origin × genotype kernels, offspring multinomial distributions, within-individual A/B arbitration, population-wide first-expression labels and every DP generation must sum to one.
- **Absorbing input boundary:** A route absent from both local and source populations cannot appear without the disabled mutation/input operator. No molecular re-expression is spontaneously fabricated.
- **Operator specificity:** With local q=0, assurance=1 and self viability=1, viable locally recruited offspring can be selfed but **pollen cannot enter via outcross**. Immigrant seeds may bring route alleles, but expression of route A remains zero at q=0 under the declared A law. Present genetic support does not imply visible return.
- **Matched photo control:** Choose seed input to match **the effective number of carriers**, not the nominal migration parameter. Confirm identical generation-1 expression and support, and different genotype dosage and later retention.
- **Hybrid route:** Pollen from B-only donor and resident A-only mother may produce joint A+B offspring immediately. The matched pure seed source cannot produce joint A+B offspring until after local generations of mating.
- **Observer comparison:** Photo-only binary distributions coincide under the matched pulse, but full genotype census distinguishes the two interventions by exact total variation.
- **Route tie and absence:** If both routes are activatable within a plant, probabilistic arbitration must be exclusive; if both first appear in different plants at the same generation, record `simultaneous` instead of inventing a unique temporal winner.
- **Timing and environment:** Cumulative first-event probabilities do not decrease with additional model generations; future q is declared rather than inferred from history after observing the outcome.
- **Fails closed:** Invalid genotype counts, unsupported donor/recipient sizes, probability out of bounds, missing viable local reproduction or an inconsistent forcing schedule raise errors.

### Reproduce

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_landscape.py --output /tmp/sdea_v4.json
python scripts/verify_fixture.py results/v4_landscape_exact.json /tmp/sdea_v4.json
```

The runnable source is [`src/sdea/landscape.py`](../src/sdea/landscape.py). Frozen readout: [`results/v4_landscape_exact.json`](../results/v4_landscape_exact.json). This standard-library generator enumerates all genotype states exactly for **two recipient diploids**, conditional on the specified life cycle. It includes **neither free-flying animal movement, genuine two-patch metapopulation population feedback, sequence-based ancestry reconstruction, species-level adaptation, nor externally calibrated model parameters**. The maintained donor population and fixed recipient N are deliberately narrow assumptions, not general landscape truth.

## 7. How this becomes a biological theory test

To test the general theory rather than this construction in real flowering plants, prospective independent **within-lineage** histories require: (i) separately validated A/B floral genetic modules with functional assays; (ii) genotype, parentage and migrant source assignment at the same cohort; (iii) measured incoming seed and pollen dispersal and local viable seed production; (iv) before/after environmental partner-function measurements; (v) descendant phenotypes over at least several generations; (vi) fitness establishment and survival rather than a photos-only endpoint; (vii) competing immigration, introgression, retention, plasticity and de novo-route hypotheses; and (viii) whole-population/whole-clade held-out prediction. The Cirsium application and original biological observations remain under [aza3](https://github.com/zuizui0223/aza3).

**Claim ceiling:** SDEA v4 constructs a falsifiable mechanism where *matched visible reappearance fails to identify pollen vs seed history and genetic outcomes*, and demonstrates that source-provenance and joint genotype competence depend on movement operator. It does **not** establish a universal gain of seed flow, ancestral trait resurrection in wild plants, or a new empirical source of pigment.
