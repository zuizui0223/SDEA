# SDEA v1 — preread exact theoretical validation protocol

**Frozen choices: 2026-10-26.** This specification governs the initial constructed exact benchmark; it is *not* a prospective field preregistration or a new biological simulation cohort. All future changes must version the model rather than replacing unfavorable results.

## Question and conceptual separation

Can history, latent functional toolkit composition, spatial ecological alignment and finite allele transmission distinguish the future of systems with the same currently visible phenotype? The required distinctions are: visible absence; functional joint capacity; first expression; continued genetic support; one-step occupancy/survival; *long-term* ecological establishment (not tested).

The model evaluates a **one-step, conditional probability**, not a multigeneration lineage process. Current phenotype is white/absent for all constructed founders. An integer mask 0–7 encodes three latent module-presence bits in the *same* lineage. Full expression requires mask 7. In the assay experiment, histories A/B set different prior mask distributions. The model explicitly assumes complete observed state plus future forcing law is Markov-sufficient; whether real phenotypes satisfy that assumption is open.

## Exact model and assumptions

See [`src/sdea/benchmark.py`](../src/sdea/benchmark.py) for executable functions. All endpoints require stated local q in [0,1], activation α, assurance a, diploid size N and target allele frequency p. Expression E=1{mask=7} αq. Survival gate S=q+(1-q)a. Target allele's relative weight is w=exp[β(q−0.5)−ca]. Conditional target allele fraction after selection is pw/(1−p+pw); optional one-way mutation and migrant replacement follow; G=1−(1−p_eff)^(2N). The joint event J=ESG assumes independence of the three gates conditional on this synthetic full state. This construction **does not implement Mendelian selfing, breeding, genotype density, source ecological feedback or branch survival**. The q index is not raw visitor count. The genetic-support gate concerns a *target allele*, while the mask concerns a distinct *functional toolkit*: one does not stand in for the other.

Default parameters: α=0.65, a=0.30, p=0.12, N=4, β=0.50, c=1.00, mutation and migrant fraction 0. Assess only the explicitly prespecified cells; no statistical significance/uncertainty is appropriate for exact values.

## E1: hidden toolkit / historical proxy

History A (even parity): 000,011,101,110, each 1/4. History B (odd parity): 001,010,100,111, each 1/4. Label priors 1/2 each. Single and pairwise module marginals must match (1/2 and 1/4); triple-functional toolkit probability must be 0 vs 1/4. All founders express no visible trait at time t. A constant q=.65 is shared. The future one-step event is **joint**, not a historical re-evolution. Use Bayes-optimal expected Brier risk under visible-only, history label, noisy three-module assay, assay+history, full true mask, mask+history. Noisy assay assumes independent sensitivity .90/specificity .95 for each bit; these are stipulated, not validated.

**Pass / falsifiers:** history must improve visible-only one-step prediction; noisy assay+history must improve over noisy assay; full mask+history must match full mask to numerical tolerance. Under the null that both history labels have the same uniform 8-mask distribution, history must add zero information. With sensitivity=specificity=.5, assay must be uninformative; perfect assays must equal full-mask conditioning.

## E2: ecological-genetic spatial alignment

Two patches have q=.90 and q=.20, one mask 7 and one mask 0. Both compared systems have exactly the same pooled mask, observed phenotype, q, site number, allele p, N and assurance. Only the mask–q association differs (complete toolkit at q_high vs q_low). Report the one-step J at each assignment. Neither a total-marginal nor partner-richness state may substitute for the joint patch mapping.

**Pass / falsifier:** aligned and misaligned differ, while the equal-q (.60/.60) negative control yields exactly equal outcomes. This does not guarantee an eventual functional-loss effect in a longer-lived system.

## E3: present occupancy vs future allelic support

For q in {.20,.50,.80} and tradeoff penalty c in {0,.5,1,2}, contrast assurance .10 vs .90, with all other variables fixed. Report separately S, G and J. The tradeoff penalty represents an **assumed direct relative-transmission penalty on the high-trait allele**, introduced solely to exhibit a mechanism; it is not inferred from izu-core and must not be read as measured inbreeding depression or selfing-associated loss.

**Pass / falsifiers:** at c=0, more assurance raises S and J while G is unchanged. At q=.80, c=2, more assurance raises S but lowers G and J. Other cells can differ in sign: a coexistence demonstration, not a universal demographic tradeoff. A real test would intervene on assurance independently of genetic support and separate parentage, offspring number and selection.

## E3b: exploratory conditional cost crossover (added after E3 inspection)

This is **a transparently post-E3 analytical extension**, not a preregistered new confirmatory test. Under 0<q<1, 0<p<1 and no mutation or migration, for assurance aL=.1 and aH=.9 define c*(q) as a root of J(aH,c)-J(aL,c)=0. The model uses a fixed N=4, α=.65 and other default parameters. Run a bracketing bisection on c∈[0,8] for q in {.2,.5,.8}; explicitly fail if no change of sign exists. At c=0, S_H>S_L and G_H=G_L. When c is sufficiently large, G_H/G_L→0, so at least one crossing exists for these interior parameter settings. Check J_H-J_L signs immediately below/above each computed c*. This identifies a *conditional* tradeoff frontier, not a universal threshold or the evolution of assurance.

## E4: no-input closure and restoration route

Set initial target allele frequency p=0. Without mutation or immigration, G must be exactly zero for any positive finite N. With independent input mutation_rate=.02 or immigrant_fraction=.02 (incoming frequency 1), G must be positive. A missing target allele is not identical to loss of the three-module toolkit, so these tests are about allele-support *logic* only, not ancestral gene-module recovery.

## Identification, generality and independence gates

1. The history label is **predictive proxy**, not intrinsically causal. Any apparent residual historical effect after full functional state conditioning prompts a check for missing measured coordinates or changed future forcing—not proof of non-Markov evolutionary memory.
2. `joint` is not a verified stable multi-generation phenotype, a lineage nonextinction probability or a proxy for actual flower colour. Model gates are assumed independent and can be causally wrong in organisms.
3. The EGWE matched-summary result and ANTO parity example motivated constructions, but this new fixture is **not an independent empirical confirmation** of those papers. izu-core's 160×160 registered outcomes are not fitted, pooled or overwritten.
4. For biological transport: preregister taxa/populations; infer dated independent histories; verify same-individual expression, functional module compatibility and nonfloral regulation; measure partner interaction and mating/parentage; genotype descendants over >1 generation; test new clades/regions with strict held-out units. No admitted biological datasets currently satisfy this entire contract in SDEA.
5. **Reject universal claim** if any independent natural system violates a universal sign or if state representation remains insufficient after the proposed full-state assays. Do not add post-hoc qualifiers until every discrepancy vanishes. Use explicit model comparison and uncertain measurement corrections.

## Reproduction and audit

Run `PYTHONPATH=src python -m unittest discover -s tests -v` then `python scripts/run_benchmarks.py --output results/v1_exact_benchmarks.json`. CI verifies matching structure and checks finite floating values to a tight cross-Python tolerance (relative 1e-12, absolute 1e-13). Exact probabilities have no sampling confidence interval; numbers are not estimated ecological quantities. At present v1 covers one step and one constructed law only. The untested future target is a **multi-generation, spatially explicit plant lineage state model** with independently validated biological priors and a lineage-level held-out empirical design.
