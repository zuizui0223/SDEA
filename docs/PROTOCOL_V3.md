# SDEA v3 — Parent-level diploid inheritance and predictive state sufficiency

**2026-10-26, synthetic mechanism challenge (exploratory).** Version 3 is a deliberately **different generative closure** from v1's independent allele-copy binomial gate and v2's pooled-frequency binomial recursion. Neither their results nor preread claims are retroactively reclassified. No empirical data, genetic fitness or flower-colour reversion is measured. The [single-locus exact operator](../src/sdea/mendelian.py) is not `izu-core` Model 3.

## Central challenge to portability

> Is population-wide allele frequency a sufficient present descriptor of future genetic accessibility when autonomous self-fertilisation couples gametes from the **same diploid parent**?

**Two equally silent, complete-toolkit synthetic populations:** four diploid founders; one target locus, a/a, a/A and A/A; eight total allele slots; target allele A total copies k=2. Two placements:

- **HET**: `(aa=2, Aa=2, AA=0)`.
- **HOM**: `(aa=3, Aa=0, AA=1)`.

They match census, allele frequency 2/8, initial visibly absent phenotype, all abstract functional modules, baseline q, founder genotype counts' marginal *copy* total, and future stochastic forcing laws. They differ in **within-individual** allele arrangement. The stated history labels are merely names for these present genotypes, not actual ancestor reconstructions.

## Model construction

For each occupied generation, N=4 parents have genotypes `g∈{0,1,2}` target copies. A genotype has reproductive weight `exp[g(β(q−0.5)−c a)/2]` (an explicit assumption, not a measured selection coefficient). Sample each contributing parent with replacement according to normalized counts times weights. Let `u_g = g/2` be probability that a Mendelian gamete carries A. A one-way mutation rate μ (default zero) changes it to `u'_g=u_g+(1-u_g)μ`.

A synthetic potential offspring attempts autonomous selfing with probability `s=a(1−q)`; otherwise it attempts outcrossing. An outcross produces a viable seed with probability `q`; an autonomous selfed offspring survives a declared viability factor `1−δ` (δ=.2 by default). Thus the *single-stage population reproduction gate* is

\[
S=(1-s)q+s(1-\delta),\qquad \eta=\frac{s(1-\delta)}{S}\quad(S>0),
\]

where η is the conditional selfing fraction **among viable offspring**, not the raw attempt fraction. No viable offspring gate (`S=0`) means extinction. **Conditional on this one reproduction gate succeeding**, exactly N new diploids are drawn independently from the viable offspring genotype distribution; hence finite multinomial recruitment. The same conditional recruitment process can sample a parent repeatedly. There is no explicit individual lifespan, full birth–death abundance, biparental self-incompatibility, linkage, diploid mate-pair exclusion, gamete limitation beyond q, or population-dependent recruitment ceiling. Independent founders/recruits **are not** all independent evolutionary lineages.

The focal zygote distribution mixes (i) **outcross**: two independent weighted-pool gametes (`Binom(2, mean u'_g)`), and (ii) **selfed**: two gametes from the same sampled parent (`mixture of Binom(2,u'_g)`), with probabilities 1−η and η. Optionally, each recruited zygote is wholly replaced by an A/A immigrant with probability `i` (default zero). **Immigration occurs only in occupied systems following a successful reproduction gate**; there is no colonization of extinct populations. It is an explicit genotype input, not pollen dispersal calibrated from `izu-core`.

Local q switches between .25 and .75 with `P(high'|high)=.85`, `P(high'|low)=.15`, optionally +0.12 following a successful visible expression event (clipped at 1). Expression is an independent Bernoulli gate of probability `0.7q` if toolkit present and ≥1 A copy; initially all groups are *not visibly expressing*. A first expression may occur even if the subsequent reproduction gate fails; 3 consecutive **survived** expression updates are separately tracked. Extinction is absorbing, and conditional population size N=4 is not demographic growth.

## Derived one-generation identity (reusable theory, restricted to this kernel)

With weighted parental gamete mean \(\bar u=E_\pi[u_g']\), variance \(V=\mathrm{Var}_\pi(u_g')\), immigration probability \(i\), and viable selfed fraction \(\eta\), exact Mendelian segregation and the mixture law yield

\[
\begin{aligned}
P(aa\mid\text{successful reproduction}) &= (1-i)\big[(1-\bar u)^2+\eta V\big],\\
P(k_{t+1}=0\mid\text{successful reproduction}) &= P(aa)^N.
\end{aligned}
\]

*Proof:* Outcross gametes are independent from the same weighted parent pool, hence `P(aa|outcross)=(1−E[u'])²`. Selfing samples one parent before the two gametes, hence `P(aa|self)=E[(1−u')²]=(1−E[u'])²+Var(u')`. Weight by η, then multiply by `1−i` because immigrants here are always A/A. With conditional IID N recruits, total allele absence is all N genotype-0 offspring. This identity requires this genotype/fitness/mating/immigrant closure and is **not** a general biological extinction law.

The important identification result is: **holding \(\bar u\) fixed does not hold future allele-loss risk fixed when \(\eta V\) differs**. For this one-generation endpoint, reporting weighted gamete mean *and* parent-level variance and the realized viable self fraction closes the missing state coordinate. Higher-order multilocus assortment, individual pollen limitation and long-term spatial history are not settled by that sufficiency result.

## Exact test cases and negative controls

The fixed q=.75 and δ=.2 case has the following **conditional upon successfully producing recruits** one-generation allele disappearance probabilities:

| Assurance a | Viable self fraction η | Two Aa founders (HET) | One AA founder (HOM) | HOM − HET |
|---:|---:|---:|---:|---:|
| 0.0 | 0.000000 | 0.100112915 | 0.100112915 | 0 |
| 0.5 | 0.132231 | 0.106127443 | 0.118965339 | +0.012837896 |
| 1.0 | 0.262295 | 0.112303818 | 0.139990335 | +0.027686517 |

The **zero-selfing negative control** removes the genotype-configuration difference when weights are neutral. Under full selfing at q=0, δ=0, the conditional disappearance probabilities become `(5/8)^4 = 0.152587890625` versus `(3/4)^4 = 0.31640625`. Thus the sign is not a consequence of manually writing a different allele count: both founder copies remain k=2.

At a=0 to a=1 with **no direct target-allele cost**, the reproduction-gate success in the q=.75, δ=.2 cell increases 0.7500→0.7625, but the probability of **conditional genetic support** in HOM decreases 0.899887085→0.860009665. This contradicts any universal inference from v1–v2 that a cost-free assurance policy must leave allele support unchanged: those versions explicitly held independent-copy inheritance fixed when imposing the policy. Here assurance also changes mating structure. It does **not** demonstrate a negative *unconditional* demographic benefit, because S rises, nor natural selfing depression.

A posthoc **selection perturbation** (β=.6, assurance=.7; q=.75) reverses the HET/HOM conditional disappearance ordering (HET ≈0.09837856 > HOM ≈0.09607110). The differing sign arises because adult genotypes receive different assumed reproductive weights; the neutral comparison is not a species-invariant rule. This treatment was **not** an independently preregistered confirmatory test.

The exact H=6 extension tests non-equivalence of (i) expression at least once; (ii) ever 3 survived consecutive expressions; (iii) still occupied at H after that streak; (iv) living with ≥1 A; (v) living fixed for A. They are separate events and the latter is not guaranteed by an earlier expression event. Other adverse controls: no genotype A and no input makes target loss absorbing; explicit gamete mutation or A/A immigrant input can restore allele support before total extinction; absent latent toolkit blocks expression but not genetic persistence; q=0 with assurance=0 forces immediate extinction; genotype probabilities and occupied/extinct mass conserve to numerical tolerance; neutral inheritance preserves the expected **number** of A copies conditional on successful recruitment, though fixation and loss remain possible.

## State prediction check (EGWE interface)

With 50:50 HET/HOM prior, equal current observed phenotype, allele count k=2, q and a, forecast the **conditional next-generation A-loss** outcome. The Bayes-optimal Brier risk is

\[ B_{\mathrm{coarse}}=p(1-p),\quad B_{\mathrm{genotype}}=\tfrac12[p_{HET}(1-p_{HET})+p_{HOM}(1-p_{HOM})]. \]

With a=1, the observed genotype arrangement reduces expected Brier risk by `0.0001916358128` (an exact hypothetical forecasting gain, not empirical accuracy). With a=0 and neutral weights, gain is exactly zero. Once genotype is known, the artificial history label adds no further information under this fully specified transition law. This is **constructive** state insufficiency of frequency-only summaries, not a field demonstration of persistent evolutionary memory.

## Source ownership and falsifiers

- [ANTO](https://github.com/zuizui0223/anto): joint latent functional modules and route-specific re-expression; v3 keeps a **single boolean toolkit**, not ancestral vs alternative metabolic pathways or failure damage.
- [EGWE](https://github.com/zuizui0223/egwe): state representation and cross-layer alignment; v3 explicitly tests within-parent association and a q feedback **hypothesis**, not the EGC simulator operator.
- [izu-core](https://github.com/zuizui0223/izu-core): finite inherited reproduction, assurance, trait reachability vs occupancy. v3 has real **one-locus Mendelian segregation**, but not full Model 3 assortative matching, three-locus genetics, individual flower allocation or pollen dynamics.
- [aza3](https://github.com/zuizui0223/aza3): sole active source for original plant molecular/seed experiments; none are copied into this synthetic fixture.

**Reject/qualify the mechanism if:** (a) the one-step identity disagrees with direct genotype enumeration; (b) probability mass is not conserved; (c) matched-copy pure outcross founders differ under neutral weights; (d) parent-history continues to improve forecast after parent genotypes and future forcing are supplied; (e) observed natural taxa fail precommitted independent-lineage holdouts, or the required parentage and genotype assays are unavailable. Simulation replication is *not* independent empirical corroboration.

## Reproduce

Python 3.11+, only standard library:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python scripts/run_mendelian.py --output results/v3_mendelian_exact.json
python scripts/verify_fixture.py results/v3_mendelian_exact.json results/v3_mendelian_exact.json
```

GitHub Actions recomputes to a temporary path and compares with the checked-in frozen JSON using cross-Python floating tolerance. v1/v2 fixtures and the existing Draft PR remain unchanged by v3.
