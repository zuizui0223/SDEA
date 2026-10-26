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

Results are in [`results/v1_exact_benchmarks.json`](results/v1_exact_benchmarks.json). Cross-Python checks compare exact structure and floating-point agreement within specified numerical tolerance. Full model decisions, falsifiers and evidence tiers are in [`docs/PROTOCOL_V1.md`](docs/PROTOCOL_V1.md). Source ownership and pinned documents are in [`docs/SOURCE_REGISTRY.md`](docs/SOURCE_REGISTRY.md).

## Boundary / next step

SDEA owns **cross-model definitions, matched-state counterexamples, validation contracts and eventual held-out cross-system tests**. It does not silently copy, recalibrate or reinterpret the parent repositories' simulations. Actual Cirsium molecular/experimental work belongs in [aza3](https://github.com/zuizui0223/aza3); natural fragmentation synthesis belongs in [egwee](https://github.com/zuizui0223/egwee). Next research gate: a common empirical observation design that independently measures functional module compatibility, ecological partner state, finite allelic support and descendant expression across time. Until then, *universal, causal and field-confirmed* claims are **NOT established**.
