# SDEA v1 analytical extension — conditional assurance/accessibility frontier

**Scope:** Exploratory follow-up to E3, added after inspecting the first synthetic E3 cells. Not a preregistered empirical finding, a natural selective cost, or a cross-species constant.

## Scientific question

At which assumed cost of reproductive assurance can a near-term occupancy benefit become insufficient to preserve the one-step probability of **both** expressed phenotype and future target-allele support? This separates present survival from later response options in one common kernel without reusing izu-core's fitted parameters (there are none here).

Let `q` be local interaction, `a_L=.1`, `a_H=.9`, `p=.12`, `N=4`, `α=.65`, `β=.5`, and assumed direct allele-weight cost `c≥0`. For a retained complete toolkit:

\[
S(a,q)=q+(1-q)a,\quad
p_s(a,c,q)=\frac{p e^{\beta(q-0.5)-ca}}{1-p+p e^{\beta(q-0.5)-ca}},
\]
\[
G(a,c,q)=1-[1-p_s(a,c,q)]^{2N},\quad
J(a,c,q)=\alpha q S(a,q)G(a,c,q).
\]

No mutation or migration, no demographic feedback, no multi-generation lineage model, and no change to q, p or N when varying a. The exact numerical procedure solves

\[
J(a_H,c^*;q)=J(a_L,c^*;q).
\]

At c=0, the two policies have equal genetic support while `S_H>S_L` for q<1, so the high-assurance J is greater. For sufficiently large c, the allele's selected fraction is asymptotically proportional to `exp(-ca)` at both assurance levels; with a_H>a_L, `G_H/G_L → 0`, eventually overwhelming the bounded positive survival ratio. **Continuity gives at least one crossing** for interior q,p and fixed finite N. This does not prove uniqueness for all other biological transition laws or with gene flow; the implementation checks numerical bracketing instead of silently assuming it.

## Computed v1 cells

| q | Assumed cost crossover c* (constructed units) | Direction for c below / above c* |
|---:|---:|---|
| 0.20 | 1.86097 | higher assurance raises / lowers J |
| 0.50 | 0.99945 | higher assurance raises / lowers J |
| 0.80 | 0.40124 | higher assurance raises / lowers J |

The threshold decreases across **these three specified q cells**. The reason is not that all reproductive assurance is harmful at high q, but that where ecological interaction already provides strong short-term survival, the incremental survival benefit of a higher a is smaller, so a fixed assumed allele-weight cost can dominate sooner. q=0 has no phenotype expression by definition; q=1 has no assurance survival advantage. Neither is assigned a positive universal crossover.

## Falsification / scope

- Change c to zero: higher assurance must increase J for q∈(0,1). If not, implementation/model accounting is incorrect.
- Below versus above each computed c* must give opposite signs and the root must equalize J to strict numerical tolerance.
- Introduce mutation/immigration or nontrivial correlated survival/trait segregation: **this version's frontier is no longer licensed as a theorem** and must be recalculated or rejected.
- A real plant's selfing cost, flower-colour genotypes, inbreeding depression, effective offspring number and ecology are unmeasured here. **This cannot estimate an empirical selection threshold.**

Reproduce: `PYTHONPATH=src python -m unittest discover -s tests -v`, `python scripts/run_benchmarks.py`; inspect `conditional_cost_frontier` in JSON. All values are exact-model numerical root calculations, not sample estimates.
