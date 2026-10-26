"""SDEA v3: exact finite diploid single-locus mating-structure counterexample.

All rates are constructed and uncalibrated. This is NOT izu-core Model 3 or
Mendelian multi-locus genomic simulation. Adult reproductive weights, gamete
transmission, autonomous selfing, pollen-dependent outcrossing, viable-seed
conditioning and finite recruitment are explicit. Recruitment is with
replacement conditional on successful reproduction and fixed population size.
"""

from collections import defaultdict
from dataclasses import dataclass, replace
from functools import lru_cache
from math import comb, exp, isfinite


def _unit(name, value):
    if type(value) not in (int, float) or not isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{name} must be finite and within [0, 1]")


@dataclass(frozen=True)
class Genetics:
    diploids: int = 4
    q_low: float = .25
    q_high: float = .75
    high_stays_high: float = .85
    low_becomes_high: float = .15
    expression_feedback: float = 0.
    activation: float = .70
    assurance: float = .50
    inbreeding_depression: float = .20
    selection_strength: float = .0
    direct_allele_cost: float = .0
    mutation: float = 0.
    immigration: float = 0.
    consecutive_required: int = 3

    def __post_init__(self):
        if type(self.diploids) is not int or not 1 <= self.diploids <= 12:
            raise ValueError("diploids must be 1..12")
        if type(self.consecutive_required) is not int or not 2 <= self.consecutive_required <= 8:
            raise ValueError("consecutive_required must be 2..8")
        for name in ("q_low", "q_high", "high_stays_high", "low_becomes_high",
                     "expression_feedback", "activation", "assurance", "inbreeding_depression",
                     "mutation", "immigration"):
            _unit(name, getattr(self, name))
        if self.q_low > self.q_high:
            raise ValueError("q_low cannot be above q_high")
        for name in ("selection_strength", "direct_allele_cost"):
            value = getattr(self, name)
            if type(value) not in (int, float) or not isfinite(value) or not 0 <= value <= 100:
                raise ValueError(f"{name} must be finite within [0, 100]")


def _adult(counts, par):
    if not isinstance(counts, tuple) or len(counts) != 3 or any(type(n) is not int or n < 0 for n in counts):
        raise ValueError("genotype counts must be a tuple of 3 nonnegative integers")
    if sum(counts) != par.diploids:
        raise ValueError("genotype counts must sum to diploids")


def allele_copies(counts):
    return counts[1] + 2 * counts[2]


def all_genotype_counts(n):
    return tuple((n0, n1, n - n0 - n1) for n0 in range(n + 1) for n1 in range(n - n0 + 1))


def gamete_offspring_distribution(counts, qi, par=Genetics()):
    """(viable reproduction gate, child genotype 0/1/2 probabilities,
    gametic mean, effective selfing among viable produced seeds).

    Adult weight exp(g/2 * [beta(q-.5)-c*a]); parent sampled with replacement.
    In outcrossing, paternal and maternal gametes are independent draws from
    that weighted adult pool. Selfing uses two meioses of the SAME selected
    parent, i.e. the induced zygote distribution is a mixture of Binom(2,g/2).
    Autonomous selfing attempt share s=a(1-q); outcross success q, and
    selfed-seed viability 1-delta. At most one conditional successful family
    event occurs per generation; then N adults are sampled IID among viable
    zygotes. Population counts change by finite multinomial recruitment.
    """
    _adult(counts, par)
    if qi not in (0, 1):
        raise ValueError("ecology index must be 0 or 1")
    q = (par.q_low, par.q_high)[qi]
    direction = par.selection_strength * (q - .5) - par.direct_allele_cost * par.assurance
    weights = [counts[g] * exp(.5 * g * direction) for g in range(3)]
    total = sum(weights)
    mother = [w / total for w in weights]
    allele_freq = sum(mother[g] * g / 2 for g in range(3))
    mutant_freq = allele_freq + (1 - allele_freq) * par.mutation
    outcross = ((1 - mutant_freq)**2, 2*mutant_freq*(1 - mutant_freq), mutant_freq**2)
    selfed = [0., 0., 0.]
    for g in range(3):
        u = g / 2
        mutated_gamete = u + (1 - u) * par.mutation
        selfed[0] += mother[g] * (1 - mutated_gamete)**2
        selfed[1] += mother[g] * 2 * mutated_gamete * (1 - mutated_gamete)
        selfed[2] += mother[g] * mutated_gamete**2
    self_attempt = par.assurance * (1-q)
    outcross_viable = (1-self_attempt) * q
    self_viable = self_attempt * (1-par.inbreeding_depression)
    survival = outcross_viable + self_viable
    if survival:
        child = [(outcross_viable*outcross[g] + self_viable*selfed[g])/survival
                 for g in range(3)]
        # Immigrants have homozygous target genotype, conditional on occupancy;
        # they cannot recolonize the absorbing extinct state.
        child = [(1-par.immigration)*child[g] + (par.immigration if g == 2 else 0)
                 for g in range(3)]
    else:
        child = [0., 0., 0.]
    if survival and abs(sum(child)-1) > 1e-12:
        raise AssertionError("offspring genotype kernel not normalized")
    return survival, tuple(child), allele_freq, self_viable / survival if survival else 0.


@lru_cache(maxsize=15000)
def _next_counts(counts, qi, par):
    survival, child, _, _ = gamete_offspring_distribution(counts, qi, par)
    if survival == 0:
        return ()
    N = par.diploids
    out = []
    for new in all_genotype_counts(N):
        n0,n1,n2 = new
        prob = comb(N,n0)*comb(N-n0,n1)*child[0]**n0*child[1]**n1*child[2]**n2
        if prob:
            out.append((new,prob))
    if abs(sum(p for _,p in out)-1) > 5e-12:
        raise AssertionError("finite recruitment kernel not normalized")
    return tuple(out)


def _report(states, extinct, extinct_seen, extinct_three, generation, par):
    alive = sum(states.values())
    def mass(fn):
        return sum(w for (n,qi,run,seen,three),w in states.items() if fn(n,qi,run,seen,three))
    return {
        "generation": generation,
        "alive": alive,
        "extinct": extinct,
        "first_expression_ever": extinct_seen + mass(lambda n,q,r,s,t: s),
        "three_consecutive_survived_ever": extinct_three + mass(lambda n,q,r,s,t: t),
        "alive_and_ever_three": mass(lambda n,q,r,s,t: t),
        "alive_and_current_three": mass(lambda n,q,r,s,t: r==par.consecutive_required),
        "alive_and_allele_supported": mass(lambda n,q,r,s,t: allele_copies(n)>0),
        "alive_and_fixed_target": mass(lambda n,q,r,s,t: n[2]==par.diploids),
        "alive_and_polymorphic": mass(lambda n,q,r,s,t: 0 < allele_copies(n) < 2*par.diploids),
    }


def evolve_genotypes(initial_counts=(2,2,0), start_q="high", horizon=6,
                     par=Genetics(), toolkit=True):
    """Exact finite-multinomial propagation, unconditional on surviving.

    A fully capable founder has *no observed expression at time zero*, with
    initial run=0, seen=False. Expression requires >=1 target allele copy
    and a functional toolkit. Expression and offspring-fate are separate
    Bernoulli events conditional on the same present state; the latter
    follows the independently specified mating/seed viability mechanism.
    """
    _adult(initial_counts, par)
    if start_q not in ("high", "low"):
        raise ValueError("initial ecology must be high or low")
    if type(horizon) is not int or not 1 <= horizon <= 30:
        raise ValueError("horizon must be an integer 1..30")
    if type(toolkit) is not bool:
        raise ValueError("toolkit must be boolean")
    qi = 1 if start_q == "high" else 0
    states = {(initial_counts, qi, 0, False, False): 1.}
    extinct = extinct_seen = extinct_three = 0.
    rows = []
    for t in range(1, horizon+1):
        nxt = defaultdict(float)
        for (counts, qi, run, seen, three), weight in states.items():
            q = (par.q_low, par.q_high)[qi]
            pe = par.activation*q if toolkit and allele_copies(counts)>0 else 0.
            survival, child, _, _ = gamete_offspring_distribution(counts, qi, par)
            for expression, expr_weight in ((False,1-pe), (True,pe)):
                if expr_weight == 0:
                    continue
                seen_next = seen or expression
                dead = weight * expr_weight*(1-survival)
                extinct += dead
                if seen_next:
                    extinct_seen += dead
                if three:
                    extinct_three += dead
                surviving = weight * expr_weight * survival
                if surviving == 0:
                    continue
                new_run = min(par.consecutive_required, run+1) if expression else 0
                three_next = three or new_run == par.consecutive_required
                high_p = (par.high_stays_high if qi else par.low_becomes_high)
                high_p = min(1., high_p + (par.expression_feedback if expression else 0.))
                for child_counts,pchild in _next_counts(counts,qi,par):
                    for nqi,pq in ((0,1-high_p),(1,high_p)):
                        if pq:
                            nxt[(child_counts,nqi,new_run,seen_next,three_next)] += surviving*pchild*pq
        states = nxt
        rows.append(_report(states,extinct,extinct_seen,extinct_three,t,par))
    return {"initial_genotypes":list(initial_counts),"initial_q":start_q,"horizon":horizon,
            "toolkit":toolkit,"trajectory":rows,"final":rows[-1]}


def selfing_structure_comparison(q="high", par=Genetics(), horizon=6):
    """Matched total target-copy count k=2, but paired as two heterozygotes
    versus one homozygote. Both initial phenotypes are unexpressed.
    """
    if par.diploids != 4:
        raise ValueError("comparison requires N=4")
    heterozygote = (2,2,0)
    homozygote = (3,0,1)
    gene_k = allele_copies(heterozygote)
    assert gene_k == allele_copies(homozygote) == 2
    outputs = {}
    for name, counts in (("two_heterozygotes",heterozygote),("one_homozygote",homozygote)):
        survival, pchild, allele_freq, realized_self = gamete_offspring_distribution(counts,0 if q=="low" else 1,par)
        f = evolve_genotypes(counts,q,horizon,par)["final"]
        outputs[name] = {"first_step":{
            "survival":survival,"offspring_aa":pchild[0],"offspring_Aa":pchild[1],
            "offspring_AA":pchild[2],"offspring_target_absence":pchild[0]**par.diploids,
            "expected_target_copies":2*par.diploids*(pchild[2]+pchild[1]/2),
            "maternal_gamete_frequency":allele_freq,"realized_self_fraction":realized_self},
            "horizon":f}
    return {"initial_target_copies":gene_k,"diploids":par.diploids,"initial_q":q,
            "parameters":vars(par),"by_founder":outputs,
            "difference_absence_homo_minus_hetero":outputs["one_homozygote"]["first_step"]["offspring_target_absence"]-outputs["two_heterozygotes"]["first_step"]["offspring_target_absence"]}


def prediction_risks(par=Genetics(),q="high",horizon=6):
    """Two histories proxy genotype arrangement; both matched in visible
    phenotype, allele copy count, q, N, no ancestral genotype measurement.
    These are exact constructed Bayes-Brier risks for *first-step allele loss*
    and future alive + 3 survived expressions.
    """
    founders=((2,2,0),(3,0,1))
    cells={}
    for i,counts in enumerate(founders):
        immediate = gamete_offspring_distribution(counts, 1 if q=="high" else 0,par)[1][0]**par.diploids
        final = evolve_genotypes(counts,q,horizon,par)["final"]["alive_and_ever_three"]
        cells[i]={"allele_loss_given_successful_reproduction":immediate,"alive_and_ever_three":final}
    results={}
    for outcome in ("allele_loss_given_successful_reproduction","alive_and_ever_three"):
        a,b=(cells[0][outcome],cells[1][outcome])
        pooled=(a+b)/2
        brier_pooled=pooled*(1-pooled)
        brier_full=(a*(1-a)+b*(1-b))/2
        results[outcome]={"coarse_state_brier":brier_pooled,
                          "founder_history_brier":brier_full,
                          "full_genotype_brier":brier_full,
                          "full_genotype_plus_history_brier":brier_full,
                          "history_gain":brier_pooled-brier_full,
                          "founder_probabilities":list((a,b))}
    return results


def single_generation_loss_identity(counts, qi, par=Genetics()):
    """A sufficient *one-generation* statistic for target-allele loss.

    At the maternal pool, gamete allele probability u_g is inherited from
    parental diploid genotype g (plus one-way mutation). The weighted first
    two moments of u summarize the selfing-induced gamete correlation.

    P(offspring aa) = (1-i) [(1-u_bar)^2 + eta Var(u)]
    P(all N offspring aa | reproduction succeeds) = P(aa)^N.

    Exact within this conditional IID recruitment and all-AA immigration
    closure. It is NOT a general species-level extinction theorem.
    """
    _adult(counts,par)
    if qi not in (0,1):
        raise ValueError("ecology index must be 0 or 1")
    q=(par.q_low,par.q_high)[qi]
    direction=par.selection_strength*(q-.5)-par.direct_allele_cost*par.assurance
    w=[counts[g]*exp(.5*g*direction) for g in range(3)]
    masses=[v/sum(w) for v in w]
    u=[g/2+(1-g/2)*par.mutation for g in range(3)]
    mean=sum(a*b for a,b in zip(masses,u))
    variance=sum(a*(b-mean)**2 for a,b in zip(masses,u))
    survival,actual,_,eta=gamete_offspring_distribution(counts,qi,par)
    predicted=(1-par.immigration)*((1-mean)**2+eta*variance)
    if survival and abs(predicted-actual[0])>1e-12:
        raise AssertionError("gamete-variance identity failed")
    return {"survival":survival,"selected_gamete_mean":mean,
            "within_parent_gamete_variance":variance,
            "realized_viable_self_fraction":eta,
            "selfing_excess_non_target_zygotes":(1-par.immigration)*eta*variance,
            "offspring_aa_by_identity":predicted if survival else None,
            "offspring_aa_direct":actual[0] if survival else None,
            "conditional_complete_target_loss":predicted**par.diploids if survival else None}


def benchmark():
    """Frozen v3 cases and negative controls; all values are constructed."""
    neutral=Genetics(selection_strength=0.,direct_allele_cost=0.,inbreeding_depression=.2)
    across={}
    for assurance in (0., .5, 1.):
        p=replace(neutral,assurance=assurance)
        across[str(assurance)]=selfing_structure_comparison(par=p,horizon=6)
    p=replace(neutral,assurance=.7,selection_strength=.6)
    selection=selfing_structure_comparison(par=p,horizon=6)
    return {"protocol":"SDEA v3 constructed Mendelian selfing-and-outcrossing bridge (2026-10-26)",
            "status":"EXACT SYNTHETIC, not a natural calibration or reproduction of Model 3",
            "matched_founder_structure":across,
            "selection_perturbation":selection,
        "gamete_variance_identity":{
            str(a):{name:single_generation_loss_identity(counts,1,replace(neutral,assurance=a))
                    for name,counts in (("two_heterozygotes",(2,2,0)),
                                        ("one_homozygote",(3,0,1)))}
            for a in (0.,.5,1.)},
            "history_brier_neutral_assurance_0":prediction_risks(replace(neutral,assurance=0.)),
            "history_brier_neutral_assurance_1":prediction_risks(replace(neutral,assurance=1.)),
            "history_brier_with_selection":prediction_risks(p),
            "reopening":{
                "closed_no_target":evolve_genotypes((4,0,0),par=neutral,horizon=6)["final"],
                "one_way_mutation":evolve_genotypes((4,0,0),par=replace(neutral,mutation=.03),horizon=6)["final"],
                "immigrant_zygotes":evolve_genotypes((4,0,0),par=replace(neutral,immigration=.03),horizon=6)["final"]},
            "feedback":{
                "none":evolve_genotypes((2,2,0),par=p,horizon=6)["final"],
                "expression_conditional":evolve_genotypes((2,2,0),par=replace(p,expression_feedback=.12),horizon=6)["final"]}}
