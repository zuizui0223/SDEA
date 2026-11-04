"""SDEA v4: exact *constructed* two-patch reproductive route/flow comparison.

A small recipient population receives pollen or whole zygote genotypes from a
STATIC donor population. Genetic routes A and B are separate, unlinked, binary
loci with Mendelian segregation and independent assortment. There is no
explicit ancestral genealogy, sampled floral pigment, or molecular reversion.

All probabilities are *conditional on local viable reproduction* and a fixed
number of recruits, and the donor remains a stipulated source reservoir.
This is not izu-core Model 3 and is not a coupled landscape metapopulation.
"""
from collections import defaultdict
from dataclasses import dataclass, replace
from functools import lru_cache
from itertools import product
from math import comb, isfinite

GENOTYPES = tuple(product(range(3), repeat=2))  # (copies_ancestral, copies_alternative)
INDEX = {g: i for i, g in enumerate(GENOTYPES)}
ROUTES = ("none", "ancestral", "alternative", "simultaneous")
ORIGINS = ("selfed_local", "outcross_local", "immigrant_pollen", "immigrant_seed")


def _unit(name, value):
    if type(value) not in (int, float) or not isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{name} must be finite in [0,1]")


@dataclass(frozen=True)
class Flow:
    recruits: int = 2
    recipient_q: float = .6
    donor_q: float = .8
    assurance: float = .25
    inbreeding_depression: float = .15
    pollen_rate: float = .0
    seed_rate: float = .0
    activation_ancestral: float = .8
    activation_alternative: float = .8
    tie_ancestral: float = .5

    def __post_init__(self):
        if type(self.recruits) is not int or not 1 <= self.recruits <= 5:
            raise ValueError("recruits must be an integer from 1 to 5")
        for name in ("recipient_q", "donor_q", "assurance", "inbreeding_depression",
                     "pollen_rate", "seed_rate", "activation_ancestral",
                     "activation_alternative", "tie_ancestral"):
            _unit(name, getattr(self, name))
        for q in (self.recipient_q, self.donor_q):
            if q == 0 and (self.assurance == 0 or self.inbreeding_depression == 1):
                raise ValueError("no viable locally produced zygotes; source/recipient conditional kernel undefined")


def counts(*genotypes, n=2):
    """Immutable adult counts given a list of 2-locus diploid genotypes."""
    if len(genotypes) != n:
        raise ValueError("need exactly n diploid adults")
    result = [0] * len(GENOTYPES)
    for g in genotypes:
        if g not in INDEX:
            raise ValueError("invalid diploid genotype")
        result[INDEX[g]] += 1
    return tuple(result)


def _counts_check(state, n):
    if not isinstance(state, tuple) or len(state) != 9 or any(type(x) is not int or x < 0 for x in state):
        raise ValueError("state must be a tuple of nine nonnegative genotype counts")
    if sum(state) != n:
        raise ValueError("genotype counts must sum to recruits")


def route_copies(state, route):
    if route not in (0, 1):
        raise ValueError("route index must be 0 or 1")
    return sum(g[route] * state[i] for i, g in enumerate(GENOTYPES))


def _mix_add(target, weight, distribution):
    for j, p in enumerate(distribution):
        target[j] += weight * p


@lru_cache(maxsize=None)
def cross(maternal, paternal):
    """Independent meioses and independent assortment at 2 unlinked loci."""
    if maternal not in INDEX or paternal not in INDEX:
        raise ValueError("unknown parental genotype")
    probabilities = []
    for locus in range(2):
        p, q = maternal[locus] / 2, paternal[locus] / 2
        probabilities.append(((1-p)*(1-q), p*(1-q)+(1-p)*q, p*q))
    return tuple(probabilities[0][g[0]] * probabilities[1][g[1]] for g in GENOTYPES)


def reproduction_weights(q, assurance, depression):
    """V3-compatible viable selfing and outcrossing weights; conditional on viability."""
    s = assurance * (1-q)
    out = (1-s) * q
    selfed = s * (1-depression)
    total = selfed+out
    if total <= 0:
        raise ValueError("no viable local reproduction")
    return {"viability_gate": total, "selfed": selfed/total, "outcross": out/total}


def _adult_probs(state):
    n = sum(state)
    if n == 0:
        raise ValueError("empty donor or recipient not supported")
    return tuple(x/n for x in state)


def _zygote(state, q, par):
    """Locally produced viable selfed and outcrossed zygotes; source seeds use these."""
    moms = _adult_probs(state)
    w = reproduction_weights(q, par.assurance, par.inbreeding_depression)
    own = [0.] * 9
    crossed = [0.] * 9
    for i, mi in enumerate(moms):
        if mi == 0:
            continue
        _mix_add(own, mi, cross(GENOTYPES[i], GENOTYPES[i]))
        for j, fj in enumerate(moms):
            if fj:
                _mix_add(crossed, mi*fj, cross(GENOTYPES[i],GENOTYPES[j]))
    return (tuple(own), tuple(crossed), w)


@lru_cache(maxsize=30000)
def tagged_offspring(recipient, donor, par=Flow()):
    """Per recruited offspring: P(genotype, parental/propagule origin).

    Pollen migration replaces PATERNAL gamete source in viable outcrossing,
    never the mother and never selfed gametes. Seed migration replaces an
    entire offspring genotype with a zygote from donor's *local reproduction*.
    This is a conditional fixed-recruit mixture: migration fractions are not
    rates per km, and an unoccupied recipient cannot be recolonized here.
    """
    _counts_check(recipient, par.recruits)
    _counts_check(donor, par.recruits)
    focal_p = _adult_probs(recipient)
    donor_p = _adult_probs(donor)
    own, local, w = _zygote(recipient, par.recipient_q, par)
    donor_own, donor_out, donor_w = _zygote(donor, par.donor_q, par)
    donor_seed = tuple(donor_w["selfed"]*donor_own[g] + donor_w["outcross"]*donor_out[g] for g in range(9))
    migrant_pollen = [0.] * 9
    for i, m in enumerate(focal_p):
        for j, d in enumerate(donor_p):
            if m and d:
                _mix_add(migrant_pollen, m*d, cross(GENOTYPES[i],GENOTYPES[j]))
    wlocal = 1-par.seed_rate
    source = {
        "selfed_local": tuple(wlocal*w["selfed"]*x for x in own),
        "outcross_local": tuple(wlocal*w["outcross"]*(1-par.pollen_rate)*x for x in local),
        "immigrant_pollen": tuple(wlocal*w["outcross"]*par.pollen_rate*x for x in migrant_pollen),
        "immigrant_seed": tuple(par.seed_rate*x for x in donor_seed),
    }
    result = tuple(sum(source[k][j] for k in ORIGINS) for j in range(9))
    if abs(sum(result) - 1) > 1e-12 or min(result) < -1e-14:
        raise AssertionError("source-tagged offspring not normalized")
    return {"origin": source, "genotype": result,
            "local_viability_gate": w["viability_gate"],
            "effective_local_selfing": w["selfed"],
            "donor_viability_gate": donor_w["viability_gate"]}


@lru_cache(maxsize=None)
def all_count_states(n):
    if type(n) is not int or n < 1 or n > 5:
        raise ValueError("n must be an integer 1..5")
    def build(rem, parts):
        if len(parts) == 8:
            yield tuple(parts + [rem])
        else:
            for i in range(rem + 1):
                yield from build(rem-i, parts + [i])
    return tuple(build(n, []))


@lru_cache(maxsize=30000)
def _recruits(recipient, donor, par):
    p = tagged_offspring(recipient,donor,par)["genotype"]
    out=[]
    for state in all_count_states(par.recruits):
        prob=1.
        remaining=par.recruits
        for n, pp in zip(state, p):
            prob *= comb(remaining,n)*(pp**n)
            remaining-=n
        if prob:
            out.append((state,prob))
    if abs(sum(w for _,w in out)-1)>2e-12:
        raise AssertionError("finite multinomial recruitment not normalized")
    return tuple(out)


def route_probabilities(genotype, q, par=Flow()):
    """Exclusive *individual* expression routes, with a declared tie rule."""
    if genotype not in INDEX:
        raise ValueError("unknown genotype")
    _unit("q",q)
    a = par.activation_ancestral*q if genotype[0] else 0.
    b = par.activation_alternative*(1-q) if genotype[1] else 0.
    both=a*b
    return ( (1-a)*(1-b),
             a*(1-b) + par.tie_ancestral*both,
             b*(1-a) + (1-par.tie_ancestral)*both )


@lru_cache(maxsize=20000)
def population_expression(state, q, par=Flow()):
    """Same-generation labels none, A-only, B-only, A-and-B simultaneously.

    Independent flowering expression among conditional recruited individuals;
    earliest observed generation can contain both types, so forcing a single
    winning route would overinterpret contemporaneous flowers.
    """
    _counts_check(state,par.recruits)
    _unit("q",q)
    none=no_a=no_b=1.
    for i, n in enumerate(state):
        if n:
            p0, pa, pb = route_probabilities(GENOTYPES[i],q,par)
            none *= p0**n
            no_a *= (1-pa)**n
            no_b *= (1-pb)**n
    vals=(none, no_b-none, no_a-none, 1-no_a-no_b+none)
    if min(vals) < -1e-12 or abs(sum(vals)-1)>2e-12:
        raise AssertionError("expression categories not normalized")
    return tuple(max(0.,x) for x in vals)


def one_generation(recipient, donor, par=Flow()):
    """Exact one-generation genetic accessibility, entry provenance and visible routes."""
    tags=tagged_offspring(recipient,donor,par)
    n=par.recruits
    offspring=tags["genotype"]
    route_support={key: 1-(1-sum(offspring[i] for i,g in enumerate(GENOTYPES)
                                     if g[j]>0))**n
                   for j,key in ((0,"ancestral"),(1,"alternative"))}
    gen_joint=sum(offspring[i] for i,g in enumerate(GENOTYPES) if g[0] and g[1])
    first=[0.]*4
    for new,w in _recruits(recipient,donor,par):
        for i,x in enumerate(population_expression(new,par.recipient_q,par)):
            first[i]+=w*x
    mean_copies=[sum(offspring[i]*g[j] for i,g in enumerate(GENOTYPES))*n for j in (0,1)]
    return {"route_support_any":route_support,
            "at_least_one_dual_competent":1-(1-gen_joint)**n,
            "first_generation_expression":dict(zip(ROUTES,first)),
            "mean_allele_copies":dict(zip(("ancestral","alternative"),mean_copies)),
            "effective_local_selfing":tags["effective_local_selfing"],
            "local_viability_gate":tags["local_viability_gate"],
            "genotype_distribution":{f"{g[0]}-{g[1]}":offspring[i] for i,g in enumerate(GENOTYPES)}}


def new_route_source(recipient, donor, route=1, par=Flow()):
    """Exact disjoint provenance when recipient initially has *zero* route copies.

    Both-source is possible in the same cohort but never for a single recruit.
    This records source of a genetic route, NOT observed expressed ancestry.
    """
    _counts_check(recipient,par.recruits)
    _counts_check(donor,par.recruits)
    if route not in (0,1) or route_copies(recipient,route):
        raise ValueError("recipient must initially lack chosen route")
    tag=tagged_offspring(recipient,donor,par)["origin"]
    p=sum(tag["immigrant_pollen"][i] for i,g in enumerate(GENOTYPES) if g[route]>0)
    s=sum(tag["immigrant_seed"][i] for i,g in enumerate(GENOTYPES) if g[route]>0)
    n=par.recruits
    none=(1-p-s)**n
    pollen_only=(1-s)**n-none
    seed_only=(1-p)**n-none
    both=1-none-pollen_only-seed_only
    if abs(none+pollen_only+seed_only+both-1)>1e-12:
        raise AssertionError("provenance categories do not partition input")
    return {"route":("ancestral","alternative")[route],
            "per_recruit_from_pollen":p,"per_recruit_from_seed":s,
            "probabilities":{"none":none,"pollen_only":pollen_only,
                             "seed_only":seed_only,"both":both},
            "any_route_reopened":1-none}


def evolve_recipient(recipient, donor, par=Flow(), horizon=4, pulse=False,
                     q_schedule=None):
    """Exact receiver DP, first observed A/B route, static donor reservoir.

    Genotype counts of n recruits evolve as a finite multinomial. The source is
    fixed over time; no return migration or source demographic fluctuations.
    When pulse=True, both flows occur at t=1 only, then stop. Expression is
    generated in recruited offspring at the end of each generation.
    """
    _counts_check(recipient,par.recruits)
    _counts_check(donor,par.recruits)
    if type(horizon) is not int or not 1<=horizon<=12:
        raise ValueError("horizon must be integer 1..12")
    if type(pulse) is not bool:
        raise ValueError("pulse must be bool")
    if q_schedule is not None:
        if not isinstance(q_schedule,(tuple,list)) or len(q_schedule)!=horizon:
            raise ValueError("q schedule must have exactly horizon entries")
        for q in q_schedule:
            _unit("future forcing q",q)
    states={(recipient,0):1.}
    trajectory=[]
    for step in range(1,horizon+1):
        pp=par if not pulse or step==1 else replace(par,pollen_rate=0.,seed_rate=0.)
        q=par.recipient_q if q_schedule is None else q_schedule[step-1]
        pp=replace(pp,recipient_q=q)
        new=defaultdict(float)
        for (adults, label), mass in states.items():
            for recruited, count_prob in _recruits(adults,donor,pp):
                p=mass*count_prob
                if label:
                    new[(recruited,label)]+=p
                else:
                    for l,prob in enumerate(population_expression(recruited,q,pp)):
                        if prob:
                            new[(recruited,l)]+=p*prob
        states=new
        cats=[0.]*4
        a=b=dual=0.
        for (adult,label),p in states.items():
            cats[label]+=p
            if route_copies(adult,0):a+=p
            if route_copies(adult,1):b+=p
            if any(n and g[0] and g[1] for n,g in zip(adult,GENOTYPES)):dual+=p
        row={"generation":step,"first_route_by_horizon":dict(zip(ROUTES,cats)),
             "allele_supported_ancestral":a,"allele_supported_alternative":b,
             "joint_in_individual":dual,"mass":sum(cats)}
        if abs(row["mass"]-1)>2e-11:
            raise AssertionError("DP state mass is not conserved")
        trajectory.append(row)
    return {"horizon":horizon,"pulse":pulse,"q_schedule":q_schedule,
            "trajectory":trajectory,"final":trajectory[-1]}



def history_identifiability(recipient, donor, pollen_flow, seed_flow):
    """Can one infer *the movement process* from one post-arrival census?

    For two equally probable interventions, binary presence of any visible
    phenotype may have the same distribution even though full diploid-genotype
    counts do not. This is an exact total-variation classifier diagnostic,
    NOT a natural historical-ancestry assay or a test of historical causation.
    """
    _counts_check(recipient,pollen_flow.recruits)
    _counts_check(donor,pollen_flow.recruits)
    from dataclasses import fields
    for field in fields(Flow):
        if field.name not in ("pollen_rate","seed_rate") and getattr(pollen_flow,field.name)!=getattr(seed_flow,field.name):
            raise ValueError("observer comparison requires identical non-migration parameters")
    o1=one_generation(recipient,donor,pollen_flow)
    o2=one_generation(recipient,donor,seed_flow)
    appear1=1-o1["first_generation_expression"]["none"]
    appear2=1-o2["first_generation_expression"]["none"]
    p1=dict(_recruits(recipient,donor,pollen_flow))
    p2=dict(_recruits(recipient,donor,seed_flow))
    tv=.5*sum(abs(p1.get(state,0)-p2.get(state,0)) for state in all_count_states(pollen_flow.recruits))
    return {"visible_any_expression_pollen":appear1,
            "visible_any_expression_seed":appear2,
            "visible_binary_total_variation":abs(appear1-appear2),
            "full_genotype_count_total_variation":tv,
            "best_equal_prior_origin_accuracy_from_visible":.5+.5*abs(appear1-appear2),
            "best_equal_prior_origin_accuracy_from_full_genotypes":.5+.5*tv,
            "note":"classification identifies constructed movement assignment, not natural ancestral regain"}


def benchmark():
    """Fixed, auditable synthetic comparisons (no natural parameter tuning)."""
    none=(0,0); aa=(2,0); bb=(0,2)
    focal_empty=counts(none,none)
    donor_ancestor=counts(aa,aa)
    source_alternative=counts(bb,bb)
    focal_ancestor=counts(aa,aa)
    donor_mixed=counts(aa,bb)
    noflow=Flow(recruits=2,recipient_q=.6,assurance=.35)
    pollen=replace(noflow,pollen_rate=.4)
    # Matched immediate entrant-carrier fraction (not equal gene-copy flow):
    # donor is AA and recipient is aa, hence pollen carrier = outcross_fraction*pollen_rate
    # whereas a whole donor seed carries the genotype with probability seed_rate.
    matched_seed_rate=.4*reproduction_weights(noflow.recipient_q,noflow.assurance,
                                              noflow.inbreeding_depression)["outcross"]
    seed=replace(noflow,seed_rate=matched_seed_rate)
    seed_nominal=replace(noflow,seed_rate=.4)
    both=replace(noflow,pollen_rate=.4,seed_rate=.4)
    restored={}
    for label,par in (("closed",noflow),("pollen",pollen),("seed_matched_carriers",seed),
                      ("seed_equal_nominal_rate",seed_nominal)):
        restored[label]={"immediate":one_generation(focal_empty,donor_ancestor,par),
                         "horizon_4_pulse":evolve_recipient(focal_empty,donor_ancestor,par,4,True)["final"],
                         "route_provenance":new_route_source(focal_empty,donor_ancestor,0,par)}
    hybrid={}
    for label,par in (("closed",noflow),("pollen",pollen),("seed",seed_nominal),("both",both)):
        hybrid[label]={"immediate":one_generation(focal_ancestor,source_alternative,par),
                       "horizon_4_persistent":evolve_recipient(focal_ancestor,source_alternative,par,4,False)["final"],
                       "alternative_provenance":new_route_source(focal_ancestor,source_alternative,1,par)}
    regimens={}
    for q in (.2,.5,.8):
        par=replace(noflow,recipient_q=q,pollen_rate=.4,seed_rate=.4)
        regimens[str(q)]=one_generation(focal_empty,donor_mixed,par)
    neutral=replace(noflow,recipient_q=.5,pollen_rate=.4,seed_rate=.4)
    q_shift=replace(noflow,recipient_q=.2,pollen_rate=.4,seed_rate=.4)
    return {"protocol":"SDEA v4 exact source-recipient two-route flow (2026-10-26)",
            "status":"SYNTHETIC; independent unlinked loci, fixed donor, conditional recruitment, not re-evolution",
            "defaults":vars(noflow),
            "matched_seed_rate_for_carrier_equivalence":matched_seed_rate,
            "phenotype_equivalent_but_genetically_distinct_history":history_identifiability(focal_empty,donor_ancestor,pollen,seed),
            "ancestral_reopening_pulse":restored,
            "alternative_invasion_and_joint_competence":hybrid,
            "symmetrical_route_donor_q_gradient":regimens,
            "neutral_environment_symmetry":one_generation(focal_empty,donor_mixed,neutral),
            "environment_switch_alternative":evolve_recipient(focal_empty,donor_mixed,q_shift,4,True,(.2,.2,.8,.8))["final"],
            "constant_low_environment":evolve_recipient(focal_empty,donor_mixed,q_shift,4,True)["final"]}
