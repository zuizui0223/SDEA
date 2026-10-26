"""SDEA v2: exact finite-state multi-generation *synthetic* accessibility model.

This is NOT Mendelian Model 3, ancestral floral restoration, or ecological
establishment. Every probability is conditional on declared transition rules.
No external datasets or empirical calibration are used.
"""
from collections import defaultdict
from dataclasses import dataclass, replace
from functools import lru_cache
from math import comb, exp, isfinite

from .benchmark import HISTORIES


def _prob(name, x):
    if not isinstance(x, (int, float)) or not isfinite(x) or not 0 <= x <= 1:
        raise ValueError(f"{name} must be finite and in [0,1]")


@dataclass(frozen=True)
class Dynamics:
    """Explicit, conditional kernel; N constant only *given* occupancy.

    Allele copies are independently resampled (binomial), not modeled as
    Mendelian offspring genotypes. Expressed phenotype requires both all
    three modules within a lineage and >=1 target allele. Expression does
    not imply fitness benefit; a distinct q/assurance survival gate applies.
    """
    diploids: int = 3
    q_low: float = .2
    q_high: float = .8
    high_stays_high: float = .85
    low_becomes_high: float = .15
    expression_recoupling: float = .0
    activation: float = .75
    assurance: float = .30
    selection_strength: float = .50
    assurance_cost: float = 1.0
    mutation: float = .0
    immigrant_fraction: float = .0
    immigrant_allele_frequency: float = 1.0
    module_retention: float = 1.0
    module_reacquisition: float = .0
    post_expression_damage: float = .0
    consecutive_required: int = 3

    def __post_init__(self):
        if type(self.diploids) is not int or not 1 <= self.diploids <= 16:
            raise ValueError("diploids must be an integer in [1,16]")
        if type(self.consecutive_required) is not int or not 2 <= self.consecutive_required <= 12:
            raise ValueError("consecutive_required must be an integer in [2,12]")
        for name in ("q_low", "q_high", "high_stays_high", "low_becomes_high",
                     "expression_recoupling", "activation", "assurance", "mutation",
                     "immigrant_fraction", "immigrant_allele_frequency", "module_retention",
                     "module_reacquisition", "post_expression_damage"):
            _prob(name, getattr(self, name))
        if self.q_low > self.q_high:
            raise ValueError("q_low cannot exceed q_high")
        for name in ("selection_strength", "assurance_cost"):
            value = getattr(self, name)
            if not isfinite(value) or value < 0 or value > 200:
                raise ValueError(f"{name} must be finite and between 0 and 200")


def history_masks(label):
    if label not in HISTORIES:
        raise ValueError("unknown parity history")
    return {mask: 1 / len(HISTORIES[label]) for mask in HISTORIES[label]}


def _binom_probs(n, p):
    if p == 0:
        return ((0, 1.),)
    if p == 1:
        return ((n, 1.),)
    return tuple((k, comb(n, k) * p**k * (1-p)**(n-k)) for k in range(n+1))


@lru_cache(maxsize=20000)
def _transition(mask, copies, qi, expressed, par):
    """Next state kernel conditional on this generation surviving.

    Order: visible expression -> survival -> module/allele transmission ->
    next ecological q. Post-expression damage destroys the *whole* toolkit,
    not a molecularly observed pathway; optional input can reopen pathways.
    """
    q = (par.q_low, par.q_high)[qi]
    p = copies / (2 * par.diploids)
    w = exp(par.selection_strength * (q-.5) - par.assurance_cost * par.assurance)
    ps = p*w / (1-p+p*w)
    ps = ps + (1-ps)*par.mutation
    ps = (1-par.immigrant_fraction)*ps + par.immigrant_fraction*par.immigrant_allele_frequency
    allele = _binom_probs(2 * par.diploids, ps)
    transition_high = (par.high_stays_high if qi else par.low_becomes_high)
    transition_high = min(1.0, transition_high + par.expression_recoupling * expressed)
    ecologies = ((0, 1-transition_high), (1, transition_high))
    bit_probs = tuple(par.module_retention if mask & (1<<bit) else par.module_reacquisition
                      for bit in range(3))
    masks = []
    for newmask in range(8):
        prob = 1.
        for bit in range(3):
            t = bit_probs[bit]
            prob *= t if newmask & (1<<bit) else 1-t
        if prob:
            prob *= (1-par.post_expression_damage) if expressed else 1.
            if prob:
                masks.append((newmask, prob))
    if expressed and par.post_expression_damage:
        masks.append((0, par.post_expression_damage))
    result = defaultdict(float)
    for newmask, pm in masks:
        for newcopies, pa in allele:
            for newqi, pq in ecologies:
                if pm and pa and pq:
                    result[(newmask, newcopies, newqi)] += pm * pa * pq
    return tuple(result.items())


def _snapshot(states, extinct, extinct_seen, extinct_three, generation, par):
    alive = sum(states.values())
    def total(condition):
        return sum(prob for (m,k,qi,streak,seen,three),prob in states.items()
                   if condition(m,k,qi,streak,seen,three))
    return {
        "generation": generation,
        "alive": alive,
        "extinct": extinct,
        "first_expression_ever": extinct_seen + total(lambda m,k,q,r,s,t: s),
        "three_consecutive_survived_ever": extinct_three + total(lambda m,k,q,r,s,t: t),
        "alive_and_ever_three": total(lambda m,k,q,r,s,t: t),
        "alive_and_current_three": total(lambda m,k,q,r,s,t: r == par.consecutive_required),
        "alive_and_allele_supported": total(lambda m,k,q,r,s,t: k > 0),
        "alive_and_joint_toolkit_supported": total(lambda m,k,q,r,s,t: m == 7 and k > 0),
    }


def evolve_exact(mask_prior, initial_copies=2, start_q="high", horizon=8, par=Dynamics()):
    """Enumerate all finite-state descendants; report *unconditional* risks.

    Initial observation is phenotype-absent even for complete latent toolkits.
    First-expression events count even if the lineage dies in that generation;
    a three-consecutive event requires expression AND survival in each of the
    three consecutive reproduction updates. It is NOT indefinite stability.
    """
    if start_q not in ("low", "high"):
        raise ValueError("start_q must be low or high")
    if type(initial_copies) is not int or not 0 <= initial_copies <= 2*par.diploids:
        raise ValueError("initial copies outside available allele slots")
    if type(horizon) is not int or not 1 <= horizon <= 50:
        raise ValueError("horizon must be integer 1..50")
    if not mask_prior or any(type(m) is not int or not 0 <= m <= 7 for m in mask_prior):
        raise ValueError("invalid mask prior")
    for mass in mask_prior.values():
        _prob("mask prior weight", mass)
    if abs(sum(mask_prior.values())-1.) > 1e-12:
        raise ValueError("mask prior must sum to one")
    qi = 1 if start_q == "high" else 0
    states = defaultdict(float, {(m, initial_copies, qi, 0, False, False):mass
                                 for m,mass in mask_prior.items() if mass})
    extinct = extinct_seen = extinct_three = 0.
    rows = []
    for generation in range(1, horizon + 1):
        next_states = defaultdict(float)
        for (mask, copies, ecological_index, run, seen, three), mass in states.items():
            q = (par.q_low, par.q_high)[ecological_index]
            expr_p = par.activation*q if mask == 7 and copies > 0 else 0.
            survival = q + (1-q)*par.assurance
            for expression, pe in ((False,1-expr_p),(True,expr_p)):
                if pe == 0:
                    continue
                reached_seen = seen or expression
                extinct_here = mass * pe * (1-survival)
                extinct += extinct_here
                if reached_seen:
                    extinct_seen += extinct_here
                if three:
                    extinct_three += extinct_here
                live_mass = mass * pe * survival
                if live_mass == 0:
                    continue
                new_run = min(par.consecutive_required, run+1) if expression else 0
                reached_three = three or new_run == par.consecutive_required
                for (nm,nc,nq),pt in _transition(mask,copies,ecological_index,expression,par):
                    next_states[(nm,nc,nq,new_run,reached_seen,reached_three)] += live_mass*pt
        states = next_states
        rows.append(_snapshot(states,extinct,extinct_seen,extinct_three,generation,par))
    return {"horizon": horizon, "initial_copies":initial_copies,"initial_q":start_q,
            "trajectory":rows,"final":rows[-1]}


def patch_alignment(horizon=8, par=Dynamics(), initial_copies=2):
    """Independent patches only; same coarse marginals, reversed mask/q map."""
    def patch(m,q):
        return evolve_exact({m:1.}, initial_copies, q, horizon, par)["final"]
    ah, zl = patch(7,"high"), patch(0,"low")
    al, zh = patch(7,"low"), patch(0,"high")
    def union(a,b,key):
        return 1-(1-a[key])*(1-b[key])
    key="three_consecutive_survived_ever"
    return {"aligned":union(ah,zl,key),"misaligned":union(al,zh,key),
            "difference":union(ah,zl,key)-union(al,zh,key),
            "endpoint":key,"independent_patches":True}


def run_multigeneration():
    """Frozen v2 baseline, controls and exploratory tradeoff; no data fitting."""
    p = Dynamics()
    empty = evolve_exact({0:1.},horizon=8,par=p)["final"]
    complete = evolve_exact({7:1.},horizon=8,par=p)["final"]
    history = {label:evolve_exact(history_masks(label),horizon=8,par=p)["final"]
               for label in ("even","odd")}
    closed = evolve_exact({7:1.},initial_copies=0,horizon=8,par=p)["final"]
    reopened_allele = evolve_exact({7:1.},initial_copies=0,horizon=8,
                                    par=replace(p,mutation=.05))["final"]
    reopened_modules = evolve_exact({0:1.},horizon=8,
                                     par=replace(p,module_reacquisition=.12))["final"]
    damage = evolve_exact({7:1.},horizon=8,
                           par=replace(p,post_expression_damage=1.))["final"]
    costs = {}
    for c in (0.,3.):
        for a in (.1,.9):
            pp = replace(p, assurance=a, assurance_cost=c)
            costs[f"c{c}_a{a}"] = evolve_exact({7:1.},horizon=8,par=pp)["final"]
    return {
        "protocol":"SDEA multigeneration synthetic exact v2 (2026-10-26)",
        "source":"new constructed finite transition kernel; not reexecution of ANTO, EGWE, izu-core",
        "parameters":vars(p),
        "baseline":{"complete":complete,"incomplete":empty,"histories":history},
        "alignment":patch_alignment(horizon=8,par=p),
        "negative_controls":{"allele_absent_closed":closed,"allele_mutation_input":reopened_allele,
                             "modules_absent_closed":empty,"module_reacquisition_input":reopened_modules,
                             "complete_with_certain_post_expression_damage":damage},
        "assurance_cost_conditions":costs,
        "time_course_complete":evolve_exact({7:1.},horizon=12,par=p)["trajectory"],
        "prediction_sufficiency":horizon_brier_risks(horizon=8,par=p),
        "prediction_null_history":horizon_brier_risks(horizon=8,par=p,null_history=True),
        "finite_horizon_frontiers":[multigen_cost_crossover(endpoint=e,horizon=8,par=p)
             for e in ("alive_and_ever_three","three_consecutive_survived_ever",
                       "alive_and_joint_toolkit_supported")],
    }


def horizon_brier_risks(horizon=8, par=Dynamics(), initial_copies=2,
                        sensitivity=.9, specificity=.95, null_history=False):
    """Exact Bayes risk for a future *alive-and-ever-3* event; no fitted data.

    Initial observed phenotype is negative in both parity histories. The
    hidden full biological state is (mask, copies, q); copies and q are held
    fixed in this comparison. Assays have declared independent errors.
    """
    from .benchmark import assay_likelihood
    _prob("sensitivity", sensitivity)
    _prob("specificity", specificity)
    out = {label:defaultdict(lambda:[0.,0.]) for label in
           ("visible", "history", "assay", "assay_history", "full", "full_history")}
    futures = {mask:evolve_exact({mask:1.},initial_copies,"high",horizon,par)
               ["final"]["alive_and_ever_three"] for mask in range(8)}
    for history in HISTORIES:
        masks = range(8) if null_history else HISTORIES[history]
        for mask in masks:
            prior=.5 / len(masks)
            for observed in range(8):
                mass=prior*assay_likelihood(observed,mask,sensitivity,specificity)
                keys={"visible":(),"history":(history,),"assay":(observed,),
                      "assay_history":(history,observed),"full":(mask,),
                      "full_history":(mask,history)}
                for scheme,key in keys.items():
                    bucket=out[scheme][key]
                    bucket[0]+=mass
                    bucket[1]+=mass*futures[mask]
    return {name:sum(s-s*s/w for w,s in groups.values() if w)
            for name,groups in out.items()}


def multigen_cost_crossover(endpoint="alive_and_ever_three",horizon=8,
                             par=Dynamics(), low=.1, high=.9,
                             bracket_upper=8., iterations=48):
    """Exploratory *finite-horizon* crossover, not a universal critical cost.

    Sign brackets one numerical root; existence/uniqueness beyond that range
    are not assumed. The endpoint is explicitly chosen before bisection.
    """
    allowed=("alive_and_ever_three", "three_consecutive_survived_ever",
             "alive_and_joint_toolkit_supported")
    if endpoint not in allowed:
        raise ValueError("endpoint cannot be used to define the crossover")
    if not (0 <= low < high <= 1):
        raise ValueError("assurance policies must be ordered in [0,1]")
    if not isfinite(bracket_upper) or bracket_upper <= 0:
        raise ValueError("invalid bracket_upper")
    if type(iterations) is not int or iterations < 1:
        raise ValueError("iterations must be positive integer")
    if par.mutation or par.immigrant_fraction:
        raise ValueError("frontier only for closed allele input")
    def difference(cost):
        pf=replace(par,assurance_cost=cost)
        a=evolve_exact({7:1.}, horizon=horizon,par=replace(pf,assurance=low))["final"][endpoint]
        b=evolve_exact({7:1.}, horizon=horizon,par=replace(pf,assurance=high))["final"][endpoint]
        return b-a
    v0, vupper=difference(0.), difference(bracket_upper)
    if not v0>0 or not vupper<0:
        raise ValueError("crossover not bracketed; do not invent one")
    lo,hi=0.,bracket_upper
    for _ in range(iterations):
        mid=(lo+hi)/2.
        if difference(mid)>0:
            lo=mid
        else:
            hi=mid
    root=(lo+hi)/2
    return {"endpoint":endpoint,"horizon":horizon,"cost_star":root,
            "diff_at_cost_zero":v0,"diff_at_upper":vupper,"root_residual":difference(root)}
