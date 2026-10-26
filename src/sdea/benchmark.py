"""Exact counterexample benchmarks for state-dependent evolutionary accessibility.

All mechanisms and numerical parameters are constructed for this benchmark.
These are NOT fitted to, or independent reproductions of, ANTO/EGWE/izu-core.
No natural trait-regain, individual fitness, or demographic rate is estimated.
"""

from collections import defaultdict
from dataclasses import dataclass, replace
from math import exp, isfinite

HISTORIES = {"even": (0, 3, 5, 6), "odd": (1, 2, 4, 7)}
OBSERVATIONS = ("visible", "history", "assay", "assay_history", "full", "full_history")
ENDPOINTS = ("expression", "survival", "genetic_support", "joint")


def _unit(name, value):
    if not isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{name} must be finite and in [0, 1]")


@dataclass(frozen=True)
class Parameters:
    """One-step model; allele copies are independent conditional binomial draws.

    `assurance_cost` is an explicitly assumed penalty on a target allele's
    relative transmission weight, not an estimate of inbreeding depression.
    """

    activation: float = 0.65
    assurance: float = 0.30
    allele_frequency: float = 0.12
    offspring_diploids: int = 4
    selection_strength: float = 0.50
    assurance_cost: float = 1.00
    mutation_rate: float = 0.0
    immigrant_fraction: float = 0.0
    immigrant_allele_frequency: float = 1.0

    def __post_init__(self):
        for name in ("activation", "assurance", "allele_frequency", "mutation_rate",
                     "immigrant_fraction", "immigrant_allele_frequency"):
            _unit(name, getattr(self, name))
        for name in ("selection_strength", "assurance_cost"):
            val = getattr(self, name)
            if not isfinite(val) or val < 0:
                raise ValueError(f"{name} must be nonnegative and finite")
        if type(self.offspring_diploids) is not int or self.offspring_diploids < 1:
            raise ValueError("offspring_diploids must be a positive integer")


def _validate(mask, q):
    if type(mask) is not int or not 0 <= mask <= 7:
        raise ValueError("module mask must be an integer from 0 to 7")
    _unit("ecological interaction q", q)


def transmitted_frequency(q, par):
    """Target allele's modeled frequency after selection, mutation and arrivals.

    w = exp(strength*(q-0.5) - cost*assurance); p' = pw/(1-p+pw).
    Mutation is a one-way 0->1 event; migration replaces a declared fraction.
    """
    _unit("ecological interaction q", q)
    w = exp(par.selection_strength * (q - 0.5) - par.assurance_cost * par.assurance)
    p = par.allele_frequency
    selected = p * w / (1 - p + p * w)
    mutated = selected + (1 - selected) * par.mutation_rate
    return (1 - par.immigrant_fraction) * mutated + par.immigrant_fraction * par.immigrant_allele_frequency


def outcome_probability(mask, q, par=Parameters(), endpoint="joint"):
    """Exact one-step probability for a chosen endpoint under *independent gates*.

    expression: first visible expression, given all 3 modules co-occur.
    survival: a constructed one-step occupancy/reproductive gate.
    genetic_support: at least one target allele in 2N binomial draws.
    joint: all three gates; NOT indefinite establishment or actual fixation.
    """
    _validate(mask, q)
    if endpoint not in ENDPOINTS:
        raise ValueError("unsupported endpoint")
    expression = par.activation * q if mask == 7 else 0.0
    survival = q + (1 - q) * par.assurance
    p = transmitted_frequency(q, par)
    support = 1 - (1 - p) ** (2 * par.offspring_diploids)
    return {
        "expression": expression,
        "survival": survival,
        "genetic_support": support,
        "joint": expression * survival * support,
    }[endpoint]


def module_marginals(history):
    """First, pairwise and joint retained-toolkit probabilities."""
    masks = HISTORIES[history]
    first = [sum(bool(m & (1 << i)) for m in masks) / len(masks) for i in range(3)]
    pair = [sum(bool(m & (1 << i)) and bool(m & (1 << j)) for m in masks) / len(masks)
            for i, j in ((0, 1), (0, 2), (1, 2))]
    return {"individual": first, "pairs": pair, "all_three": masks.count(7) / len(masks)}


def assay_likelihood(observed_mask, true_mask, sensitivity, specificity):
    """Independent binary assays, sensitivity/specificity specified rather than measured."""
    if not (0 <= observed_mask <= 7 and 0 <= true_mask <= 7):
        raise ValueError("invalid assay mask")
    _unit("assay sensitivity", sensitivity)
    _unit("assay specificity", specificity)
    likelihood = 1.0
    for bit in range(3):
        truth = bool(true_mask & (1 << bit))
        observed = bool(observed_mask & (1 << bit))
        positive = sensitivity if truth else 1 - specificity
        likelihood *= positive if observed else (1 - positive)
    return likelihood


def brier_risks(par=Parameters(), q=0.65, sensitivity=0.9, specificity=0.95,
                null_history=False):
    """Exact Bayes-optimal Brier risk under each prospective observation scheme.

    Population: equally likely histories; each history has a specified mask law.
    White visible state is constant. Event is the next-step `joint` endpoint.
    `null_history`: the two labels have identical uniform toolkit distributions.
    The future forcing q and kernel par are held fixed across both histories.
    """
    _unit("q", q)
    _unit("sensitivity", sensitivity)
    _unit("specificity", specificity)
    buckets = {name: defaultdict(lambda: [0., 0.]) for name in OBSERVATIONS}
    for history in HISTORIES:
        masks = tuple(range(8)) if null_history else HISTORIES[history]
        prior = 0.5 / len(masks)
        for mask in masks:
            event = outcome_probability(mask, q, par, "joint")
            for observed in range(8):
                mass = prior * assay_likelihood(observed, mask, sensitivity, specificity)
                keys = {
                    "visible": (), "history": (history,), "assay": (observed,),
                    "assay_history": (observed, history), "full": (mask,),
                    "full_history": (mask, history),
                }
                for scheme, key in keys.items():
                    bucket = buckets[scheme][key]
                    bucket[0] += mass
                    bucket[1] += mass * event
    # For Y~Bernoulli, E[(Y - E[Y|Z])^2] = sum_z (S_z - S_z^2/W_z).
    return {scheme: sum(s - s*s/w for w, s in groups.values() if w > 0)
            for scheme, groups in buckets.items()}


def aligned_patch_test(q_high=0.9, q_low=0.2, par=Parameters()):
    """Two matched coarse states; swap the complete toolkit across two sites.

    Both states: same two q values, same module/global phenotype marginals,
    one complete and one absent toolkit, and both currently white.
    """
    _unit("q_high", q_high)
    _unit("q_low", q_low)
    if q_high < q_low:
        raise ValueError("q_high must be >= q_low")
    values = {
        "aligned": 1 - (1 - outcome_probability(7, q_high, par)),
        "misaligned": 1 - (1 - outcome_probability(7, q_low, par)),
    }
    return {"q_high": q_high, "q_low": q_low, **values,
            "difference": values["aligned"] - values["misaligned"]}


def assurance_contrast(q, cost, low=0.1, high=0.9, par=Parameters()):
    """Survival/support/joint contrasts holding q, initial allele p and N fixed."""
    _unit("q", q)
    if high < low:
        raise ValueError("high assurance must be >= low assurance")
    results = {}
    for name, a in (("low", low), ("high", high)):
        p = replace(par, assurance=a, assurance_cost=cost)
        results[name] = {e: outcome_probability(7, q, p, e)
                         for e in ("survival", "genetic_support", "joint")}
    return {"q": q, "cost": cost, "low_assurance": low, "high_assurance": high,
            **results,
            "difference_high_minus_low": {k: results["high"][k] - results["low"][k]
                                          for k in results["low"]}}


def run_benchmarks():
    """Fixed protocol: no optimization or threshold selection from results."""
    p = Parameters()
    risks = brier_risks(par=p)
    null = brier_risks(par=p, null_history=True)
    neutral_assay = brier_risks(par=p, sensitivity=.5, specificity=.5)
    perfect_assay = brier_risks(par=p, sensitivity=1, specificity=1)
    aligned = aligned_patch_test(par=p)
    uniform_q = aligned_patch_test(q_high=.6, q_low=.6, par=p)
    assurance = [assurance_contrast(q, cost, par=p)
                 for cost in (0., .5, 1., 2.) for q in (.2, .5, .8)]
    closed = outcome_probability(7, .6, replace(p, allele_frequency=0.), "genetic_support")
    reopened_mutation = outcome_probability(7, .6, replace(p, allele_frequency=0., mutation_rate=.02), "genetic_support")
    reopened_migration = outcome_probability(7, .6, replace(p, allele_frequency=0., immigrant_fraction=.02), "genetic_support")
    return {
        "protocol": "SDEA exact synthetic benchmarks v1 (2026-10-26)",
        "scope": "CONSTRUCTED; no observed evolutionary regain or natural-island calibration",
        "parameters": vars(p),
        "history": {"even": module_marginals("even"), "odd": module_marginals("odd"),
                    "brier": risks, "null_history": null,
                    "uninformative_assay": neutral_assay, "perfect_assay": perfect_assay},
        "patch_alignment": {"comparison": aligned, "equal_q_negative_control": uniform_q},
        "assurance_tradeoff": assurance,
        "allele_reopening": {"closed_without_mutation_or_migration": closed,
                             "mutation_0_02": reopened_mutation,
                             "immigration_0_02": reopened_migration},
    }
