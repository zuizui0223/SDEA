"""SDEA v7b: three-patch directed-seed-reachability stress test.

A deliberately restricted, independently implemented subcase of SDEA v6's
finite-allocated Mendelian-seed generator: fixed q in every patch and year,
three patches, two unlinked diploid candidate route loci, finite source seed
attempts and capacity, no mutation, no genotype selection, no catastrophe.

This complement to the six-patch v7 alignment test compares seed routing OPERATORS, not inferred natural traits.
The first-step matched-state identity is proven analytically; the source-
recipient reachability result is a structural condition of seed routing.
Nothing here estimates ancestral floral regain, fitness or natural gene flow.
"""

from collections import defaultdict
from dataclasses import dataclass, replace
from itertools import product
from math import comb, isfinite, sqrt
from random import Random

A = (2, 0)
B = (0, 2)
O = (0, 0)
IDENTITY = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
CYCLE = ((.6, .4, 0.0), (0.0, .6, .4), (.4, 0.0, .6))
RETURN = ((.6, .4, 0.0), (.4, .6, 0.0), (.4, 0.0, .6))
REVERSE = ((.6, 0.0, .4), (.4, .6, 0.0), (0.0, .4, .6))
BALANCED = ((.6, .2, .2), (.2, .6, .2), (.2, .2, .6))
NETWORKS = {"cycle": CYCLE, "return": RETURN, "reverse": REVERSE,
            "balanced": BALANCED}


def _unit(name, value):
    if type(value) not in (int, float) or not isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{name} must be finite within [0,1]")


def _routing(name, rows):
    if type(rows) is not tuple or len(rows) != 3:
        raise ValueError(f"{name} must have three rows")
    for row in rows:
        if type(row) is not tuple or len(row) != 3:
            raise ValueError(f"{name} requires exactly 3 values per row")
        for x in row:
            _unit(name, x)
        if abs(sum(row) - 1) > 1e-12:
            raise ValueError(f"{name} row must sum to one")


@dataclass(frozen=True)
class Scenario:
    seed_routing: tuple = CYCLE
    pollen_routing: tuple = IDENTITY
    capacity: int = 2
    seed_attempts: int = 3
    q: float = .75
    assurance: float = .45
    inbreeding_depression: float = .20

    def __post_init__(self):
        _routing("seed_routing", self.seed_routing)
        _routing("pollen_routing", self.pollen_routing)
        if type(self.capacity) is not int or not 1 <= self.capacity <= 5:
            raise ValueError("capacity must be integer 1..5")
        if type(self.seed_attempts) is not int or not 0 <= self.seed_attempts <= 40:
            raise ValueError("seed_attempts must be integer 0..40")
        for name in ("q", "assurance", "inbreeding_depression"):
            _unit(name, getattr(self, name))


def _validate_population(pops, par):
    if type(pops) is not tuple or len(pops) != 3:
        raise ValueError("exactly three population tuples are required")
    for patch in pops:
        if type(patch) is not tuple or len(patch) > par.capacity:
            raise ValueError("patch population must be a tuple within capacity")
        for genotype in patch:
            if (type(genotype) is not tuple or len(genotype) != 2 or
                any(type(g) is not int or not 0 <= g <= 2 for g in genotype)):
                raise ValueError("genotype must be a two-locus diploid pair")


def _choose(row, rng):
    u = rng.random()
    for j, p in enumerate(row):
        if u < p:
            return j
        u -= p
    return len(row) - 1  # floating-point terminal guard


def _zygote(mother, father, rng):
    return tuple((rng.random() < a/2.) + (rng.random() < b/2.)
                 for a, b in zip(mother, father))


def viability(population_size, par=Scenario()):
    """Per maternal seed attempt, fixed-q viability with local pollen closure.

    Local outcrossing requires two distinct adults. If the paternal network
    is nonlocal this shortcut is invalid, hence exact tests require identity.
    """
    if type(population_size) is not int or population_size < 0:
        raise ValueError("population_size must be nonnegative integer")
    s = par.assurance * (1-par.q)
    return (s*(1-par.inbreeding_depression) if population_size >= 1 else 0.) + \
           ((1-s)*par.q if population_size >= 2 else 0.)


def step(populations, par, rng):
    """All outgoing propagules use previous adult populations synchronously.

    Shared seed budget: one successful seed is routed to exactly ONE patch.
    Failed paternal choices waste the corresponding non-selfing attempt.
    Capacity competition is unbiased sampling without replacement.
    """
    _validate_population(populations, par)
    arrivals = [[], [], []]
    attempted = [0, 0, 0]
    viable = [0, 0, 0]
    for source, mothers in enumerate(populations):
        if not mothers:
            continue
        s = par.assurance*(1-par.q)
        for _ in range(par.seed_attempts):
            attempted[source] += 1
            maternal_i = rng.randrange(len(mothers))
            mother = mothers[maternal_i]
            if rng.random() < s:
                if rng.random() >= 1 - par.inbreeding_depression:
                    continue
                father = mother
            else:
                father_site = _choose(par.pollen_routing[source], rng)
                fathers = populations[father_site]
                options = (tuple(f for i,f in enumerate(fathers) if i != maternal_i)
                           if father_site == source else fathers)
                if not options:
                    continue
                father = options[rng.randrange(len(options))]
                if rng.random() >= par.q:
                    continue
            child = _zygote(mother, father, rng)
            viable[source] += 1
            destination = _choose(par.seed_routing[source], rng)
            arrivals[destination].append(child)
    next_pops = []
    for pool in arrivals:
        if len(pool) > par.capacity:
            pool = rng.sample(pool, par.capacity)
        next_pops.append(tuple(sorted(pool)))
    ledger = {
        "attempted": tuple(attempted),
        "viable": tuple(viable),
        "arrived": tuple(len(x) for x in arrivals),
    }
    if sum(ledger["viable"]) != sum(ledger["arrived"]):
        raise AssertionError("a seed was lost or duplicated between source and destination")
    return tuple(next_pops), ledger


def reachable_sources(routing, initially_occupied):
    """Directed seed reachability, allowing self loops; no pollen rescue."""
    _routing("routing", routing)
    if not initially_occupied or any(i not in (0,1,2) for i in initially_occupied):
        raise ValueError("supply nonempty legitimate initially occupied set")
    reached = set(initially_occupied)
    agenda = list(reached)
    while agenda:
        source = agenda.pop()
        for j, m in enumerate(routing[source]):
            if m > 0 and j not in reached:
                reached.add(j)
                agenda.append(j)
    return tuple(sorted(reached))


def exact_matched_two_step(par=Scenario()):
    """Exact t1 focal and t2 patch2 colonization, one A source at patch0.

    Scope: initial A/A,A/A at patch0; patch1 & patch2 empty; q fixed in all
    patches; no external sources, mutation, selection, catastrophes, pollen
    migration; 0->2=0. At t1 only source0 can reproduce. At t2 only patch1
    can deliver a seed to patch2 via 1->2. 1->2=0 is a structural zero.
    """
    if par.capacity != 2 or par.pollen_routing != IDENTITY or par.seed_routing[0][2] != 0:
        raise ValueError("exact t2 expression requires K=2, local pollen and no 0→2 route")
    budget = par.seed_attempts
    S2 = viability(2, par)
    p01 = S2 * par.seed_routing[0][1]
    p00 = S2 * par.seed_routing[0][0]
    t1_p1 = 1-(1-p01)**budget
    t1_p0 = 1-(1-p00)**budget
    m12 = par.seed_routing[1][2]
    total = 0.
    for x in range(budget+1):
        prob_x = comb(budget,x)*p01**x*(1-p01)**(budget-x)
        n1 = min(x,par.capacity)
        source_success = viability(n1,par)
        t2_p2_given_x = 1-(1-m12*source_success)**budget
        total += prob_x*t2_p2_given_x
    return {"seed_budget":budget,"t1_patch0_occupied":t1_p0,
            "t1_patch1_colonized":t1_p1,"t1_patch2_colonized":0.,
            "t2_patch2_colonized":total,"source_viability":S2,
            "single_adult_viability":viability(1,par),
            "two_adult_viability":S2,
            "reachability":list(reachable_sources(par.seed_routing,(0,)))}



def operator_brier_gain(first=Scenario(seed_routing=CYCLE),
                        second=Scenario(seed_routing=RETURN)):
    """Exact Bayes Brier value of observing future seed routing at generation 1.

    Prior over the two imposed *future network laws* is 50/50. The complete
    observed current genotype state at generation 1 has the same distribution
    under both laws; n1 adult A/A count is sufficient for the t2 patch2
    occupancy event under the restricted one-source/locus/no-pollen setup.
    Because route0->patch2 is zero, only patch1 can colonize patch2 at t2.
    The network is an EXTERNAL OPERATOR: this does not violate Markov state
    sufficiency conditional on a known transition law.
    """
    for par in (first, second):
        if par.capacity != 2 or par.pollen_routing != IDENTITY or par.seed_routing[0][2] != 0:
            raise ValueError("Brier exercise restricted to K=2, local pollen, no 0→2")
    if (replace(first,seed_routing=second.seed_routing) != second or
        first.seed_routing[0] != second.seed_routing[0]):
        raise ValueError("hold all state and reproductive parameters and t1 routing fixed")
    budget=first.seed_attempts
    incoming=viability(2,first)*first.seed_routing[0][1]
    coarse=full=0.
    conditional=[]
    for n in (0,1,2):
        mass=sum(comb(budget,x)*incoming**x*(1-incoming)**(budget-x)
                 for x in range(budget+1) if min(2,x)==n)
        s_n=viability(n,first)
        y1=1-(1-first.seed_routing[1][2]*s_n)**budget
        y2=1-(1-second.seed_routing[1][2]*s_n)**budget
        pooled=(y1+y2)/2
        coarse+=mass*pooled*(1-pooled)
        full+=mass*.5*(y1*(1-y1)+y2*(1-y2))
        conditional.append({"patch1_diploids_at_t1":n,"probability":mass,
                            "t2_colonization_first":y1,"t2_colonization_second":y2})
    return {"operator_hidden_brier":coarse,"operator_observed_brier":full,
            "brier_gain_from_observing_operator":coarse-full,
            "prior_operator_weights":[.5,.5],"conditional_states":conditional,
            "note":"known full organism/genotype state is not future-sufficient if the external dispersal law is unknown"}


def exhaustive_two_step(par=Scenario()):
    """Independent 3^B categorical enumeration of t1 seed destinations.

    Enumeration is intentionally algebraically different from the binomial
    reduction above, for a diagnostic cross-check with B<=12.
    """
    if par.seed_attempts > 12 or par.capacity != 2 or par.pollen_routing != IDENTITY or par.seed_routing[0][2]:
        raise ValueError("exhaustive control only for bounded exact t2 cells")
    S = viability(2,par)
    probs = (1-S, S*par.seed_routing[0][0], S*par.seed_routing[0][1])
    ans = 0.
    for allocation in product(range(3),repeat=par.seed_attempts):
        prob = 1.
        arrivals1 = 0
        for outcome in allocation:
            prob *= probs[outcome]
            arrivals1 += (outcome == 2)
        if not prob:
            continue
        n=min(arrivals1,par.capacity)
        ans += prob*(1-(1-par.seed_routing[1][2]*viability(n,par))**par.seed_attempts)
    return ans


def run_history(initial, par, horizon=7, seed=0):
    """One independent whole-landscape history; no latent history invented."""
    if type(seed) is not int or type(horizon) is not int or not 1 <= horizon <= 40:
        raise ValueError("seed must be integer and horizon 1..40")
    _validate_population(initial,par)
    populations = tuple(tuple(sorted(p)) for p in initial)
    rng = Random(seed)
    ever = [bool(p) for p in populations]
    ever_joint = [any(a and b for a,b in p) for p in populations]
    history = []
    for g in range(1,horizon+1):
        populations, ledger = step(populations,par,rng)
        for j in range(3):
            ever[j] |= bool(populations[j])
            ever_joint[j] |= any(a and b for a,b in populations[j])
        history.append({"generation":g,"sizes":tuple(map(len,populations)),
                        "A":tuple(sum(x[0] for x in p) for p in populations),
                        "B":tuple(sum(x[1] for x in p) for p in populations),
                        "seed_attempted":ledger["attempted"],
                        "viable_produced":ledger["viable"],
                        "seeds_arrived":ledger["arrived"]})
    return {"final":populations,"ever_occupied":tuple(ever),
            "ever_joint":tuple(ever_joint),"trajectory":history}


def _wilson(k,n,z=1.959963984540054):
    p = k/n
    d = 1+z*z/n
    mid = (p+z*z/(2*n))/d
    radius=z*sqrt((p*(1-p)+z*z/(4*n))/n)/d
    return [max(0.,mid-radius),min(1.,mid+radius)]


def ensemble(initial, par, horizon=7, replicates=4000, seed=26071027):
    """Independent whole-landscape Monte Carlo; intervals are simulation only."""
    if type(replicates) is not int or not 50 <= replicates <= 200000:
        raise ValueError("replicates must be integer 50..200000")
    count = defaultdict(int)
    sums = defaultdict(float)
    for i in range(replicates):
        h = run_history(initial, par, horizon, seed=seed+1000003*i)
        final = h["final"]
        t1, t2 = h["trajectory"][0], h["trajectory"][min(1,horizon-1)]
        metrics = {
            "t1_patch0_occupied":bool(t1["sizes"][0]),
            "t1_patch1_colonized":bool(t1["sizes"][1]),
            "t1_patch2_colonized":bool(t1["sizes"][2]),
            "t2_patch2_occupied":bool(t2["sizes"][2]),
            "patch1_ever_occupied":h["ever_occupied"][1],
            "patch2_ever_occupied":h["ever_occupied"][2],
            "patch0_occupied_at_h":bool(final[0]),
            "patch1_occupied_at_h":bool(final[1]),
            "patch2_occupied_at_h":bool(final[2]),
            "all_three_occupied_at_h":all(bool(p) for p in final),
            "any_patch_occupied_at_h":any(bool(p) for p in final),
            "patch1_A_supported_at_h":any(a>0 for a,b in final[1]),
            "patch1_B_supported_at_h":any(b>0 for a,b in final[1]),
            "patch1_joint_individual_ever":h["ever_joint"][1],
            "patch1_joint_individual_at_h":any(a>0 and b>0 for a,b in final[1]),
            "global_A_absent_at_h":not any(a>0 for p in final for a,b in p),
            "global_B_absent_at_h":not any(b>0 for p in final for a,b in p),
        }
        for k,v in metrics.items():
            count[k] += bool(v)
        sums["realized_seed_attempts"] += sum(sum(y["seed_attempted"]) for y in h["trajectory"])
        sums["realized_viable_seeds"] += sum(sum(y["viable_produced"]) for y in h["trajectory"])
    return {"replicates":replicates,"horizon":horizon,"seed":seed,
            "initial":[[list(g) for g in p] for p in initial],
            "seed_routing":[list(v) for v in par.seed_routing],
            "pollen_routing":[list(v) for v in par.pollen_routing],
            "fixed_settings":{"capacity":par.capacity,"seed_attempts":par.seed_attempts,
                              "q":par.q,"assurance":par.assurance,
                              "inbreeding_depression":par.inbreeding_depression},
            "event_counts":dict(sorted(count.items())),
            "probabilities":{k:v/replicates for k,v in sorted(count.items())},
            "mc_wilson_95":{k:_wilson(v,replicates) for k,v in sorted(count.items())},
            "mean_endogenous":{"realized_seed_attempts":sums["realized_seed_attempts"]/replicates,
                               "realized_viable_seeds":sums["realized_viable_seeds"]/replicates}}


def mc_difference(cells, left, right, endpoint):
    a,b = cells[left], cells[right]
    x,y=a["probabilities"][endpoint],b["probabilities"][endpoint]
    na,nb=a["replicates"],b["replicates"]
    d=y-x
    se=sqrt(x*(1-x)/na+y*(1-y)/nb)
    return {"left":left,"right":right,"endpoint":endpoint,
            "difference_right_minus_left":d,"simulation_standard_error":se,
            "simulation_95_interval":[d-1.95996398454*se,d+1.95996398454*se]}


def benchmark(replicates=4000):
    """Declared route interventions at fixed genotypes, seed budgets and q.

    All network-only treatment families vary only seed_routing. Separate
    initial-population families are *not* pooled as ecological replicates.
    """
    single=((A,A),(),())
    distinct=((A,A),(),(B,B))
    null=((A,A),(),(A,A))
    configs={
        "single_cycle":(single,CYCLE),
        "single_return":(single,RETURN),
        "single_balanced":(single,BALANCED),
        "distinct_cycle":(distinct,CYCLE),
        "distinct_return":(distinct,RETURN),
        "distinct_balanced":(distinct,BALANCED),
        "neutral_cycle":(null,CYCLE),
        "neutral_reverse":(null,REVERSE),
    }
    cases={}
    for j,(name,(initial,routing)) in enumerate(configs.items()):
        cases[name]=ensemble(initial,Scenario(seed_routing=routing),horizon=7,
                             replicates=replicates,seed=27102026+j*101000003)
    exact={"cycle":exact_matched_two_step(Scenario(seed_routing=CYCLE)),
           "return":exact_matched_two_step(Scenario(seed_routing=RETURN)),
           "exhaustive_cycle_t2":exhaustive_two_step(Scenario(seed_routing=CYCLE)),
           "exhaustive_return_t2":exhaustive_two_step(Scenario(seed_routing=RETURN)),
           "network_operator_brier":operator_brier_gain(),
           "network_same_law_negative_control":operator_brier_gain(Scenario(seed_routing=CYCLE),Scenario(seed_routing=CYCLE)),
           "row_export":[1-CYCLE[i][i] for i in range(3)],
           "fixed_environment_q":.75,
           "initial_single_source":[[list(g) for g in p] for p in single]}
    comparisons=[mc_difference(cases,*c) for c in (
        ("single_return","single_cycle","t2_patch2_occupied"),
        ("single_return","single_cycle","patch2_ever_occupied"),
        ("single_return","single_cycle","all_three_occupied_at_h"),
        ("distinct_cycle","distinct_return","patch1_joint_individual_ever"),
        ("distinct_cycle","distinct_balanced","patch1_joint_individual_ever"),
        ("neutral_cycle","neutral_reverse","patch1_ever_occupied"),
    )]
    return {"protocol":"SDEA v7b matched first-census directed reachability (2026-10-27)",
            "status":"EXACT conditional theorems plus synthetic Monte Carlo; no field fit or broad mechanism replication",
            "replicate_unit":"whole 3-patch history; separate founder families not pooled",
            "exact":exact,"scenarios":cases,"contrasts":comparisons}
