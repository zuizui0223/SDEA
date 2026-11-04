"""SDEA v5: stochastic, COUPLED two-patch diploid gene-flow stress test.

Both patches are finite, change genotype and ecological state, may go extinct,
and can be recolonized by donor-produced seed; pollen requires resident mothers.
This is a wholly new, explicitly constructed kernel, not the ANTO/EGWE/izu-core
code, nor a calibrated natural plant population or measured re-evolution.

Each patch holds 0..2 diploid plants. Loci A and B independently assort;
A/B are candidate *genetic capacities*, not proven ancestral pigment routes.
A recruited plant may express either exclusive route (or neither), subject to
hypothetical, q-dependent activation and within-plant tie arbitration.
"""
from collections import defaultdict
from dataclasses import dataclass, replace
from functools import lru_cache
from itertools import product
from math import comb, isfinite, sqrt
from random import Random

GENOTYPES = tuple(product(range(3), repeat=2))
ROUTE_NAMES = ('none', 'A', 'B')
ORIGINS = ('empty', 'self_local', 'outcross_local', 'pollen_from_other', 'seed_from_other')


def _unit(name, p):
    if type(p) not in (int, float) or not isfinite(p) or not 0 <= p <= 1:
        raise ValueError(f'{name} must be finite and between 0 and 1')


@dataclass(frozen=True)
class Meta:
    capacity: int = 2
    pollen_rates: tuple = (.0, .0)  # per destination (other -> this patch), only for outcross
    seed_rates: tuple = (.0, .0)  # per destination, viable donor-seed pressure per slot
    q_low: float = .25
    q_high: float = .75
    high_stays_high: float = .80
    low_becomes_high: float = .12
    a_expression_q_feedback: float = .0
    b_expression_q_feedback: float = .0
    assurance: float = .45
    inbreeding_depression: float = .20
    activation_a: float = .75
    activation_b: float = .75
    tie_a: float = .5
    catastrophe_rates: tuple = (.0, .0)  # after local+seed recruitment, local shocks

    def __post_init__(self):
        if type(self.capacity) is not int or not 1 <= self.capacity <= 3:
            raise ValueError('capacity must be integer 1..3')
        for n in ('pollen_rates', 'seed_rates', 'catastrophe_rates'):
            v = getattr(self, n)
            if type(v) is not tuple or len(v) != 2:
                raise ValueError(f'{n} must be a pair of probabilities')
            for q in v:
                _unit(n, q)
        for n in ('q_low','q_high','high_stays_high','low_becomes_high',
                  'assurance','inbreeding_depression','activation_a','activation_b','tie_a'):
            _unit(n, getattr(self, n))
        if self.q_low > self.q_high:
            raise ValueError('q_low cannot exceed q_high')
        for n in ('a_expression_q_feedback','b_expression_q_feedback'):
            x = getattr(self, n)
            if type(x) not in (float, int) or not isfinite(x) or not -1 <= x <= 1:
                raise ValueError(f'{n} must be finite in [-1,1]')


def validate_population(pop, capacity):
    if type(pop) is not tuple or len(pop)>capacity:
        raise ValueError('population must be a tuple of <=capacity diploids')
    for g in pop:
        if type(g) is not tuple or g not in GENOTYPES or any(type(n) is not int for n in g):
            raise ValueError('invalid two-locus diploid genotype')


def allele_copies(pop, locus):
    if locus not in (0, 1):
        raise ValueError('locus must be A=0 or B=1')
    return sum(g[locus] for g in pop)


@lru_cache(maxsize=81)
def cross(m, f):
    """Exact, independent Mendelian assortment of the two unlinked loci."""
    if m not in GENOTYPES or f not in GENOTYPES:
        raise ValueError('invalid parental genotype')
    z = []
    for g in (0,1):
        p = m[g]/2
        q = f[g]/2
        z.append(((1-p)*(1-q), p*(1-q)+(1-p)*q, p*q))
    return tuple(z[0][g[0]]*z[1][g[1]] for g in GENOTYPES)


def weights(pop, q, assurance, depression, external_pollen, other_occupied):
    """UNCONDITIONAL probability weights for one local recruit-slot attempt.

    Distinct locally outcrossing parents required: 1 surviving adult alone may
    self or receive external pollen, but cannot have a distinct *local* father.
    The remaining probability is failure to produce a locally viable zygote.
    """
    _unit('q',q)
    if not pop:
        return (0.,0.,0.)
    s = assurance*(1-q)
    self_weight = s*(1-depression)
    immigrant_weight = (1-s)*q*external_pollen if other_occupied else 0.
    local_weight = (1-s)*q*(1-external_pollen if other_occupied else 1.) if len(pop)>=2 else 0.
    if self_weight+local_weight+immigrant_weight > 1.+1e-12:
        raise AssertionError('viability gate exceeds one')
    return self_weight, local_weight, immigrant_weight


def _zygote_terms(adults, fathers, category):
    """Enumerate one source's zygotes; exclude the mother from local outcross."""
    n=len(adults)
    out=[0.]*len(GENOTYPES)
    if category == 'self_local':
        for mom in adults:
            for i,x in enumerate(cross(mom,mom)):
                out[i]+=x/n
    elif category == 'outcross_local':
        if n<2:
            raise ValueError('local outcross requires at least 2 adults')
        for j,mom in enumerate(adults):
            for k,dad in enumerate(adults):
                if j != k:
                    for i,x in enumerate(cross(mom,dad)):
                        out[i]+=x/(n*(n-1))
    elif category == 'pollen_from_other':
        for mom in adults:
            for dad in fathers:
                for i,x in enumerate(cross(mom,dad)):
                    out[i]+=x/(n*len(fathers))
    else:
        raise ValueError('unrecognized zygote category')
    if abs(sum(out)-1)>1e-12:
        raise AssertionError('zygote distribution not normalized')
    return tuple(out)


@lru_cache(maxsize=20000)
def slot_distribution(recipient, donor, qi, qj, destination, par=Meta()):
    """Exact per-slot (origin, genotype) probability from PREVIOUS patches.

    A viable donor seed is first offered with probability m_s*S_donor_local,
    where S_donor_local is source reproduction without imported pollen. If no
    donor seed arrives, a local reproductive attempt occurs. This precedence
    prevents a migrant from requiring a living resident mother; it also makes
    a seed rescue possible if the recipient is empty. No external reservoir.
    """
    validate_population(recipient,par.capacity)
    validate_population(donor,par.capacity)
    if qi not in (0,1) or qj not in (0,1) or destination not in (0,1):
        raise ValueError('invalid ecology index or target patch')
    q=(par.q_low,par.q_high)[qi]
    qdon=(par.q_low,par.q_high)[qj]
    rate_p=par.pollen_rates[destination]
    rate_s=par.seed_rates[destination]
    weight_d=weights(donor,qdon,par.assurance,par.inbreeding_depression,0.,False)
    source_S=sum(weight_d)
    seed_probability=rate_s*source_S
    probs=defaultdict(float)
    if seed_probability:
        # donor's local gametes only: no circular dependence on recipient pollen
        for category,weight in zip(('self_local','outcross_local'),weight_d[:2]):
            if weight:
                dist=_zygote_terms(donor,donor,category)
                for g,p in zip(GENOTYPES,dist):
                    probs[('seed_from_other',g)]+=rate_s*weight*p
    local=weights(recipient,q,par.assurance,par.inbreeding_depression,rate_p,bool(donor))
    for category,weight in zip(('self_local','outcross_local','pollen_from_other'),local):
        if weight:
            dist=_zygote_terms(recipient,donor,category)
            for g,p in zip(GENOTYPES,dist):
                probs[(category,g)]+=(1-seed_probability)*weight*p
    probs[('empty',None)]+=(1-seed_probability)*(1-sum(local))
    result=tuple((origin,g,p) for (origin,g),p in sorted(probs.items()) if p>0)
    if abs(sum(p for _,_,p in result)-1)>5e-12:
        raise AssertionError('slot probabilities do not sum to one')
    return result


def _sample_distribution(dist, rng):
    u=rng.random()
    running=0.
    for origin,g,w in dist:
        running+=w
        if u<running:
            return origin,g
    # Last-entry floating point guard, not a new probability mass.
    return dist[-1][0],dist[-1][1]


def expression_probs(genotype, q, par=Meta()):
    """Exclusive visible A, B, neither. Tie resolved with par.tie_a."""
    if genotype not in GENOTYPES:
        raise ValueError('invalid genotype')
    _unit('q',q)
    a=par.activation_a*q if genotype[0]>0 else 0.
    b=par.activation_b*(1-q) if genotype[1]>0 else 0.
    both=a*b
    result=((1-a)*(1-b),a*(1-b)+par.tie_a*both,b*(1-a)+(1-par.tie_a)*both)
    if abs(sum(result)-1)>1e-12:
        raise AssertionError('visible genotype probabilities not normalized')
    return result


def _sample_visible(genotype,q,par,rng):
    p=expression_probs(genotype,q,par)
    u=rng.random()
    if u<p[0]:
        return 'none'
    if u<p[0]+p[1]:
        return 'A'
    return 'B'


def ecological_next(qi,visible,par,rng):
    if qi not in (0,1):
        raise ValueError('invalid qi')
    base=par.high_stays_high if qi else par.low_becomes_high
    p=min(1.,max(0.,base+(par.a_expression_q_feedback if 'A' in visible else 0.)+
                      (par.b_expression_q_feedback if 'B' in visible else 0.)))
    return int(rng.random()<p)


def generation(populations, ecology, visible, par, rng):
    """Synchronous one-generation update of BOTH source and destination patches."""
    if not isinstance(populations,tuple) or len(populations)!=2:
        raise ValueError('two patch populations required')
    for adults in populations:
        validate_population(adults,par.capacity)
    if not isinstance(ecology,tuple) or len(ecology)!=2 or any(x not in (0,1) for x in ecology):
        raise ValueError('two q-state indices required')
    if not isinstance(visible,tuple) or len(visible)!=2:
        raise ValueError('invalid visible-state flags')
    news=[];origin_cohorts=[];reentries=[]
    for i in range(2):
        j=1-i
        distribution=slot_distribution(populations[i],populations[j],ecology[i],ecology[j],i,par)
        offspring=[]; origins=[]; new_gene_entry=[]
        for _ in range(par.capacity):
            origin, genotype=_sample_distribution(distribution,rng)
            if genotype is not None:
                offspring.append(genotype)
                origins.append(origin)
                for locus, route in enumerate(('A','B')):
                    if (allele_copies(populations[i],locus)==0 and genotype[locus]>0):
                        if allele_copies(populations[j],locus)==0:
                            raise AssertionError('novel route appeared without donor or mutation')
                        if origin not in ('pollen_from_other','seed_from_other'):
                            raise AssertionError('route reentry must have external germline origin')
                        new_gene_entry.append((route,origin))
        if rng.random()<par.catastrophe_rates[i]:
            offspring=[];origins=[];new_gene_entry=[]
        news.append(tuple(sorted(offspring)))
        origin_cohorts.append(tuple(origins))
        reentries.append(tuple(new_gene_entry))
    qi_next=tuple(ecological_next(ecology[i],visible[i],par,rng) for i in range(2))
    expression=[]
    for i in range(2):
        q=(par.q_low,par.q_high)[qi_next[i]]
        expression.append(tuple(_sample_visible(g,q,par,rng) for g in news[i]))
    return (tuple(news),qi_next,tuple(expression),tuple(origin_cohorts),tuple(reentries))


def run_trajectory(initial, initial_q=(0,1), horizon=7, par=Meta(), seed=0):
    """One seeded landscape history; all endpoints use one independent replicate."""
    if type(horizon) is not int or not 1<=horizon<=30:
        raise ValueError('horizon must be 1..30')
    if type(seed) is not int:
        raise ValueError('seed must be integer')
    if not isinstance(initial,tuple) or len(initial)!=2:
        raise ValueError('two initial populations required')
    for pop in initial:
        validate_population(pop,par.capacity)
    rng=Random(seed)
    populations=tuple(tuple(sorted(p)) for p in initial)
    qi=initial_q
    visible=((),())  # no observed expression at t=0
    recolon_count=[0,0]; colonization_count=[0,0]; extinction_count=[0,0]
    ever_lost_after_occupancy=[False,False]
    entry={f'{r}_{o}_{i}':False for i in range(2) for r in ('A','B') for o in ('pollen','seed')}
    first_visible=['none','none']
    ever_A=[False,False];ever_B=[False,False];ever_joint=[False,False]
    route_A_when_absent=[False,False];route_B_when_absent=[False,False]
    path=[]
    for t in range(1,horizon+1):
        old=populations
        populations,qi,visible,origins,events=generation(populations,qi,visible,par,rng)
        for i in range(2):
            if old[i] and not populations[i]:
                extinction_count[i]+=1
                ever_lost_after_occupancy[i]=True
            if not old[i] and populations[i]:
                if ever_lost_after_occupancy[i]:
                    recolon_count[i]+=1
                else:
                    colonization_count[i]+=1  # first occupancy of initially empty patch
                if not any(o=='seed_from_other' for o in origins[i]):
                    raise AssertionError('colonization requires seed from other')
            route_a='A' in visible[i]
            route_b='B' in visible[i]
            if first_visible[i]=='none':
                first_visible[i]='both' if route_a and route_b else ('A' if route_a else ('B' if route_b else 'none'))
            ever_A[i]|=route_a;ever_B[i]|=route_b
            ever_joint[i]|=any(g[0] and g[1] for g in populations[i])
            for route,origin in events[i]:
                entry[f'{route}_{"pollen" if origin=="pollen_from_other" else "seed"}_{i}']=True
            if route_a and not allele_copies(old[i],0):route_A_when_absent[i]=True
            if route_b and not allele_copies(old[i],1):route_B_when_absent[i]=True
        path.append({
            'generation':t,
            'population_sizes':tuple(len(p) for p in populations),
            'q_high':qi,
            'allele_A':tuple(allele_copies(p,0) for p in populations),
            'allele_B':tuple(allele_copies(p,1) for p in populations),
            'visible':visible,
        })
    return {'horizon':horizon,'populations':populations,'ecology':qi,'first_visible':tuple(first_visible),
            'ever_A':tuple(ever_A),'ever_B':tuple(ever_B),'ever_joint_individual':tuple(ever_joint),
            'first_visible_route_newly_supported_A':tuple(route_A_when_absent),
            'first_visible_route_newly_supported_B':tuple(route_B_when_absent),
            'colonized_initial_vacancy':tuple(bool(n) for n in colonization_count),
            'recolonized_after_extinction':tuple(bool(n) for n in recolon_count),
            'colonization_count':tuple(colonization_count),'recolon_count':tuple(recolon_count),
            'extinction_count':tuple(extinction_count), 'entry':entry,'path':path}


def _wilson(num,total,z=1.959963984540054):
    ph=num/total
    den=1+z*z/total
    ctr=(ph+z*z/(2*total))/den
    half=z*sqrt((ph*(1-ph)+z*z/(4*total))/total)/den
    return [max(0.,min(ph,ctr-half)),min(1.,max(ph,ctr+half))]


def ensemble(initial,par=Meta(),initial_q=(0,1),horizon=7,replicates=5000,seed=31062026):
    """Independent seeded finite-population trajectories, not exact long-run DP.

    Per-hypothesis Wilson intervals describe repeated RNG trajectories only;
    they are NOT empirical uncertainty or posterior evolutionary rates.
    """
    if type(replicates) is not int or not 10<=replicates<=100000:
        raise ValueError('replicates must be integer 10..100000')
    if type(seed) is not int:
        raise ValueError('seed must be integer')
    hits=defaultdict(int)
    sums=defaultdict(float)
    for r in range(replicates):
        v=run_trajectory(initial,initial_q,horizon,par,seed=seed+1000003*r)
        a=v['populations']
        f0,f1=v['first_visible']
        metrics={
            'both_patches_occupied':bool(a[0]) and bool(a[1]),
            'both_patches_extinct':not a[0] and not a[1],
            'patch0_occupied':bool(a[0]),
            'patch1_occupied':bool(a[1]),
            'patch0_colonized_initial_vacancy':v['colonized_initial_vacancy'][0],
            'patch1_colonized_initial_vacancy':v['colonized_initial_vacancy'][1],
            'patch0_recolonized_after_extinction':v['recolonized_after_extinction'][0],
            'patch1_recolonized_after_extinction':v['recolonized_after_extinction'][1],
            'patch0_A_supported':allele_copies(a[0],0)>0,
            'patch0_B_supported':allele_copies(a[0],1)>0,
            'patch0_joint_in_individual':any(g[0]>0 and g[1]>0 for g in a[0]),
            'patch1_joint_in_individual':any(g[0]>0 and g[1]>0 for g in a[1]),
            'both_patches_joint_in_individual':all(any(g[0]>0 and g[1]>0 for g in p) for p in a),
            'global_A_absent':not any(allele_copies(p,0)>0 for p in a),
            'global_B_absent':not any(allele_copies(p,1)>0 for p in a),
            'patch0_any_A_ever':v['ever_A'][0],
            'patch0_any_B_ever':v['ever_B'][0],
            'patch0_joint_ever':v['ever_joint_individual'][0],
            'patch1_joint_ever':v['ever_joint_individual'][1],
            'patch0_A_first':f0=='A',
            'patch0_B_first':f0=='B',
            'patch0_simultaneous_first':f0=='both',
            'patch0_ancestral_gene_entered_via_pollen':v['entry']['A_pollen_0'],
            'patch0_ancestral_gene_entered_via_seed':v['entry']['A_seed_0'],
            'patch0_alternative_gene_entered_via_pollen':v['entry']['B_pollen_0'],
            'patch0_alternative_gene_entered_via_seed':v['entry']['B_seed_0'],
            'patch0_q_high_at_horizon':bool(v['ecology'][0]),
        }
        for key,value in metrics.items():
            hits[key]+=bool(value)
        sums['patch0_extinctions']+=v['extinction_count'][0]
        sums['patch0_initial_colonizations']+=v['colonization_count'][0]
        sums['patch0_recolonizations']+=v['recolon_count'][0]
    values={key:count/replicates for key,count in sorted(hits.items())}
    return {
        'replicates':replicates,'horizon':horizon,'seed':seed,
        'initial_q':list(initial_q),'parameters':vars(par),
        'initial_populations':[[list(g) for g in patch] for patch in initial],
        'probabilities':values,
        'wilson_95':{key:_wilson(hits[key],replicates) for key in sorted(hits)},
        'event_counts':dict(sorted(hits.items())),
        'mean_events':{k:v/replicates for k,v in sorted(sums.items())},
    }


def exact_recolonization_probability_empty(donor,qi_donor,seed_rate,par=Meta()):
    """Exact one-generation probability of seed-only recolonization of empty patch.

    Recruits are drawn independently with seed probability m_s*S_donor_local;
    a local catastrophe then can erase all recruits.
    """
    validate_population(donor,par.capacity)
    _unit('seed_rate',seed_rate)
    if qi_donor not in (0,1):
        raise ValueError('invalid donor q index')
    q=(par.q_low,par.q_high)[qi_donor]
    S=sum(weights(donor,q,par.assurance,par.inbreeding_depression,0.,False))
    return (1-(1-seed_rate*S)**par.capacity)*(1-par.catastrophe_rates[0])



def monte_carlo_difference(cases, first, second, metric):
    """Two independent scenario trajectories: second - first with 95% Wald MC CI.

    Not an inferential CI for any *biological population* or a multiple-testing-
    corrected finding. Binomial outcomes count independent landscape histories.
    """
    a,b=cases[first],cases[second]
    pa,pb=a['probabilities'][metric],b['probabilities'][metric]
    na,nb=a['replicates'],b['replicates']
    diff=pb-pa
    standard_error=sqrt(pa*(1-pa)/na + pb*(1-pb)/nb)
    return {'first':first,'second':second,'endpoint':metric,
            'p_first':pa,'p_second':pb,'difference_second_minus_first':diff,
            'mc_standard_error':standard_error,
            'mc_95_difference_interval':[max(-1.,diff-1.95996398454*standard_error),
                                         min(1.,diff+1.95996398454*standard_error)],
            'replicates_first':na,'replicates_second':nb,
            'method':'independent landscape replicates; approximate Wald interval; exploratory synthetic comparison'}


def exact_first_step_contrasts(par=Meta()):
    """Analytical operator counterexamples, no Monte Carlo or natural calibration."""
    A,B=(2,0),(0,2)
    mothers=(A,A);fathers=(B,B)
    cp=replace(par,pollen_rates=(.45,.45))
    cs=replace(par,seed_rates=(.35,.35))
    # On a low-q ancestral resident, alternative-only source can supply
    # migrant pollen that yields joint-capable individuals at generation 1.
    pollen=slot_distribution(mothers,fathers,0,1,0,cp)
    seed=slot_distribution(mothers,fathers,0,1,0,cs)
    base=slot_distribution(mothers,fathers,0,1,0,par)
    def outcome(d):
        joint=sum(w for o,g,w in d if g is not None and g[0] and g[1])
        ancestral=sum(w for o,g,w in d if g is not None and g[0])
        occupied=sum(w for o,g,w in d if g is not None)
        return {'per_recruit_joint_competence':joint,
                'any_joint_competence_generation1':1-(1-joint)**par.capacity,
                'any_ancestral_support_generation1':1-(1-ancestral)**par.capacity,
                'any_occupancy_generation1':1-(1-occupied)**par.capacity,
                'per_recruit_viability':occupied}
    return {'no_flow':outcome(base),'pollen_only':outcome(pollen),'seed_only':outcome(seed),
            'identity':'first-generation genotype/recruitment access only; independent recruit slots; source and receiver A/B homozygotes'}

def benchmark(replicates=9000):
    """Prospectively declared v5 synthetic stress cells; no empirical fitting."""
    O=(0,0); A=(2,0); B=(0,2)
    empty=((),(A,A))
    distinct=((A,A),(B,B))
    p=Meta()
    scenarios={
        'empty_no_flow': (empty,p,(0,1)),
        'empty_pollen_only': (empty,replace(p,pollen_rates=(.45,.45)),(0,1)),
        'empty_seed_only': (empty,replace(p,seed_rates=(.35,.35)),(0,1)),
        'empty_seed_source_fragile': (empty,replace(p,seed_rates=(.35,.35),catastrophe_rates=(0.,.32)),(0,1)),
        'distinct_no_flow': (distinct,p,(0,1)),
        'distinct_pollen_only': (distinct,replace(p,pollen_rates=(.45,.45)),(0,1)),
        'distinct_seed_only': (distinct,replace(p,seed_rates=(.35,.35)),(0,1)),
        'distinct_both_flows': (distinct,replace(p,pollen_rates=(.45,.45),seed_rates=(.35,.35)),(0,1)),
        'distinct_both_with_q_feedback': (distinct,replace(p,pollen_rates=(.45,.45),seed_rates=(.35,.35),a_expression_q_feedback=.12,b_expression_q_feedback=-.12),(0,1)),
    }
    outputs={}
    for i,(name,(initial,param,initial_q)) in enumerate(scenarios.items()):
        outputs[name]=ensemble(initial,param,initial_q,horizon=7,replicates=replicates,seed=719303+i*100000000)
    contrasts=[
        monte_carlo_difference(outputs,'empty_pollen_only','empty_seed_only','patch0_colonized_initial_vacancy'),
        monte_carlo_difference(outputs,'empty_seed_source_fragile','empty_seed_only','patch0_colonized_initial_vacancy'),
        monte_carlo_difference(outputs,'empty_seed_source_fragile','empty_seed_only','patch0_A_supported'),
        monte_carlo_difference(outputs,'distinct_no_flow','distinct_pollen_only','patch0_joint_ever'),
        monte_carlo_difference(outputs,'distinct_seed_only','distinct_both_flows','both_patches_occupied'),
        monte_carlo_difference(outputs,'distinct_seed_only','distinct_both_flows','patch0_joint_in_individual'),
        monte_carlo_difference(outputs,'distinct_both_flows','distinct_both_with_q_feedback','patch0_q_high_at_horizon'),
    ]
    return {
        'protocol':'SDEA v5 (2026-10-26) dynamic bidirectional 2-patch synthetic bridge',
        'status':'CONSTRUCTED Monte Carlo, NOT cross-source replication, empirical fit, or long-run establishment',
        'source':'both finite patches update simultaneously; no static source reservoir, no external inputs',
        'recolonization_exact':{
            'donor_two_A_high_q_seed_035':exact_recolonization_probability_empty((A,A),1,.35,p),
            'donor_two_A_high_q_seed_0':exact_recolonization_probability_empty((A,A),1,0.,p),
            'no_viable_donor':exact_recolonization_probability_empty((),1,.35,p),
        },
        'scenarios':outputs,
        'exact_first_step':exact_first_step_contrasts(p),
        'exploratory_mc_contrasts':contrasts,
    }
