"""SDEA v6: finite propagules in a coupled three-patch synthetic metapopulation.

Genotype alleles (A/B) are hypothesized capacities, not phylogenetically
confirmed traits. Each occupied source has a HARD per-generation seed attempt
budget shared among local recruitment and ALL destination patches. Pollen
moves fathers only and cannot settle in a patch without a maternal plant.
Three independently assorting patches, two unlinked biallelic diploid loci.
No empirical calibration, maternal age, mutation, seed bank, or phenotype fitness.
"""
from collections import defaultdict
from dataclasses import dataclass, replace
from math import exp, isfinite, sqrt
from random import Random

GENOTYPES = tuple((a, b) for a in range(3) for b in range(3))
O = (0, 0)
A = (2, 0)
B = (0, 2)
IDENTITY = ((1., 0., 0.), (0., 1., 0.), (0., 0., 1.))


def _unit(name, value):
    if type(value) not in (float, int) or not isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{name} must lie in [0,1]")


def _matrix(name, value):
    if type(value) is not tuple or len(value) != 3:
        raise ValueError(f"{name} requires three rows")
    for row in value:
        if type(row) is not tuple or len(row) != 3:
            raise ValueError(f"{name} requires three probabilities per row")
        for p in row:
            _unit(name, p)
        if abs(sum(row) - 1) > 1e-12:
            raise ValueError(f"{name} row must sum to 1")


@dataclass(frozen=True)
class Budget:
    capacity: int = 2
    seed_attempts: int = 3
    seed_routing: tuple = IDENTITY  # each row is one maternal source -> destinations
    pollen_routing: tuple = IDENTITY  # each row is one maternal patch -> fathers
    q_low: float = .25
    q_high: float = .75
    high_stays_high: float = .80
    low_becomes_high: float = .12
    assurance: float = .45
    inbreeding_depression: float = .20
    catastrophe_rates: tuple = (.0, .0, .0)
    activation_a: float = .75
    activation_b: float = .75
    tie_a: float = .5
    expression_feedback_a: float = .0
    expression_feedback_b: float = .0
    selection_a: float = .0
    selection_b: float = .0

    def __post_init__(self):
        if type(self.capacity) is not int or not 1 <= self.capacity <= 5:
            raise ValueError("capacity must be integer 1..5")
        if type(self.seed_attempts) is not int or not 0 <= self.seed_attempts <= 32:
            raise ValueError("seed_attempts must be integer 0..32")
        _matrix('seed_routing', self.seed_routing)
        _matrix('pollen_routing', self.pollen_routing)
        if type(self.catastrophe_rates) is not tuple or len(self.catastrophe_rates) != 3:
            raise ValueError('catastrophe_rates must be triple')
        for n in ('q_low', 'q_high','high_stays_high','low_becomes_high',
                  'assurance','inbreeding_depression','activation_a','activation_b','tie_a'):
            _unit(n, getattr(self, n))
        for p in self.catastrophe_rates:
            _unit('catastrophe', p)
        if self.q_low > self.q_high:
            raise ValueError("low q must not exceed high q")
        for n in ('expression_feedback_a','expression_feedback_b'):
            v = getattr(self, n)
            if type(v) not in (float,int) or not isfinite(v) or abs(v)>1:
                raise ValueError('feedback must be finite in [-1,1]')
        for n in ('selection_a','selection_b'):
            v = getattr(self, n)
            if type(v) not in (float,int) or not isfinite(v) or abs(v)>20:
                raise ValueError('selection coefficient outside [-20,20]')


def _population(pop, capacity):
    if type(pop) is not tuple or len(pop) > capacity:
        raise ValueError('patch must be tuple of <=capacity diploid genotypes')
    if any(type(g) is not tuple or len(g)!=2 or any(type(x) is not int for x in g) or g not in GENOTYPES for g in pop):
        raise ValueError('genotypes must be two-locus diploid pairs')


def allele_count(pop, locus):
    if locus not in (0,1):
        raise ValueError('A/B locus only')
    return sum(g[locus] for g in pop)


def _choose(weights, rng):
    """Choose one index from a normalized nonnegative probability sequence."""
    t = rng.random()
    running = 0.
    for i, p in enumerate(weights):
        running += p
        if t < running:
            return i
    if abs(running-1)>1e-10:
        raise AssertionError('invalid categorical weights')
    return len(weights)-1


def _gamete(genotype, rng):
    return tuple(int(rng.random() < x/2.) for x in genotype)


def _cross(mom, dad, rng):
    a, b = _gamete(mom,rng),_gamete(dad,rng)
    return (a[0]+b[0],a[1]+b[1])


def _father(populations, mother_patch, maternal_index, routing, rng):
    """Draw paternal location, exclude the mother for local outcross.

    Selecting a father location with no eligible individual wastes this
    outcross attempt (there is no reweighting to another patch).
    """
    j=_choose(routing[mother_patch],rng)
    dads=populations[j]
    if j == mother_patch:
        idx=[k for k in range(len(dads)) if k!=maternal_index]
    else:
        idx=list(range(len(dads)))
    if not idx:
        return None,None
    return j,dads[idx[rng.randrange(len(idx))]]


def successful_seed_probability(populations, qi, source, par=Budget()):
    """Exact per-attempt viability for this given maternal-patch snapshot.

    Route of a paternal outcross is drawn first. Failed local and empty-
    source father choices do not get retried. Mean over maternal parents.
    """
    if type(source) is not int or source not in range(3):
        raise ValueError('source must be one of three patches')
    if type(populations) is not tuple or len(populations)!=3:
        raise ValueError('populations must have three patches')
    for patch in populations:
        _population(patch,par.capacity)
    if type(qi) is not tuple or len(qi)!=3 or any(type(x) is not int or x not in (0,1) for x in qi):
        raise ValueError('invalid q indices')
    if not populations[source]:
        return 0.
    q=(par.q_low,par.q_high)[qi[source]]
    s=par.assurance*(1-q)
    admissible=0.
    for j,weight in enumerate(par.pollen_routing[source]):
        admissible += weight*int(len(populations[j])>= (2 if j==source else 1))
    return s*(1-par.inbreeding_depression)+(1-s)*q*admissible


def exact_colonization_from_one_source(par, donor_patch=0, recipient_patch=1,
                                        donor_q='high', count_donor=2):
    """Exact one-generation colonization when only one patch is occupied.

    Assumes all source genotypes A/A at two loci or any viable genotype, no
    selection or other source, recipient vacant, and independent seed attempts.
    Capacity>=1 ensures any arriving viable seed can occupy the destination.
    """
    if donor_patch == recipient_patch or donor_patch not in range(3) or recipient_patch not in range(3):
        raise ValueError('distinct three-patch sites required')
    if type(count_donor) is not int or not 1<=count_donor<=par.capacity:
        raise ValueError('invalid donor size')
    qi=tuple(int(donor_q=='high') if i==donor_patch else 0 for i in range(3))
    if donor_q not in ('low','high'):
        raise ValueError('donor_q must be high or low')
    pops=tuple((A,)*count_donor if i==donor_patch else () for i in range(3))
    viable=successful_seed_probability(pops,qi,donor_patch,par)
    m=par.seed_routing[donor_patch][recipient_patch]
    arrival=1-(1-viable*m)**par.seed_attempts
    return {'source_viability_per_attempt':viable,
            'per_attempt_arrival_probability':viable*m,
            'colonized_before_catastrophe':arrival,
            'colonized_after_catastrophe':arrival*(1-par.catastrophe_rates[recipient_patch])}



def exact_shared_budget_colonization(par=Budget(), donor_patch=0, recipients=(1,2),
                                     donor_q='high', count_donor=2):
    """Exact two initially vacant recipients competing for ONE finite seed pool.

    A single maternal budget slot cannot provide a seed to both patches.
    Independent-catastrophe thinning is imposed after arrivals. The formula
    is exact conditional on one occupied source and no seed bank/immigration.
    """
    if type(recipients) is not tuple or len(recipients)!=2:
        raise ValueError('two recipient indices required')
    a,b=recipients
    if any(type(x) is not int or x not in range(3) or x==donor_patch for x in recipients) or a==b:
        raise ValueError('recipients must be distinct and exclude donor')
    first=exact_colonization_from_one_source(par,donor_patch,a,donor_q,count_donor)
    second=exact_colonization_from_one_source(par,donor_patch,b,donor_q,count_donor)
    ua=first['per_attempt_arrival_probability']
    ub=second['per_attempt_arrival_probability']
    if ua+ub>1+1e-12:
        raise AssertionError('overallocated source seed attempt')
    k=par.seed_attempts
    initial_joint=1-(1-ua)**k-(1-ub)**k+(1-ua-ub)**k
    both=max(0.,initial_joint)*(1-par.catastrophe_rates[a])*(1-par.catastrophe_rates[b])
    pa=first['colonized_after_catastrophe']
    pb=second['colonized_after_catastrophe']
    naive_independent=pa*pb
    return {'budget':k,'source_seed_success_probability':first['source_viability_per_attempt'],
            'per_attempt_receiver_arrival':[ua,ub],
            'recipient_colonization':[pa,pb],
            'joint_colonization':both,
            'naive_independent_joint':naive_independent,
            'joint_minus_independent':both-naive_independent,
            'source_retention':(1-(1-first['source_viability_per_attempt']*par.seed_routing[donor_patch][donor_patch])**k)*(1-par.catastrophe_rates[donor_patch])}


def _visible(genotypes, q, par, rng):
    result=[]
    for g in genotypes:
        ap=par.activation_a*q if g[0]>0 else 0.
        bp=par.activation_b*(1-q) if g[1]>0 else 0.
        has_a=rng.random()<ap
        has_b=rng.random()<bp
        if has_a and has_b:
            result.append('A' if rng.random()<par.tie_a else 'B')
        elif has_a:
            result.append('A')
        elif has_b:
            result.append('B')
        else:
            result.append('none')
    return tuple(result)


def _settle(pool, q, capacity, par, rng):
    """Capacity-limited weighted lottery without replacing seed individuals."""
    if len(pool)<=capacity:
        return tuple(pool)
    remaining=list(pool)
    selected=[]
    for _ in range(capacity):
        weights=[exp(par.selection_a*(q-.5)*g[0] +
                     par.selection_b*(.5-q)*g[1]) for (g,_,_) in remaining]
        total=sum(weights)
        i=_choose([w/total for w in weights],rng)
        selected.append(remaining.pop(i))
    return tuple(selected)


def generation(populations, qi, visible, par, rng):
    """Synchronous 3-patch update with FINITE source seed budgets.

    Returns new populations, new q, new visible phenotypes and seed accounting.
    A seed has exactly ONE maternal patch and ONE destination. Its paternal
    pollen may come from another patch independently, never replacing mother.
    """
    if type(populations) is not tuple or len(populations)!=3:
        raise ValueError('three populations required')
    for pop in populations:
        _population(pop,par.capacity)
    if type(qi) is not tuple or len(qi)!=3 or any(type(x) is not int or x not in (0,1) for x in qi):
        raise ValueError('three binary q states required')
    if type(visible) is not tuple or len(visible)!=3 or any(type(v) is not tuple or any(z not in ('none','A','B') for z in v) for v in visible):
        raise ValueError('three expression vectors required')
    seed_pools=[[],[],[]]
    attempts=[0,0,0]; viable=[0,0,0]; export=[0,0,0]; kept=[0,0,0]
    for i,adults in enumerate(populations):
        if not adults:
            continue
        q=(par.q_low,par.q_high)[qi[i]]
        s=par.assurance*(1-q)
        for _ in range(par.seed_attempts):
            attempts[i]+=1
            mother_index=rng.randrange(len(adults))
            mom=adults[mother_index]
            if rng.random()<s:
                if rng.random()>=1-par.inbreeding_depression:
                    continue
                dad=mom; father_site=i
            else:
                father_site,dad=_father(populations,i,mother_index,par.pollen_routing,rng)
                if dad is None or rng.random()>=q:
                    continue
            zygote=_cross(mom,dad,rng)
            viable[i]+=1
            j=_choose(par.seed_routing[i],rng)
            if j!=i:
                export[i]+=1
            else:
                kept[i]+=1
            seed_pools[j].append((zygote,i,father_site))
    survived=[]; admissions=[];visible_next=[]
    qnext=[]
    for i in range(3):
        p=par.high_stays_high if qi[i] else par.low_becomes_high
        p=min(1.,max(0.,p + (par.expression_feedback_a if 'A' in visible[i] else 0.)+
                            (par.expression_feedback_b if 'B' in visible[i] else 0.)))
        qnext.append(int(rng.random()<p))
    for i in range(3):
        q=(par.q_low,par.q_high)[qnext[i]]
        recruited=_settle(seed_pools[i],q,par.capacity,par,rng)
        if rng.random()<par.catastrophe_rates[i]:
            recruited=()
        survived.append(tuple(sorted(g for (g,_,_) in recruited)))
        admissions.append(tuple(recruited))
        visible_next.append(_visible(survived[-1],q,par,rng))
    return (tuple(survived),tuple(qnext),tuple(visible_next),
            {'attempts':tuple(attempts),'viable_produced':tuple(viable),
             'exports':tuple(export),'retained':tuple(kept),
             'pre_recruit_arrivals':tuple(len(p) for p in seed_pools),
             'after_recruit_admissions':tuple(admissions)})


def run_trajectory(initial, initial_q=(1,0,0), horizon=8,
                   par=Budget(), seed=0):
    if type(horizon) is not int or not 1<=horizon<=40:
        raise ValueError('horizon must be 1..40')
    if type(seed) is not int:
        raise ValueError('seed must be integer')
    if type(initial) is not tuple or len(initial)!=3:
        raise ValueError('initial must contain three patches')
    for p in initial:
        _population(p,par.capacity)
    if type(initial_q) is not tuple or len(initial_q)!=3 or any(type(x) is not int or x not in (0,1) for x in initial_q):
        raise ValueError('initial_q must be three binary indices')
    rng=Random(seed)
    pops=tuple(tuple(sorted(p)) for p in initial)
    qi=initial_q; visible=((),(),())
    seen_ever=[bool(pop) for pop in pops]
    extinction_seen=[False]*3
    first_colonized=[False]*3; recolonized=[False]*3
    ever_joint=[False]*3
    first_A=[False]*3; first_B=[False]*3
    total_viable=0;total_exports=0;total_attempts=0
    totals=[0,0,0];path=[]
    for t in range(1,horizon+1):
        old=pops
        pops,qi,visible,ledger=generation(pops,qi,visible,par,rng)
        total_attempts+=sum(ledger['attempts'])
        total_viable+=sum(ledger['viable_produced'])
        total_exports+=sum(ledger['exports'])
        for i in range(3):
            if old[i] and not pops[i]:
                extinction_seen[i]=True
            if not old[i] and pops[i]:
                if extinction_seen[i]:
                    recolonized[i]=True
                elif not seen_ever[i]:
                    first_colonized[i]=True
                seen_ever[i]=True
                if not any(mother!=i for g,mother,father in ledger['after_recruit_admissions'][i]):
                    raise AssertionError('empty patch can only receive seeds from elsewhere')
            ever_joint[i] |= any(g[0]>0 and g[1]>0 for g in pops[i])
            first_A[i] |= 'A' in visible[i]
            first_B[i] |= 'B' in visible[i]
            totals[i]+=len(pops[i])
        path.append({'generation':t,
                     'population_sizes':tuple(len(p) for p in pops),
                     'A_copies':tuple(allele_count(p,0) for p in pops),
                     'B_copies':tuple(allele_count(p,1) for p in pops),
                     'q_high':qi,
                     'seed_attempts':ledger['attempts'],
                     'viable_produced':ledger['viable_produced'],
                     'seed_exported':ledger['exports'],
                     'seed_retained':ledger['retained'],
                     'seed_received':ledger['pre_recruit_arrivals']})
    return {'horizon':horizon,'populations':pops,'ecology':qi,
            'colonized_initial_vacancy':tuple(first_colonized),
            'recolonized_after_extinction':tuple(recolonized),
            'ever_joint':tuple(ever_joint),
            'ever_visible_A':tuple(first_A),
            'ever_visible_B':tuple(first_B),
            'total_seed_attempts':total_attempts,
            'total_viable_seed_produced':total_viable,
            'total_seed_exports':total_exports,
            'time_integrated_occupancy':tuple(totals),
            'path':path}


def _wilson(hits,n,z=1.959963984540054):
    p=hits/n; d=1+z*z/n
    center=(p+z*z/(2*n))/d
    e=z*sqrt((p*(1-p)+z*z/(4*n))/n)/d
    return [max(0.,center-e),min(1.,center+e)]


def ensemble(initial,par=Budget(),initial_q=(1,0,0),horizon=8,
             replicates=3000,seed=11061026):
    if type(replicates) is not int or not 20<=replicates<=200000:
        raise ValueError('replicates must be integer 20..200000')
    hits=defaultdict(int)
    means=defaultdict(float)
    for i in range(replicates):
        v=run_trajectory(initial,initial_q,horizon,par,seed=seed+1000003*i)
        p=v['populations']
        flags={
            'patch0_occupied':bool(p[0]),'patch1_occupied':bool(p[1]),
            'patch2_occupied':bool(p[2]),'all_three_occupied':all(bool(x) for x in p),
            'any_patch_occupied':any(bool(x) for x in p),
            'global_extinction':not any(bool(x) for x in p),
            'patch0_A_supported':allele_count(p[0],0)>0,
            'patch1_A_supported':allele_count(p[1],0)>0,
            'patch2_A_supported':allele_count(p[2],0)>0,
            'patch0_B_supported':allele_count(p[0],1)>0,
            'patch2_B_supported':allele_count(p[2],1)>0,
            'globally_A_absent':not any(allele_count(x,0)>0 for x in p),
            'globally_B_absent':not any(allele_count(x,1)>0 for x in p),
            'patch0_joint':any(g[0]>0 and g[1]>0 for g in p[0]),
            'patch1_joint':any(g[0]>0 and g[1]>0 for g in p[1]),
            'patch2_joint':any(g[0]>0 and g[1]>0 for g in p[2]),
            'patch1_first_colonization':v['colonized_initial_vacancy'][1],
            'patch2_first_colonization':v['colonized_initial_vacancy'][2],
            'patch0_recolonized':v['recolonized_after_extinction'][0],
            'patch1_recolonized':v['recolonized_after_extinction'][1],
            'patch0_ever_joint':v['ever_joint'][0],
            'patch1_ever_joint':v['ever_joint'][1],
            'patch2_ever_joint':v['ever_joint'][2],
            'patch0_ever_visible_A':v['ever_visible_A'][0],
            'patch1_ever_visible_A':v['ever_visible_A'][1],
            'patch2_ever_visible_B':v['ever_visible_B'][2],
        }
        for key,x in flags.items():
            hits[key]+=int(x)
        for j,loc in enumerate(v['time_integrated_occupancy']):
            means[f'patch{j}_time_integrated_occupancy']+=loc
        means['viable_seeds']+=v['total_viable_seed_produced']
        means['exported_seeds']+=v['total_seed_exports']
    return {'replicates':replicates,'horizon':horizon,'seed':seed,
            'initial_populations':[[list(g) for g in x] for x in initial],
            'parameters':{'capacity':par.capacity,'seed_attempts':par.seed_attempts,
                          'seed_routing':[list(r) for r in par.seed_routing],
                          'pollen_routing':[list(r) for r in par.pollen_routing],
                          'catastrophe_rates':list(par.catastrophe_rates),
                          'selection_a':par.selection_a,'selection_b':par.selection_b,
                          'assurance':par.assurance,'inbreeding_depression':par.inbreeding_depression,
                          'q_low':par.q_low,'q_high':par.q_high},
            'event_counts':dict(sorted(hits.items())),
            'probabilities':{k:v/replicates for k,v in sorted(hits.items())},
            'mc_wilson_95':{k:_wilson(v,replicates) for k,v in sorted(hits.items())},
            'mean_events':{k:v/replicates for k,v in sorted(means.items())}}


def difference(cases,first,second,metric):
    x,y=cases[first],cases[second]
    a,b=x['probabilities'][metric],y['probabilities'][metric]
    na,nb=x['replicates'],y['replicates']
    d=b-a
    se=sqrt(a*(1-a)/na+b*(1-b)/nb)
    return {'first':first,'second':second,'endpoint':metric,
            'p_first':a,'p_second':b,'difference_second_minus_first':d,
            'mc_se':se,'mc_95_difference':[d-1.95996398454*se,d+1.95996398454*se]}


def benchmark(replicates=3000):
    """Prospectively declared synthetic frontier cells; no field data fit."""
    local=IDENTITY
    moderate=((.65,.175,.175),(.175,.65,.175),(.175,.175,.65))
    strong=((.20,.40,.40),(.40,.20,.40),(.40,.40,.20))
    chain=((.70,.30,0.),(.15,.70,.15),(0.,.30,.70))
    connected=((.65,.175,.175),(.175,.65,.175),(.175,.175,.65))
    # All cells have exactly three finite patches; seed pool budget is shared.
    empty=( (A,A), (), () )
    distinct=( (A,A), (), (B,B) )
    base=Budget()
    cells={
       'single_local_B3':(empty,replace(base,seed_attempts=3,seed_routing=local), (1,0,0)),
       'single_moderate_B3':(empty,replace(base,seed_attempts=3,seed_routing=moderate),(1,0,0)),
       'single_strong_B3':(empty,replace(base,seed_attempts=3,seed_routing=strong),(1,0,0)),
       'single_strong_B8':(empty,replace(base,seed_attempts=8,seed_routing=strong),(1,0,0)),
       'distinct_seed_chain_B3':(distinct,replace(base,seed_attempts=3,seed_routing=chain),(1,0,0)),
       'distinct_seed_pollen_chain_B3':(distinct,replace(base,seed_attempts=3,seed_routing=chain,pollen_routing=chain),(1,0,0)),
       'distinct_seed_pollen_connected_B3':(distinct,replace(base,seed_attempts=3,seed_routing=connected,pollen_routing=connected),(1,0,0)),
       'distinct_seed_pollen_connected_B3_selection':(distinct,replace(base,seed_attempts=3,seed_routing=connected,pollen_routing=connected,selection_a=1.4,selection_b=1.4),(1,0,0)),
    }
    results={}
    for i,(name,(initial,settings,q)) in enumerate(cells.items()):
        results[name]=ensemble(initial,settings,q,horizon=8,replicates=replicates,seed=719401+i*10000000)
    q=replace(base,seed_attempts=3,seed_routing=moderate)
    exact={
        'one_source_strong_budget3_joint_and_retention':exact_shared_budget_colonization(replace(q,seed_routing=strong)),
        'one_source_moderate_budget3_joint_and_retention':exact_shared_budget_colonization(q),
        'one_source_strong_budget1_joint_and_retention':exact_shared_budget_colonization(replace(q,seed_routing=strong,seed_attempts=1)),
        'one_source_strong_budget8_joint_and_retention':exact_shared_budget_colonization(replace(q,seed_routing=strong,seed_attempts=8)),
        'single_source_moderate_patch1':exact_colonization_from_one_source(q,0,1),
        'single_source_moderate_patch2':exact_colonization_from_one_source(q,0,2),
        'single_source_strong_patch1':exact_colonization_from_one_source(replace(q,seed_routing=strong),0,1),
        'single_source_strong_patch0_retained_probability':1-(1-successful_seed_probability(( (A,A),(),() ),(1,0,0),0,replace(q,seed_routing=strong))*strong[0][0])**3,
    }
    contrasts=[difference(results,*t) for t in (
        ('single_local_B3','single_strong_B3','patch0_occupied'),
        ('single_local_B3','single_strong_B3','patch1_first_colonization'),
        ('single_strong_B3','single_strong_B8','patch0_occupied'),
        ('single_strong_B3','single_strong_B8','patch1_first_colonization'),
        ('distinct_seed_chain_B3','distinct_seed_pollen_chain_B3','patch1_ever_joint'),
        ('distinct_seed_pollen_chain_B3','distinct_seed_pollen_connected_B3','patch1_ever_joint'),
        ('distinct_seed_pollen_connected_B3','distinct_seed_pollen_connected_B3_selection','patch0_A_supported'),
    )]
    return {'protocol':'SDEA v6 (2026-10-26) finite seed budgets, 3-patch pollen routing and retention',
            'status':'CONSTRUCTED Monte Carlo; no natural rates or original parent-repository mechanism reproduced',
            'replicate_unit':'whole 3-patch landscape for 8 generations',
            'exact_one_step':exact,'scenarios':results,'contrasts':contrasts}
