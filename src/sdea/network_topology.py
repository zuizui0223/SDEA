"""SDEA v7: matched-budget, degree-balanced network/trait-alignment test.

Not field data or a port of ANTO, EGWE, or izu-core. The six-patch
finite-diploid generator fixes seed attempt B, q, founder genotypes, node degree,
row/column seed flow and seed export proportion across the focal topology pairs.
Realized *future* attempts may diverge because network structure changes occupancy.
Two abstract unlinked biallelic loci (A,B) represent candidate functional resources,
NOT identified ancestral versus alternative pigment pathways.
"""
from dataclasses import dataclass, replace
from math import sqrt, isfinite
from random import Random

NODES = 6
A = (2, 0)
B = (0, 2)
O = (0, 0)
ALLELES = (A, B, O)
INIT = ((A,A), (A,A), (A,A), (B,B), (B,B), (B,B))
IDENTITY = tuple(tuple(float(i==j) for j in range(NODES)) for i in range(NODES))
ORDERS = {
    'clustered': ((0,1,2,3,4,5),),
    'alternating': ((0,3,1,4,2,5),),
    'components': ((0,1,2), (3,4,5)),
}


def routing(name, retention=.25):
    """Undirected 2-regular topology, matched outgoing/incoming movement budgets.

    Every vertex has two distinct neighbors, including in disjoint triangles;
    all rows and columns sum to one, with identical diagonal retention.
    """
    if name == 'identity':
        return IDENTITY
    if name not in ORDERS:
        raise ValueError('unknown named network')
    if type(retention) not in (float, int) or not isfinite(retention) or not 0<=retention<=1:
        raise ValueError('invalid local retention')
    rows=[[0.]*NODES for _ in range(NODES)]
    for i in range(NODES):
        rows[i][i]=retention
    for cycle in ORDERS[name]:
        n=len(cycle)
        if n < 3: raise ValueError('cycles must have at least three sites')
        for k,i in enumerate(cycle):
            for j in (cycle[(k+1)%n], cycle[(k-1)%n]):
                rows[i][j]+=(1-retention)/2
    ans=tuple(tuple(v) for v in rows)
    validate_matrix(ans)
    return ans


def validate_matrix(matrix):
    if type(matrix) is not tuple or len(matrix)!=NODES:
        raise ValueError('routing matrix must be 6x6 tuple')
    for row in matrix:
        if type(row) is not tuple or len(row)!=NODES:
            raise ValueError('routing matrix must be 6x6 tuple')
        if any(type(p) not in (float,int) or not isfinite(p) or not 0<=p<=1 for p in row):
            raise ValueError('invalid routing probability')
        if abs(sum(row)-1.)>1e-12:
            raise ValueError('routing row must sum to one')


def mixing_signature(matrix):
    """Control: incoming and outgoing totals, local shares and inter-route edges."""
    validate_matrix(matrix)
    a=set(range(3)); b=set(range(3,6))
    return {
        'outgoing':[sum(row) for row in matrix],
        'incoming':[sum(matrix[i][j] for i in range(NODES)) for j in range(NODES)],
        'retention':[matrix[i][i] for i in range(NODES)],
        'cross_route_mass':sum(matrix[i][j] for i in a for j in b) + sum(matrix[i][j] for i in b for j in a),
        'cross_route_edges':sum(matrix[i][j]>0 for i in a for j in b) + sum(matrix[i][j]>0 for i in b for j in a),
        'undirected_edges':sum(matrix[i][j]>0 for i in range(NODES) for j in range(i+1,NODES)),
    }


@dataclass(frozen=True)
class Network:
    capacity: int=2
    seed_attempts: int=3
    q: float=.75
    assurance: float=.45
    inbreeding_depression: float=.20
    catastrophe: float=0.
    activation_A: float=.75
    activation_B: float=.75
    seed: tuple=IDENTITY
    pollen: tuple=IDENTITY

    def __post_init__(self):
        if type(self.capacity) is not int or not 1<=self.capacity<=5:
            raise ValueError('capacity must be integer 1..5')
        if type(self.seed_attempts) is not int or not 0<=self.seed_attempts<=32:
            raise ValueError('seed_attempts must be integer 0..32')
        for k in ('q','assurance','inbreeding_depression','catastrophe','activation_A','activation_B'):
            v=getattr(self,k)
            if type(v) not in (float,int) or not isfinite(v) or not 0<=v<=1:
                raise ValueError('invalid probability '+k)
        validate_matrix(self.seed)
        validate_matrix(self.pollen)


def validate_populations(populations, par=Network()):
    if type(populations) is not tuple or len(populations)!=NODES:
        raise ValueError('expected six populations')
    for pop in populations:
        if type(pop) is not tuple or len(pop)>par.capacity:
            raise ValueError('adult population exceeds capacity or wrong type')
        if any(type(g) is not tuple or len(g)!=2 or any(type(x) is not int or not 0<=x<=2 for x in g) for g in pop):
            raise ValueError('expected biallelic diploid at two loci')


def _choose(weights, rng):
    r=rng.random(); total=0.
    for i,p in enumerate(weights):
        total+=p
        if r<total:
            return i
    if abs(total-1.)>1e-11:
        raise AssertionError('probabilities are not normalized')
    return len(weights)-1


def _index(n,rng):
    """Use Random.random only: portable deterministic integer draw across Pythons."""
    if n<=0: raise ValueError('cannot sample an empty sequence')
    return min(n-1,int(rng.random()*n))


def _gamete(g,rng):
    return (int(rng.random()<g[0]/2),int(rng.random()<g[1]/2))


def _cross(mom,dad,rng):
    a=_gamete(mom,rng); b=_gamete(dad,rng)
    return (a[0]+b[0],a[1]+b[1])


def _first_hybrid_exact(par=Network(), initial=INIT):
    """Exact probability that >=1 A+B diploid zygote is *produced* at t=1.

    Formula is restricted to the frozen source homozygotes in INIT, two adults
    each and no route-complete parental adults. Seed routing is irrelevant to
    production (but affects recruitment); all viable cross-type A×B have (1,1).
    """
    if initial!=INIT or par.capacity!=2:
        raise ValueError('closed-form witness is restricted to the v7 founders and capacity')
    s=par.assurance*(1-par.q)
    viable_outcross=(1-s)*par.q
    hit_probs=[]
    for i in range(NODES):
        other=range(3,6) if i<3 else range(3)
        cross_mass=sum(par.pollen[i][j] for j in other)
        hit_probs.append(viable_outcross*cross_mass)
    p_none=1.
    for p in hit_probs:
        p_none *= (1-p)**par.seed_attempts
    return {
        'probability_any_hybrid_seed_produced':1-p_none,
        'expected_hybrid_seeds_produced':par.seed_attempts*sum(hit_probs),
        'per_source_viable_hybrid_per_attempt':hit_probs,
        'attempt_budget_initial_total':par.seed_attempts*NODES,
        'exact_given_all_initial_patches_occupied':True,
    }



def network_prediction_risk(budget=3, par=Network()):
    """Exact prior-balanced Brier score for clustered vs alternating ring labels.

    Observation: one binary indicator of *any viable hybrid seed produced at t=1*.
    Their observed initial genotype counts, route totals, node degrees and q
    are matched. Full knowledge of the routing policy reduces Bayes Brier risk
    by Var(two conditional event probabilities). This is a constructed model,
    not observed predictive accuracy or biological ancestral-history memory.
    """
    if type(budget) is not int or not 0<=budget<=32:
        raise ValueError('budget must be integer 0..32')
    c=_first_hybrid_exact(replace(par,seed_attempts=budget,
                                  pollen=routing('clustered')))
    a=_first_hybrid_exact(replace(par,seed_attempts=budget,
                                  pollen=routing('alternating')))
    pc=c['probability_any_hybrid_seed_produced']
    pa=a['probability_any_hybrid_seed_produced']
    mean=(pc+pa)/2
    coarse=mean*(1-mean)
    known=.5*(pc*(1-pc)+pa*(1-pa))
    return {'attempts_per_source':budget,
            'clustered_probability':pc,'alternating_probability':pa,
            'expected_hybrid_seeds_clustered':c['expected_hybrid_seeds_produced'],
            'expected_hybrid_seeds_alternating':a['expected_hybrid_seeds_produced'],
            'coarse_brier':coarse,'network_known_brier':known,
            'network_information_gain':coarse-known,
            'exact_gain_identity':(pa-pc)**2/4}


def initial_one_step_colonization_bound(par=Network(), source=0, destination=1):
    """Optional structural bound on one occupied homogeneous donor and empty receiver.

    For this model one source requires a distinct local father for local
    outcross; with >=2 homozygous adults, exact local donor seed viability
    is s(1-d)+(1-s)*q*pollen[source][source].
    """
    if source not in range(NODES) or destination not in range(NODES) or source==destination:
        raise ValueError('need distinct source and destination')
    s=par.assurance*(1-par.q)
    viable=s*(1-par.inbreeding_depression)+(1-s)*par.q*par.pollen[source][source]
    u=viable*par.seed[source][destination]
    return 1-(1-u)**par.seed_attempts


def generation(populations, par, rng):
    """One synchronous life cycle. Every source has at most B seed attempts.

    Flow proceeds from OLD occupied parents; outcross requires a distinct
    local father or an occupied alternative paternal patch. Failed attempts
    are not retried. Each viable seed has one unique destination and is
    admitted at most once. Neutral finite recruitment chooses at most K seeds.
    """
    validate_populations(populations,par)
    pools=[[] for _ in range(NODES)]
    attempted=[0]*NODES
    viable=[0]*NODES
    exported=[0]*NODES
    hybrid_made=[0]*NODES
    s=par.assurance*(1-par.q)
    for i,mothers in enumerate(populations):
        if not mothers: continue
        for _ in range(par.seed_attempts):
            attempted[i]+=1
            mi=_index(len(mothers),rng); mom=mothers[mi]
            if rng.random()<s:
                if rng.random() >= 1-par.inbreeding_depression:
                    continue
                father_patch=i; dad=mom
            else:
                father_patch=_choose(par.pollen[i],rng)
                dads=populations[father_patch]
                if father_patch==i:
                    candidates=[g for j,g in enumerate(dads) if j!=mi]
                else:
                    candidates=dads
                if not candidates or rng.random()>=par.q:
                    continue
                dad=candidates[_index(len(candidates),rng)]
            seed=_cross(mom,dad,rng)
            viable[i]+=1
            if seed[0]>0 and seed[1]>0:
                hybrid_made[i]+=1
            d=_choose(par.seed[i],rng)
            exported[i]+=int(d!=i)
            pools[d].append(seed)
    if sum(viable)!=sum(len(pool) for pool in pools):
        raise AssertionError('a seed was duplicated or lost before recruitment')
    new=[]
    for pool in pools:
        if len(pool)>par.capacity:
            remaining=list(pool)
            recruited=[]
            for _ in range(par.capacity):
                recruited.append(remaining.pop(_index(len(remaining),rng)))
        else:
            recruited=pool
        if rng.random()<par.catastrophe:
            recruited=[]
        new.append(tuple(sorted(recruited)))
    ledger={
        'attempted':tuple(attempted),
        'viable_seed_produced':tuple(viable),
        'seed_exported':tuple(exported),
        'pre_recruit_seed_pool_sizes':tuple(len(pool) for pool in pools),
        'total_hybrid_seeds_made':sum(hybrid_made),
        'hybrid_seeds_made_by_source':tuple(hybrid_made),
    }
    return tuple(new),ledger


def _visible(populations,par,rng):
    """Separate assumed expression filter, never used in reproductive selection."""
    vals=[]
    for pop in populations:
        obs=[]
        for g in pop:
            pa=par.activation_A*par.q if g[0]>0 else 0.
            pb=par.activation_B*(1-par.q) if g[1]>0 else 0.
            a=rng.random()<pa; b=rng.random()<pb
            if a and b:
                obs.append('A' if rng.random()<.5 else 'B')
            elif a: obs.append('A')
            elif b: obs.append('B')
            else: obs.append('none')
        vals.append(tuple(obs))
    return tuple(vals)


def trajectory(initial=INIT,par=Network(),horizon=7,seed=73071027):
    validate_populations(initial,par)
    if type(horizon) is not int or not 1<=horizon<=30:
        raise ValueError('horizon must be integer 1..30')
    if type(seed) is not int: raise ValueError('seed must be integer')
    rng=Random(seed)
    pop=initial
    first_hybrid_produced=False; ever_joint=False; ever_copatch=False
    first_A_visible=False;first_B_visible=False
    attempted_all=0; viable_all=0; exported_all=0
    occupancy_integral=0
    first_generation_hybrid=False
    first_generation_joint_adults=False
    history=[]
    for step in range(1,horizon+1):
        pop,ledger=generation(pop,par,rng)
        vis=_visible(pop,par,rng)
        first_hybrid_produced |= ledger['total_hybrid_seeds_made']>0
        has_joint=any(g[0]>0 and g[1]>0 for patch in pop for g in patch)
        has_copatch=any(any(g[0]>0 for g in patch) and any(g[1]>0 for g in patch) for patch in pop)
        ever_joint |=has_joint; ever_copatch|=has_copatch
        first_A_visible |=any('A' in v for v in vis)
        first_B_visible |=any('B' in v for v in vis)
        attempted_all+=sum(ledger['attempted'])
        viable_all+=sum(ledger['viable_seed_produced'])
        exported_all+=sum(ledger['seed_exported'])
        occupancy_integral+=sum(bool(p) for p in pop)
        if step==1:
            first_generation_hybrid=ledger['total_hybrid_seeds_made']>0
            first_generation_joint_adults=has_joint
        history.append({
            'generation':step,
            'occupied_patches':sum(bool(x) for x in pop),
            'joint_genotype_present':has_joint,
            'same_patch_A_B_presence':has_copatch,
            'source_attempts':ledger['attempted'],
            'viable_seed_output':ledger['viable_seed_produced'],
            'source_exports':ledger['seed_exported'],
        })
    ga=sum(g[0] for p in pop for g in p)
    gb=sum(g[1] for p in pop for g in p)
    return {
        'horizon':horizon,
        'final_populations':pop,
        'gen1_hybrid_seed_produced':first_generation_hybrid,
        'gen1_joint_individual_recruited':first_generation_joint_adults,
        'ever_hybrid_seed_produced':first_hybrid_produced,
        'ever_joint_individual_recruited':ever_joint,
        'ever_cooccurrence_same_patch':ever_copatch,
        'ever_visible_A':first_A_visible,
        'ever_visible_B':first_B_visible,
        'final_joint_individual':any(g[0]>0 and g[1]>0 for p in pop for g in p),
        'final_same_patch_A_B':any(any(g[0] for g in p) and any(g[1] for g in p) for p in pop),
        'final_A_global':ga>0,
        'final_B_global':gb>0,
        'final_both_routes_global':ga>0 and gb>0,
        'final_all_patches_occupied':all(bool(p) for p in pop),
        'final_any_patch_occupied':any(bool(p) for p in pop),
        'final_occupied_patches':sum(bool(p) for p in pop),
        'time_integrated_occupied_patches':occupancy_integral,
        'total_attempts':attempted_all,
        'total_viable_seeds':viable_all,
        'total_exports':exported_all,
        'history':history,
    }


def _wilson(success,total,z=1.959963984540054):
    p=success/total;d=1+z*z/total;center=(p+z*z/(2*total))/d
    half=z*sqrt((p*(1-p)+z*z/(4*total))/total)/d
    return [max(0.,center-half),min(1.,center+half)]


def scenarios():
    """Predeclared network-only cells at identical B, q, K and founders."""
    cases={}
    for operator in ('seed','pollen','both'):
        for topology in ('components','clustered','alternating'):
            cases[f'{operator}_{topology}']=Network(
                seed=routing(topology) if operator in ('seed','both') else IDENTITY,
                pollen=routing(topology) if operator in ('pollen','both') else IDENTITY,
            )
    cases['seedcluster_pollenalt']=Network(seed=routing('clustered'), pollen=routing('alternating'))
    cases['seedalt_pollencluster']=Network(seed=routing('alternating'), pollen=routing('clustered'))
    cases['no_flow']=Network()
    return cases


def estimate(replicates=2500,horizon=7,seed=981027):
    """Paired-by-PRNG-seed prospective synthetic landscape Monte Carlo.

    Unit of independence is one whole six-patch replicate seed. Every cell
    reuses the same replicate keys; paired differences are estimated from
    paired 0/1 outcomes (not as independent-binomial contrasts). Shared seeds
    are not a guarantee of a low-variance coupling when routing differs.
    """
    if type(replicates) is not int or not 20<=replicates<=100000:
        raise ValueError('replicates must be integer 20..100000')
    if type(horizon) is not int or not 1<=horizon<=30:
        raise ValueError('horizon must be integer 1..30')
    if type(seed) is not int: raise ValueError('seed must be integer')
    cases=scenarios()
    bool_names=(
       'gen1_hybrid_seed_produced','gen1_joint_individual_recruited',
       'ever_hybrid_seed_produced','ever_joint_individual_recruited',
       'ever_cooccurrence_same_patch','ever_visible_A','ever_visible_B',
       'final_joint_individual','final_same_patch_A_B','final_A_global',
       'final_B_global','final_both_routes_global','final_all_patches_occupied',
       'final_any_patch_occupied',
    )
    scalar_names=('final_occupied_patches','time_integrated_occupied_patches',
                  'total_attempts','total_viable_seeds','total_exports')
    sums={name:{k:0. for k in bool_names+scalar_names} for name in cases}
    groups={
      'seed_alt_minus_cluster':('seed_clustered','seed_alternating'),
      'pollen_alt_minus_cluster':('pollen_clustered','pollen_alternating'),
      'both_alt_minus_cluster':('both_clustered','both_alternating'),
      'seed_cluster_minus_components':('seed_components','seed_clustered'),
      'pollen_cluster_minus_components':('pollen_components','pollen_clustered'),
      'both_cluster_minus_components':('both_components','both_clustered'),
      'pollen_components_minus_no_flow':('no_flow','pollen_components'),
      'pollen_topology_given_clustered_seed':('both_clustered','seedcluster_pollenalt'),
      'pollen_topology_given_alternating_seed':('seedalt_pollencluster','both_alternating'),
      'seed_topology_given_clustered_pollen':('both_clustered','seedalt_pollencluster'),
      'seed_topology_given_alternating_pollen':('seedcluster_pollenalt','both_alternating'),
    }
    endpoints=('gen1_hybrid_seed_produced','ever_joint_individual_recruited',
               'final_joint_individual','final_both_routes_global',
               'final_any_patch_occupied')
    contrast_sums={name:{k:[0.,0.] for k in endpoints} for name in groups}
    interaction_sums={k:[0.,0.] for k in endpoints}
    for i in range(replicates):
        key=seed+i*1000003
        vals={name:trajectory(par=par,horizon=horizon,seed=key)
              for name,par in cases.items()}
        for name,v in vals.items():
            for k in bool_names: sums[name][k]+=int(v[k])
            for k in scalar_names:
                val=float(v[k]);sums[name][k]+=val
        for name,(base,variant) in groups.items():
            for metric in endpoints:
                difference=int(vals[variant][metric])-int(vals[base][metric])
                contrast_sums[name][metric][0]+=difference
                contrast_sums[name][metric][1]+=difference*difference
        # Matched 2x2 topology-factorial difference-in-differences. This is a
        # policy interaction, not a partition of biological evolutionary rates.
        for metric in endpoints:
            v=(int(vals['both_alternating'][metric]) -
               int(vals['seedalt_pollencluster'][metric]) -
               int(vals['seedcluster_pollenalt'][metric]) +
               int(vals['both_clustered'][metric]))
            interaction_sums[metric][0]+=v
            interaction_sums[metric][1]+=v*v
    results={}
    for name,par in cases.items():
        results[name]={
            'n':replicates,
            'seed_topology':'identity' if par.seed==IDENTITY else next(k for k in ORDERS if par.seed==routing(k)),
            'pollen_topology':'identity' if par.pollen==IDENTITY else next(k for k in ORDERS if par.pollen==routing(k)),
            'probabilities':{k:sums[name][k]/replicates for k in bool_names},
            'wilson95':{k:_wilson(round(sums[name][k]),replicates) for k in bool_names},
            'means':{k:sums[name][k]/replicates for k in scalar_names},
        }
    contrasts={}
    for name,(base,variant) in groups.items():
        contrasts[name]={}
        for metric,(s1,s2) in contrast_sums[name].items():
            avg=s1/replicates
            sample_var=max(0.,(s2-replicates*avg*avg)/(replicates-1))
            se=sqrt(sample_var/replicates)
            contrasts[name][metric]={
                'baseline':base,'variant':variant,
                'difference_variant_minus_baseline':avg,
                'paired_mc_se':se,
                'paired_mc_95':[avg-1.95996398454*se,avg+1.95996398454*se],
                'sign_resolved':avg-1.95996398454*se>0 or avg+1.95996398454*se<0,
            }
    interaction={}
    for metric,(s1,s2) in interaction_sums.items():
        avg=s1/replicates
        sample_var=max(0.,(s2-replicates*avg*avg)/(replicates-1))
        se=sqrt(sample_var/replicates)
        interaction[metric]={'difference_in_differences':avg,'paired_mc_se':se,
             'paired_mc_95':[avg-1.95996398454*se,avg+1.95996398454*se]}
    exact={name:_first_hybrid_exact(par) for name,par in cases.items()}
    matrices={key:{'signature':mixing_signature(routing(key)),'matrix':[list(r) for r in routing(key)]}
              for key in ORDERS}
    return {
       'protocol':'SDEA v7: topology vs latent genotype alignment under matched finite per-source reproductive budgets',
       'scope':'CONSTRUCTED model; no natural dispersal/trait estimates; not independent source-repo reproduction',
       'design':{'replicates':replicates,'horizon':horizon,'unit':'whole six-patch landscape',
                 'replicate_seed':seed,'same_initial_genotypes':[[list(g) for g in p] for p in INIT],
                 'same_initial_q':.75,'same_seed_attempt_budget_per_occupied_source':3,
                 'same_capacity':2,'same_seed_retention_for_active_topologies':.25,
                 'same_pollen_retention_for_active_topologies':.25,
                 'outgoing_and_incoming_movement_totals_match_within_operator':True,
                 'future_total_attempts_not_controlled_because_source_occupancy_is_endogenous':True,
                 'paired_seed_keyed_model_draws':True},
       'network_matrices':matrices,'exact_first_generation_hybrid_production':exact,
       'exploratory_binary_saturation_and_brier':
           [network_prediction_risk(budget=b) for b in (1,2,3,8)],
       'scenarios':results,'paired_contrasts':contrasts,
       'factorial_topology_interaction':interaction,
    }
