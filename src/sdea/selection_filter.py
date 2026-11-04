"""SDEA v8 — a selection-filter challenge to a reachability claim.

Synthetic, conditional six-patch diploid life-cycle. The pre-existing v7
Network, graph-routing and Mendelian kernels are reused, but post-zygotic
recruitment fitness is added as a *different* operator. The new operator
cannot affect seed production in the same generation: this is a deliberate
causal separation and a possible falsifier. No real genetic or fitness data.
"""
from collections import defaultdict
from dataclasses import dataclass
from math import isfinite, sqrt
from random import Random

from .network_topology import (
    A, B, INIT, NODES, Network, _choose, _cross, _first_hybrid_exact,
    _index, routing, validate_populations,
)


def _nonnegative_weight(value):
    if type(value) not in (int, float) or not isfinite(value) or not 0 <= value <= 30:
        raise ValueError("joint-zygote relative recruitment weight must be finite in [0,30]")


def joint(g):
    """Presence of at least one copy of each A/B route-associated allele."""
    return g[0] > 0 and g[1] > 0


def _weight(g, w):
    return w if joint(g) else 1.


def exact_pool_joint_recruitment(pool, capacity, joint_weight):
    """P(at least one joint genotype in a no-replacement weighted sample).

    The conditional seed pool is held fixed. A zero relative weight excludes
    all joint-genotype seeds, even if the pool has unfilled adult slots.
    The probability is enumerated over disjoint ordered-selection branches,
    not by multiplying independent recruit probabilities.
    """
    _nonnegative_weight(joint_weight)
    if type(capacity) is not int or capacity < 1:
        raise ValueError("capacity must be a positive integer")
    if type(pool) not in (tuple, list) or any(
        type(g) is not tuple or len(g) != 2 or any(type(x) is not int or not 0 <= x <= 2 for x in g)
        for g in pool
    ):
        raise ValueError("pool requires diploid A/B copy-count genotypes")
    admissible = tuple(g for g in pool if _weight(g, joint_weight) > 0.)
    if not admissible:
        return 0.
    if len(admissible) <= capacity:
        return float(any(joint(g) for g in admissible))

    def no_joint_draws(rem, k):
        if not rem or k == 0:
            return 1.
        denominator = sum(_weight(g, joint_weight) for g in rem)
        total = 0.
        for i, g in enumerate(rem):
            if not joint(g):
                total += (_weight(g, joint_weight) / denominator) * no_joint_draws(
                    rem[:i] + rem[i+1:], k-1
                )
        return total

    return 1. - no_joint_draws(admissible, capacity)


def exact_three_seed_witness(joint_weight):
    """A source-only exact witness: pool [A+B, A-only, B-only], K=2.

    Relative recruit weight (w,1,1). P(joint missed) =
    P(A then B)+P(B then A)=2/[(w+2)(w+1)].
    This is a conditional recruitment identity, not a population-level rate.
    """
    _nonnegative_weight(joint_weight)
    algebraic = 1. - 2. / ((joint_weight + 2.) * (joint_weight + 1.))
    enumerated = exact_pool_joint_recruitment(((1,1),A,B), 2, joint_weight)
    if abs(algebraic-enumerated)>1e-13:
        raise AssertionError("weighted-selection enumeration disagrees with analytic witness")
    return {'relative_weight':joint_weight, 'conditional_joint_recruitment':algebraic,
            'enumerated_cross_check':enumerated}


def production(populations, par, rng):
    """Viable diploid seeds formed in v7 order, BEFORE genotype selection.

    Inheritance follows v7's unlinked two-locus Mendelian gametes. Both
    parent location and seed destination are selected by the two matrices.
    Attempt counts and seed export budgets never depend on joint fitness.
    """
    validate_populations(populations, par)
    pools = [[] for _ in range(NODES)]
    attempts = [0]*NODES
    viable = [0]*NODES
    export = [0]*NODES
    made = [0]*NODES
    s = par.assurance*(1.-par.q)
    for i, mothers in enumerate(populations):
        if not mothers:
            continue
        for _ in range(par.seed_attempts):
            attempts[i]+=1
            mi=_index(len(mothers),rng)
            mom=mothers[mi]
            if rng.random()<s:
                if rng.random()>=1.-par.inbreeding_depression:
                    continue
                dad=mom
            else:
                father_patch=_choose(par.pollen[i],rng)
                dads=populations[father_patch]
                candidates=([g for j,g in enumerate(dads) if j!=mi]
                            if father_patch==i else list(dads))
                if not candidates or rng.random()>=par.q:
                    continue
                dad=candidates[_index(len(candidates),rng)]
            seed=_cross(mom,dad,rng)
            viable[i]+=1
            made[i]+=int(joint(seed))
            dest=_choose(par.seed[i],rng)
            export[i]+=int(dest!=i)
            pools[dest].append(seed)
    if sum(viable)!=sum(map(len,pools)) or sum(export)>sum(viable):
        raise AssertionError("per-generation seed production/assignment did not conserve")
    return tuple(tuple(pool) for pool in pools),{
        'source_attempts':tuple(attempts),
        'source_viable_seeds':tuple(viable),
        'source_exported_seeds':tuple(export),
        'source_joint_seeds':tuple(made),
        'destination_seed_counts':tuple(len(x) for x in pools),
    }


def settle(pool, capacity, joint_weight, rng):
    """Hard selection gate at w=0; positive w affects only crowding.

    This is *relative establishment competition*, not a source egg/seed
    viability modifier. Selection is therefore exactly downstream of seed
    production, and it is inactive if all admissible seeds fit available slots.
    """
    _nonnegative_weight(joint_weight)
    if type(capacity) is not int or capacity < 1:
        raise ValueError("capacity must be a positive integer")
    if any(type(g) is not tuple or len(g)!=2 or any(type(x) is not int or not 0<=x<=2 for x in g) for g in pool):
        raise ValueError("invalid diploid recruit genotype")
    remaining=[g for g in pool if _weight(g,joint_weight)>0]
    if len(remaining)<=capacity:
        return tuple(sorted(remaining))
    chosen=[]
    for _ in range(capacity):
        weights=[_weight(g,joint_weight) for g in remaining]
        total=sum(weights)
        i=_choose([x/total for x in weights],rng)
        chosen.append(remaining.pop(i))
    return tuple(sorted(chosen))


def generation_selected(populations, par, joint_weight, production_seed, recruitment_seed):
    """Synchronous generation; independent seeded production and selection stages.

    The independent streams provide paired same-production comparisons when
    input populations match. Across later generations initial conditions may
    differ *as a consequence* of earlier selection; no mediation fractions are
    inferred from the full trajectories.
    """
    _nonnegative_weight(joint_weight)
    pools, ledger=production(populations,par,Random(production_seed))
    rng=Random(recruitment_seed)
    offspring=[]
    recruited_hybrids=0
    for pool in pools:
        p=settle(pool,par.capacity,joint_weight,rng)
        if rng.random()<par.catastrophe:
            p=()
        recruited_hybrids+=sum(joint(g) for g in p)
        offspring.append(p)
    ledger={**ledger,
            'total_joint_seeds_produced':sum(ledger['source_joint_seeds']),
            'joint_adults_recruited':recruited_hybrids,
            'total_recruited_adults':sum(len(p) for p in offspring)}
    if ledger['total_recruited_adults']>sum(ledger['source_viable_seeds']):
        raise AssertionError("more recruits than unique produced seeds")
    return tuple(offspring),ledger


def trajectory(par, joint_weight, horizon=7, seed=8172026, initial=INIT):
    """Unconditional synthetic six-patch history, no mutation or immigration."""
    _nonnegative_weight(joint_weight)
    validate_populations(initial,par)
    if type(horizon) is not int or not 1<=horizon<=30:
        raise ValueError('horizon must be integer 1..30')
    if type(seed) is not int:
        raise ValueError('seed must be an integer')
    populations=initial
    first_hybrid_made=False
    first_joint_adult=False
    ever_joint=False
    ever_seed_made=False
    first_seed_gen=None
    viable_sum=0
    attempts_sum=0
    occupancy_integral=0
    first_gene_pool=0
    for t in range(1,horizon+1):
        # Explicit generation & stage keys: t1 production is common across w.
        populations,ledger=generation_selected(
            populations,par,joint_weight,seed+41009*t,seed+52009*t
        )
        this_made=ledger['total_joint_seeds_produced']>0
        this_recruited=ledger['joint_adults_recruited']>0
        if t==1:
            first_hybrid_made=this_made
            first_joint_adult=this_recruited
            first_gene_pool=ledger['total_joint_seeds_produced']
        if this_made and first_seed_gen is None:
            first_seed_gen=t
        ever_seed_made|=this_made
        ever_joint|=this_recruited
        viable_sum+=sum(ledger['source_viable_seeds'])
        attempts_sum+=sum(ledger['source_attempts'])
        occupancy_integral+=sum(bool(p) for p in populations)
    both_global=(any(g[0]>0 for p in populations for g in p) and
                 any(g[1]>0 for p in populations for g in p))
    return {
        'gen1_joint_seed_produced':first_hybrid_made,
        'gen1_joint_seed_count':first_gene_pool,
        'gen1_joint_adult_recruited':first_joint_adult,
        'ever_joint_seed_produced':ever_seed_made,
        'ever_joint_adult_recruited':ever_joint,
        'final_joint_adult':any(joint(g) for p in populations for g in p),
        'final_joint_adult_count':sum(joint(g) for p in populations for g in p),
        'final_both_alleles_global':both_global,
        'final_any_occupied':any(bool(p) for p in populations),
        'final_all_patches_occupied':all(bool(p) for p in populations),
        'final_occupied_patches':sum(bool(p) for p in populations),
        'total_seed_attempts':attempts_sum,
        'total_viable_seeds':viable_sum,
        'time_integrated_occupied_patches':occupancy_integral,
        'first_joint_seed_generation':first_seed_gen,
        'final_populations':populations,
    }


def _paired_summary(xs):
    """Paired MC point/normal precision interval; not natural-data CI."""
    n=len(xs)
    if n<2:
        raise ValueError('at least two independent landscape keys required')
    mean=sum(xs)/n
    var=sum((x-mean)**2 for x in xs)/(n-1)
    se=sqrt(var/n)
    return {'difference':mean,'paired_mc_se':se,
            'paired_mc_95':[mean-1.959963984540054*se, mean+1.959963984540054*se]}


def _wilson(k,n,z=1.959963984540054):
    p=k/n
    d=1+z*z/n
    c=(p+z*z/(2*n))/d
    e=z*sqrt((p*(1-p)+z*z/(4*n))/n)/d
    return [max(0.,c-e),min(1.,c+e)]


def cases():
    """Predefined 3 × 4 topology × competition table; no post-result cells."""
    result={}
    for topology in ('components','clustered','alternating'):
        mat=routing(topology)
        par=Network(seed=mat,pollen=mat)
        for key,w in [('lethal',0.),('disfavored',1./3.),('neutral',1.),('favored',3.)]:
            result[f'{topology}_{key}']=(par,w)
    return result


def benchmark(replicates=2000,horizon=7,seed=1182026):
    """Prespecified synthetic experiment; unit=whole six-patch replicate.

    Use paired keys for contrasts; stochastic readout is deterministic under
    Python 3.11 Random.random but NOT an independent empirical replication.
    """
    if type(replicates) is not int or not 50<=replicates<=100000:
        raise ValueError('replicates must be integer 50..100000')
    if type(horizon) is not int or not 1<=horizon<=30:
        raise ValueError('invalid horizon')
    if type(seed) is not int:
        raise ValueError('invalid seed')
    cells=cases()
    boolean_keys=(
        'gen1_joint_seed_produced','gen1_joint_adult_recruited',
        'ever_joint_seed_produced','ever_joint_adult_recruited',
        'final_joint_adult','final_both_alleles_global','final_any_occupied',
        'final_all_patches_occupied',
    )
    numeric_keys=('gen1_joint_seed_count','final_joint_adult_count',
                  'final_occupied_patches','total_seed_attempts',
                  'total_viable_seeds','time_integrated_occupied_patches')
    sums={name:{k:0. for k in (*boolean_keys,*numeric_keys)} for name in cells}
    hits={name:{k:0 for k in boolean_keys} for name in cells}
    comparisons={
        f'{topology}_favored_minus_neutral':(f'{topology}_neutral',f'{topology}_favored')
        for topology in ('clustered','alternating')
    }
    comparisons.update({
        f'{topology}_disfavored_minus_neutral':(f'{topology}_neutral',f'{topology}_disfavored')
        for topology in ('clustered','alternating')
    })
    comparisons.update({
        'neutral_alt_minus_cluster':('clustered_neutral','alternating_neutral'),
        'favored_alt_minus_cluster':('clustered_favored','alternating_favored'),
        'disfavored_alt_minus_cluster':('clustered_disfavored','alternating_disfavored'),
        'lethal_alt_minus_cluster':('clustered_lethal','alternating_lethal'),
    })
    endpoints=('gen1_joint_seed_produced','gen1_joint_adult_recruited',
               'ever_joint_adult_recruited','final_joint_adult',
               'final_any_occupied')
    differences={name:{key:[] for key in endpoints} for name in comparisons}
    for rep in range(replicates):
        trajectory_seed=seed+rep*1000003
        values={name:trajectory(par,w,horizon,trajectory_seed) for name,(par,w) in cells.items()}
        for name,v in values.items():
            for key in boolean_keys:
                hit=int(v[key]); hits[name][key]+=hit; sums[name][key]+=hit
            for key in numeric_keys:
                sums[name][key]+=v[key]
        for label,(base,variant) in comparisons.items():
            for key in endpoints:
                differences[label][key].append(int(values[variant][key])-int(values[base][key]))
    results={name:{
        'topology':name.rsplit('_',1)[0],
        'selection_policy':name.rsplit('_',1)[1],
        'joint_recruitment_weight':w,
        'replicate_unit':'entire six-patch landscape',
        'probabilities':{k:hits[name][k]/replicates for k in boolean_keys},
        'mc_wilson_95':{k:_wilson(hits[name][k],replicates) for k in boolean_keys},
        'means':{k:sums[name][k]/replicates for k in numeric_keys},
    } for name,(_,w) in cells.items()}
    theory={topology:_first_hybrid_exact(Network(pollen=routing(topology)))
            for topology in ('components','clustered','alternating')}
    witness={str(w):exact_three_seed_witness(w) for w in (0.,1./3.,1.,3.)}
    p0,p1=witness[str(1./3.)]['conditional_joint_recruitment'],witness[str(3.)]['conditional_joint_recruitment']
    brier_gain=(p1-p0)**2/4.
    return {
        'protocol':'SDEA v8 (2026-10-28) exact accessibility vs post-zygotic establishment selection',
        'scope':'SYNTHETIC; v7-based finite Mendelian generator; no natural trait regain or measured fitness',
        'design':{'replicates_per_cell':replicates,'horizon':horizon,'seed':seed,
                  'model_classes':['structurally_separated','connected_clustered','connected_alternating'],
                  'selection_policies':['lethal','disfavored','neutral','favored'],
                  'number_cells':len(cells),'independent_unit':'whole six-patch trajectory',
                  'fixed_q':.75,'fixed_seed_attempt_budget_per_occupied_source':3,
                  'fixed_initial_genotypes':[[list(g) for g in p] for p in INIT],
                  'selection_stage':'post-zygote competitive adult recruitment only',
                  'selection_affects_first_gen_seed_production':False,
                  'future_realized_total_seed_budget_endogenous':True,
                  'no_mutation_or_external_gene_reacquisition':True,
                  'paired_random_keys_not_proof_of_causal_mediation':True},
        'exact_first_gen_joint_seed_production':theory,
        'exact_fixed_mixed_pool_recruitment':witness,
        'exact_pool_forecast_brier_gain_favored_vs_disfavored':brier_gain,
        'scenarios':results,
        'paired_mc_contrasts':{label:{key:_paired_summary(arr) for key,arr in values.items()}
                               for label,values in differences.items()},
    }
