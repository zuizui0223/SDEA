"""SDEA v4: exact multi-route, parentage operator and non-identifiability gates."""
import json
from dataclasses import replace
from math import isclose
from pathlib import Path
import unittest

from sdea.landscape import (
    Flow, GENOTYPES, ROUTES, counts, cross, tagged_offspring, reproduction_weights,
    route_copies, route_probabilities, all_count_states, _recruits,
    population_expression, one_generation, new_route_source,
    history_identifiability, evolve_recipient, benchmark,
)
from scripts.verify_fixture import assert_same

NONE, ANC, ALT, HET_A = (0,0), (2,0), (0,2), (1,0)
FOCAL_NONE = counts(NONE,NONE)
DONOR_ANC = counts(ANC,ANC)
FOCAL_ANC = counts(ANC,ANC)
DONOR_ALT = counts(ALT,ALT)
DONOR_MIX = counts(ANC,ALT)


def matched_pair():
    par=Flow(recipient_q=.6,assurance=.35,inbreeding_depression=.15,pollen_rate=.4)
    outcross = reproduction_weights(par.recipient_q,par.assurance,par.inbreeding_depression)["outcross"]
    seed=replace(par,pollen_rate=0.,seed_rate=outcross*.4)
    return par,seed


class TestLandscape(unittest.TestCase):
    def test_independent_mendelian_segregation(self):
        p=cross((1,1),(1,1))
        self.assertAlmostEqual(sum(p),1.)
        self.assertAlmostEqual(p[GENOTYPES.index((0,0))],.0625)
        self.assertAlmostEqual(p[GENOTYPES.index((1,1))],.25)
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(sum(cross((i,j),(i,j))),1.)

    def test_genotype_counts_and_combinatorics(self):
        self.assertEqual(len(GENOTYPES),9)
        self.assertEqual(len(all_count_states(2)),45)
        self.assertEqual(len(all_count_states(3)),165)
        self.assertEqual(route_copies(counts((1,2),(2,0)),0),3)
        self.assertEqual(route_copies(counts((1,2),(2,0)),1),2)

    def test_no_external_input_no_route_reopening(self):
        r=evolve_recipient(FOCAL_NONE,DONOR_ANC,Flow(),horizon=4)
        for x in r["trajectory"]:
            self.assertEqual(x["allele_supported_ancestral"],0.)
            self.assertEqual(x["first_route_by_horizon"]["ancestral"],0.)
            self.assertAlmostEqual(x["first_route_by_horizon"]["none"],1.)

    def test_pollen_delivers_fathers_not_whole_diploids(self):
        p=Flow(pollen_rate=1.,seed_rate=0.,assurance=0.)
        d=tagged_offspring(FOCAL_NONE,DONOR_ANC,p)
        self.assertAlmostEqual(d["origin"]["immigrant_pollen"][GENOTYPES.index(HET_A)],1.)
        self.assertAlmostEqual(d["genotype"][GENOTYPES.index(ANC)],0.)
        self.assertAlmostEqual(d["origin"]["immigrant_seed"][GENOTYPES.index(ANC)],0.)

    def test_seed_delivers_whole_zygote_genotypes(self):
        p=Flow(seed_rate=1.,pollen_rate=0.)
        d=tagged_offspring(FOCAL_NONE,DONOR_ANC,p)
        self.assertAlmostEqual(d["origin"]["immigrant_seed"][GENOTYPES.index(ANC)],1.)
        self.assertAlmostEqual(d["genotype"][GENOTYPES.index(HET_A)],0.)

    def test_donor_seed_formation_is_sexual_not_adult_cloning(self):
        donor=counts(HET_A,HET_A)
        p=Flow(seed_rate=1.,pollen_rate=0.)
        d=tagged_offspring(FOCAL_NONE,donor,p)["genotype"]
        # Two independent gametes from Aa yield Mendelian quarter-half-quarter.
        self.assertAlmostEqual(d[GENOTYPES.index(NONE)],.25)
        self.assertAlmostEqual(d[GENOTYPES.index(HET_A)],.50)
        self.assertAlmostEqual(d[GENOTYPES.index(ANC)],.25)

    def test_tagged_offspring_conserves_probability(self):
        p=Flow(pollen_rate=.4,seed_rate=.3,assurance=.5,recipient_q=.65)
        d=tagged_offspring(FOCAL_ANC,DONOR_ALT,p)
        total=sum(sum(v) for v in d["origin"].values())
        self.assertAlmostEqual(total,1.,places=12)
        for j in range(9):
            self.assertAlmostEqual(d["genotype"][j],sum(v[j] for v in d["origin"].values()))
        self.assertTrue(all(v>=0 for o in d["origin"].values() for v in o))

    def test_recruited_multinomial_kernel_normalizes(self):
        p=Flow(pollen_rate=.8,seed_rate=.2,assurance=.6)
        r=_recruits(FOCAL_ANC,DONOR_ALT,p)
        self.assertAlmostEqual(sum(prob for _,prob in r),1.,places=12)
        self.assertTrue(all(sum(c)==2 for c,_ in r))
        self.assertTrue(all(prob>=0 for _,prob in r))

    def test_matched_first_expression_with_different_genetic_dosage(self):
        pollen,seed=matched_pair()
        p=one_generation(FOCAL_NONE,DONOR_ANC,pollen)
        s=one_generation(FOCAL_NONE,DONOR_ANC,seed)
        self.assertAlmostEqual(p["first_generation_expression"]["ancestral"],
                               s["first_generation_expression"]["ancestral"],places=12)
        self.assertAlmostEqual(p["route_support_any"]["ancestral"],
                               s["route_support_any"]["ancestral"],places=12)
        self.assertAlmostEqual(s["mean_allele_copies"]["ancestral"],
                               2*p["mean_allele_copies"]["ancestral"],places=12)

    def test_independent_matched_expression_identity(self):
        pollen,seed=matched_pair()
        effective_out = reproduction_weights(pollen.recipient_q,pollen.assurance,
                                             pollen.inbreeding_depression)["outcross"]
        c = effective_out*pollen.pollen_rate
        expected_carrier = 1-(1-c)**pollen.recruits
        expected_visible = 1-(1-c*pollen.activation_ancestral*pollen.recipient_q)**pollen.recruits
        for v in (pollen,seed):
            r=one_generation(FOCAL_NONE,DONOR_ANC,v)
            self.assertAlmostEqual(r["first_generation_expression"]["ancestral"],expected_visible,places=12)
            self.assertAlmostEqual(r["route_support_any"]["ancestral"],expected_carrier,places=12)

    def test_same_visible_phenotype_cannot_identify_movement_source(self):
        pollen,seed=matched_pair()
        r=history_identifiability(FOCAL_NONE,DONOR_ANC,pollen,seed)
        self.assertEqual(r["visible_binary_total_variation"],0.)
        self.assertAlmostEqual(r["best_equal_prior_origin_accuracy_from_visible"],.5)
        self.assertGreater(r["full_genotype_count_total_variation"],.5)
        self.assertGreater(r["best_equal_prior_origin_accuracy_from_full_genotypes"],.75)

    def test_pulse_divergence_after_equal_first_generation(self):
        pollen,seed=matched_pair()
        p=evolve_recipient(FOCAL_NONE,DONOR_ANC,pollen,4,True)["final"]
        s=evolve_recipient(FOCAL_NONE,DONOR_ANC,seed,4,True)["final"]
        self.assertGreater(s["allele_supported_ancestral"],p["allele_supported_ancestral"])
        self.assertGreater(s["first_route_by_horizon"]["ancestral"],p["first_route_by_horizon"]["ancestral"])
        for key in ("mass","allele_supported_ancestral"):
            self.assertTrue(-1e-12<=p[key]<=1+1e-12 and -1e-12<=s[key]<=1+1e-12)

    def test_pollen_can_combine_routes_in_same_offspring(self):
        pollen=Flow(pollen_rate=.4,seed_rate=0.)
        seed=Flow(seed_rate=.4,pollen_rate=0.)
        a=one_generation(FOCAL_ANC,DONOR_ALT,pollen)
        b=one_generation(FOCAL_ANC,DONOR_ALT,seed)
        self.assertGreater(a["at_least_one_dual_competent"],0.)
        self.assertEqual(b["at_least_one_dual_competent"],0.)
        self.assertGreater(b["route_support_any"]["alternative"],0.)

    def test_exclusive_source_provenance_and_joint_movement(self):
        p=Flow(pollen_rate=.4,seed_rate=.4)
        r=new_route_source(FOCAL_ANC,DONOR_ALT,route=1,par=p)
        self.assertGreater(r["probabilities"]["both"],0.)
        self.assertGreater(r["probabilities"]["pollen_only"],0.)
        self.assertGreater(r["probabilities"]["seed_only"],0.)
        self.assertAlmostEqual(sum(r["probabilities"].values()),1.)
        self.assertAlmostEqual(r["any_route_reopened"],1-r["probabilities"]["none"])

    def test_provenance_requires_route_absent_in_recipient(self):
        with self.assertRaises(ValueError):
            new_route_source(FOCAL_ANC,DONOR_ALT,route=0)
        with self.assertRaises(ValueError):
            new_route_source(FOCAL_NONE,DONOR_ANC,route=2)

    def test_route_arbitration_is_exclusive_and_symmetric(self):
        p=Flow(recipient_q=.5,activation_ancestral=.8,activation_alternative=.8,tie_ancestral=.5)
        none,a,b=route_probabilities((1,1),.5,p)
        self.assertAlmostEqual(none+a+b,1.)
        self.assertAlmostEqual(a,b)
        a_bias=replace(p,tie_ancestral=1.)
        self.assertGreater(route_probabilities((1,1),.5,a_bias)[1],a)
        self.assertLess(route_probabilities((1,1),.5,a_bias)[2],b)

    def test_population_expression_categories_disjoint(self):
        p=Flow()
        s=counts((1,1),(1,0))
        x=population_expression(s,p.recipient_q,p)
        self.assertAlmostEqual(sum(x),1.,places=12)
        self.assertTrue(all(0<=z<=1 for z in x))
        self.assertGreater(x[3],0.)
        self.assertEqual(population_expression(FOCAL_NONE,p.recipient_q,p),(1.,0.,0.,0.))

    def test_balanced_source_both_routes_and_environment_filter(self):
        donor=DONOR_MIX
        pars=[Flow(recipient_q=q,pollen_rate=.4,seed_rate=.4) for q in (.2,.5,.8)]
        r=[one_generation(FOCAL_NONE,donor,x)["first_generation_expression"] for x in pars]
        self.assertGreater(r[0]["alternative"],r[0]["ancestral"])
        self.assertAlmostEqual(r[1]["alternative"],r[1]["ancestral"],places=12)
        self.assertGreater(r[2]["ancestral"],r[2]["alternative"])
        self.assertGreater(r[1]["simultaneous"],0.)

    def test_q_zero_blocks_pollen_transfer_but_not_seed_gene_input(self):
        p=Flow(recipient_q=0.,assurance=1.,inbreeding_depression=0.,pollen_rate=1.)
        s=replace(p,pollen_rate=0.,seed_rate=.4)
        a=one_generation(FOCAL_NONE,DONOR_ANC,p)
        b=one_generation(FOCAL_NONE,DONOR_ANC,s)
        self.assertEqual(a["route_support_any"]["ancestral"],0.)
        self.assertGreater(b["route_support_any"]["ancestral"],0.)
        self.assertEqual(b["first_generation_expression"]["ancestral"],0.)

    def test_dynamic_recursion_conserves_mass(self):
        p=Flow(pollen_rate=.4,seed_rate=.4,assurance=.7)
        for pulse in (True,False):
            r=evolve_recipient(FOCAL_NONE,DONOR_MIX,p,4,pulse)
            for x in r["trajectory"]:
                self.assertAlmostEqual(x["mass"],1.,places=11)
                self.assertAlmostEqual(sum(x["first_route_by_horizon"].values()),1.,places=11)
                self.assertTrue(all(0<=v<=1+1e-12 for v in x["first_route_by_horizon"].values()))

    def test_first_route_cumulative_and_simultaneous_are_distinct(self):
        p=Flow(pollen_rate=.4,seed_rate=.4,recipient_q=.4)
        r=evolve_recipient(FOCAL_NONE,DONOR_MIX,p,5,True)
        rows=r["trajectory"]
        for a,b in zip(rows,rows[1:]):
            self.assertLessEqual(a["first_route_by_horizon"]["ancestral"], b["first_route_by_horizon"]["ancestral"]+1e-12)
            self.assertLessEqual(a["first_route_by_horizon"]["alternative"], b["first_route_by_horizon"]["alternative"]+1e-12)
            self.assertGreaterEqual(a["first_route_by_horizon"]["none"]+1e-12,b["first_route_by_horizon"]["none"])
        self.assertGreater(r["final"]["first_route_by_horizon"]["simultaneous"],0.)

    def test_environment_history_affects_fate_given_prior_state(self):
        p=Flow(recipient_q=.2,pollen_rate=.4,seed_rate=.4)
        shift=evolve_recipient(FOCAL_NONE,DONOR_MIX,p,4,True,(.2,.2,.8,.8))["final"]
        persist=evolve_recipient(FOCAL_NONE,DONOR_MIX,p,4,True)["final"]
        self.assertGreater(shift["first_route_by_horizon"]["ancestral"],persist["first_route_by_horizon"]["ancestral"])
        self.assertLess(shift["first_route_by_horizon"]["none"],persist["first_route_by_horizon"]["none"])

    def test_invalid_scenarios_fail_closed(self):
        with self.assertRaises(ValueError):Flow(recruits=0)
        with self.assertRaises(ValueError):Flow(pollen_rate=1.1)
        with self.assertRaises(ValueError):Flow(recipient_q=0.,assurance=0.)
        with self.assertRaises(ValueError):Flow(donor_q=0.,assurance=0.)
        with self.assertRaises(ValueError):counts((0,0))
        with self.assertRaises(ValueError):cross((3,0),(0,0))
        with self.assertRaises(ValueError):tagged_offspring(FOCAL_NONE,counts((0,0)),Flow())
        with self.assertRaises(ValueError):evolve_recipient(FOCAL_NONE,DONOR_ANC,horizon=0)
        with self.assertRaises(ValueError):evolve_recipient(FOCAL_NONE,DONOR_ANC,horizon=3,q_schedule=(.2,.3))
        with self.assertRaises(ValueError):history_identifiability(FOCAL_NONE,DONOR_ANC,Flow(recruits=2),Flow(recruits=2,recipient_q=.5))
        with self.assertRaises(ValueError):history_identifiability(FOCAL_NONE,DONOR_ANC,Flow(assurance=.1),Flow(assurance=.9))

    def test_no_retroactive_change_v1_v2_v3_and_v4_fixture(self):
        path=Path(__file__).resolve().parents[1]/"results/v4_landscape_exact.json"
        frozen=json.loads(path.read_text())
        assert_same(frozen,benchmark())
        self.assertEqual(frozen["status"].split(';')[0],"SYNTHETIC")
        self.assertAlmostEqual(frozen["matched_seed_rate_for_carrier_equivalence"],
                               .4*reproduction_weights(.6,.35,.15)["outcross"])

if __name__ == "__main__":
    unittest.main()
