"""Exact Mendelian single-locus operator and adverse inference controls."""
import json
from dataclasses import replace
from math import comb, isclose
from pathlib import Path
import unittest

from sdea.mendelian import (
    Genetics, allele_copies, all_genotype_counts, gamete_offspring_distribution,
    _next_counts, evolve_genotypes, selfing_structure_comparison, prediction_risks,
    single_generation_loss_identity, benchmark,
)


class TestMendelian(unittest.TestCase):
    def test_genotype_frequency_and_recruitment_normalization_exhaustive(self):
        p=Genetics(selection_strength=.7,direct_allele_cost=.4,assurance=.6,
                   mutation=.07, immigration=.03, inbreeding_depression=.25)
        for counts in all_genotype_counts(p.diploids):
            for qi in (0,1):
                survival, child, _, selffrac = gamete_offspring_distribution(counts, qi, p)
                self.assertTrue(0 <= survival <= 1)
                self.assertTrue(0 <= selffrac <= 1)
                self.assertAlmostEqual(sum(child),1.,places=12)
                self.assertAlmostEqual(sum(w for _,w in _next_counts(counts,qi,p)),1.,places=12)
                self.assertTrue(all(0 <= x <= 1 for x in child))

    def test_mendelian_two_gametes_from_one_heterozygote(self):
        # With all adults Aa and pure selfing at q=0, each viable zygote
        # inherits 1/4 aa, 1/2 Aa and 1/4 AA from Mendelian segregation.
        p=Genetics(q_low=0.,assurance=1.,inbreeding_depression=0.)
        success,child,_,eta=gamete_offspring_distribution((0,4,0),0,p)
        self.assertEqual(success,1.)
        self.assertEqual(child,(.25,.5,.25))
        self.assertEqual(eta,1.)

    def test_pure_selfing_homozygote_vs_heterozygote_proof(self):
        p=Genetics(q_low=0.,assurance=1.,inbreeding_depression=0.)
        a=gamete_offspring_distribution((2,2,0),0,p)[1]
        b=gamete_offspring_distribution((3,0,1),0,p)[1]
        self.assertAlmostEqual(a[0],5/8)
        self.assertAlmostEqual(b[0],3/4)
        self.assertAlmostEqual(a[0]**4,(5/8)**4)
        self.assertAlmostEqual(b[0]**4,(3/4)**4)
        self.assertGreater(b[0]**4,a[0]**4)

    def test_pure_outcross_matched_alleles_identical(self):
        p=Genetics(assurance=0.,selection_strength=0.,direct_allele_cost=0.)
        a=gamete_offspring_distribution((2,2,0),1,p)
        b=gamete_offspring_distribution((3,0,1),1,p)
        self.assertAlmostEqual(a[0],b[0])
        self.assertEqual(a[1],b[1])
        self.assertEqual(a[1],(.75**2,2*.25*.75,.25**2))
        f=selfing_structure_comparison(par=p,horizon=6)
        x=f['by_founder']['two_heterozygotes']['horizon']
        y=f['by_founder']['one_homozygote']['horizon']
        for endpoint in ('alive_and_ever_three','three_consecutive_survived_ever',
                         'alive_and_allele_supported'):
            self.assertAlmostEqual(x[endpoint],y[endpoint],places=13)

    def test_matched_copy_structure_differs_when_selfing_allowed(self):
        x=selfing_structure_comparison(par=Genetics(assurance=1.),horizon=6)
        self.assertEqual(x['initial_target_copies'],2)
        self.assertGreater(x['difference_absence_homo_minus_hetero'],0.)
        self.assertAlmostEqual(x['by_founder']['two_heterozygotes']['first_step']['offspring_target_absence'],
                               .1123038178349138,places=12)
        self.assertAlmostEqual(x['by_founder']['one_homozygote']['first_step']['offspring_target_absence'],
                               .13999033533432936,places=12)

    def test_one_generation_variance_decomposition_all_configurations(self):
        for par in (Genetics(),Genetics(selection_strength=1.,direct_allele_cost=.4,
                 mutation=.09,immigration=.07,inbreeding_depression=.32,assurance=.85)):
            for counts in all_genotype_counts(par.diploids):
                for qi in (0,1):
                    row=single_generation_loss_identity(counts,qi,par)
                    self.assertAlmostEqual(row['offspring_aa_by_identity'],
                                           row['offspring_aa_direct'],places=12)
                    self.assertAlmostEqual(row['conditional_complete_target_loss'],
                                           row['offspring_aa_direct']**par.diploids,places=12)

    def test_variance_term_zero_without_selfing(self):
        p=Genetics(assurance=0.)
        for founder in ((2,2,0),(3,0,1)):
            row=single_generation_loss_identity(founder,1,p)
            self.assertEqual(row['realized_viable_self_fraction'],0.)
            self.assertEqual(row['selfing_excess_non_target_zygotes'],0.)
            self.assertEqual(row['offspring_aa_by_identity'],(.75)**2)

    def test_neutral_expected_allele_copy_count_preserved(self):
        for assurance in (0.,.25,.75,1.):
            p=Genetics(assurance=assurance,mutation=0.,immigration=0.,
                       selection_strength=0.,direct_allele_cost=0.)
            for counts in all_genotype_counts(p.diploids):
                for qi in (0,1):
                    viability,child,_,_=gamete_offspring_distribution(counts,qi,p)
                    expected=2*p.diploids*(child[1]/2+child[2])
                    self.assertAlmostEqual(expected,allele_copies(counts),places=12)

    def test_cost_free_assurance_changes_conditional_support(self):
        lo=Genetics(assurance=0.,direct_allele_cost=0.,selection_strength=0.)
        hi=replace(lo,assurance=1.)
        founder=(3,0,1)
        s0,p0,_,_=gamete_offspring_distribution(founder,1,lo)
        s1,p1,_,_=gamete_offspring_distribution(founder,1,hi)
        self.assertGreater(s1,s0)
        self.assertGreater(p1[0]**hi.diploids,p0[0]**lo.diploids)
        self.assertAlmostEqual(1-p0[0]**lo.diploids,.8998870849609375)
        self.assertLess(1-p1[0]**hi.diploids,1-p0[0]**lo.diploids)

    def test_multigeneration_probability_conservation_with_feedback(self):
        par=Genetics(assurance=.7, selection_strength=.55,direct_allele_cost=.2,
                     mutation=.04,immigration=.02,expression_feedback=.15)
        x=evolve_genotypes((2,2,0),horizon=8,par=par)
        for row in x['trajectory']:
            self.assertAlmostEqual(row['alive']+row['extinct'],1.,places=10)
            self.assertGreaterEqual(row['first_expression_ever'],row['three_consecutive_survived_ever'])
            self.assertGreaterEqual(row['three_consecutive_survived_ever'],row['alive_and_ever_three'])
            self.assertLessEqual(row['alive_and_allele_supported'],row['alive']+1e-12)
            self.assertLessEqual(row['alive_and_fixed_target'],row['alive_and_allele_supported']+1e-12)

    def test_cumulative_events_monotonic_survival_nonincreasing(self):
        x=evolve_genotypes((2,2,0),horizon=10)
        for a,b in zip(x['trajectory'],x['trajectory'][1:]):
            self.assertLessEqual(a['first_expression_ever'],b['first_expression_ever']+1e-12)
            self.assertLessEqual(a['three_consecutive_survived_ever'],b['three_consecutive_survived_ever']+1e-12)
            self.assertGreaterEqual(a['alive']+1e-12,b['alive'])

    def test_no_target_copies_absorbing_without_input(self):
        p=Genetics()
        f=evolve_genotypes((4,0,0),horizon=6,par=p)['final']
        self.assertEqual(f['first_expression_ever'],0.)
        self.assertEqual(f['alive_and_allele_supported'],0.)
        self.assertGreater(f['alive'],0.)
        mutation=evolve_genotypes((4,0,0),horizon=6,par=replace(p,mutation=.02))['final']
        migration=evolve_genotypes((4,0,0),horizon=6,par=replace(p,immigration=.02))['final']
        self.assertGreater(mutation['alive_and_allele_supported'],0.)
        self.assertGreater(migration['alive_and_allele_supported'],0.)
        self.assertGreater(migration['first_expression_ever'],0.)

    def test_toolkit_absent_does_not_block_genetic_retention(self):
        yes=evolve_genotypes((2,2,0),horizon=6)['final']
        no=evolve_genotypes((2,2,0),horizon=6,toolkit=False)['final']
        self.assertEqual(no['first_expression_ever'],0.)
        self.assertEqual(no['alive_and_ever_three'],0.)
        self.assertAlmostEqual(yes['alive_and_allele_supported'],no['alive_and_allele_supported'],places=12)

    def test_three_run_requires_three_survived_updates(self):
        for H in (1,2):
            f=evolve_genotypes((2,2,0),horizon=H)['final']
            self.assertEqual(f['three_consecutive_survived_ever'],0.)
            self.assertGreater(f['first_expression_ever'],0.)

    def test_history_predictive_only_if_current_genotype_hidden(self):
        no_self=prediction_risks(par=Genetics(assurance=0.))
        selfing=prediction_risks(par=Genetics(assurance=1.))
        k='allele_loss_given_successful_reproduction'
        self.assertAlmostEqual(no_self[k]['history_gain'],0.,places=13)
        self.assertGreater(selfing[k]['history_gain'],0.)
        for endpoint in selfing:
            r=selfing[endpoint]
            self.assertAlmostEqual(r['full_genotype_plus_history_brier'],r['full_genotype_brier'],places=13)
            self.assertAlmostEqual(r['founder_history_brier'],r['full_genotype_brier'],places=13)

    def test_selection_can_change_direction_without_changing_initial_copy_total(self):
        # Do not universalize the neutral direction to an environmental
        # selection regime; verify a finite sign reversal witness.
        neutral=selfing_structure_comparison(par=Genetics(assurance=.7),horizon=6)
        selected=selfing_structure_comparison(par=Genetics(assurance=.7,selection_strength=.6),horizon=6)
        self.assertGreater(neutral['difference_absence_homo_minus_hetero'],0.)
        self.assertLess(selected['difference_absence_homo_minus_hetero'],0.)

    def test_positive_expression_feedback_increases_three_run_in_this_cell(self):
        p=Genetics(assurance=.7,selection_strength=.6)
        a=evolve_genotypes((2,2,0),horizon=6,par=p)['final']
        b=evolve_genotypes((2,2,0),horizon=6,par=replace(p,expression_feedback=.12))['final']
        self.assertGreater(b['three_consecutive_survived_ever'],a['three_consecutive_survived_ever'])

    def test_full_extinction_when_no_reproductive_route(self):
        par=Genetics(q_low=0.,assurance=0.)
        x=evolve_genotypes((2,2,0),start_q='low',horizon=4,par=par)['final']
        self.assertEqual(x['alive'],0.)
        self.assertEqual(x['extinct'],1.)

    def test_fails_closed_invalid_parameters(self):
        for bad in ({'diploids':0},{'diploids':True},{'q_low':.9},{'assurance':1.1},
                    {'inbreeding_depression':-0.1},{'consecutive_required':1},
                    {'direct_allele_cost':-1.}):
            with self.assertRaises(ValueError):Genetics(**bad)
        with self.assertRaises(ValueError):evolve_genotypes((4,0,1))
        with self.assertRaises(ValueError):evolve_genotypes([2,2,0])
        with self.assertRaises(ValueError):evolve_genotypes((2,2,0),horizon=0)
        with self.assertRaises(ValueError):evolve_genotypes((2,2,0),start_q='middle')
        with self.assertRaises(ValueError):evolve_genotypes((2,2,0),toolkit='yes')

    def test_exact_reproducibility_fixture(self):
        expected=json.loads((Path(__file__).resolve().parents[1]/'results/v3_mendelian_exact.json').read_text())
        actual=benchmark()
        def compare(a,b):
            self.assertIs(type(a),type(b))
            if isinstance(a,dict):
                self.assertEqual(set(a),set(b))
                for key in a:compare(a[key],b[key])
            elif isinstance(a,list):
                self.assertEqual(len(a),len(b))
                for x,y in zip(a,b):compare(x,y)
            elif isinstance(a,float):
                self.assertTrue(isclose(a,b,rel_tol=1e-12,abs_tol=1e-13),(a,b))
            else:self.assertEqual(a,b)
        compare(expected,actual)
        self.assertEqual(benchmark(),actual)


if __name__=='__main__':unittest.main()
