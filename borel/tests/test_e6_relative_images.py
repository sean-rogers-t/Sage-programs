"""E6 relative API versus independent polynomial and embedded lattice arithmetic."""
from itertools import combinations
import json
from pathlib import Path
import unittest
from sage.all import CartanType, PolynomialRing, WeylGroup, ZZ, matrix, diagonal_matrix
from sage_borel import (BorelMap, IntegralImageLattice, e6_character_source,
                        exceptional_root_weight_comparison)
from sage_borel.relative_images import ExceptionalRelativeDegree

P = PolynomialRing(ZZ, 'q')
q = P.gen()


def length_polynomial(nodes):
    out = P.one()
    if nodes:
        for degree in WeylGroup(CartanType('E6').subtype(nodes)).degrees():
            out *= sum(q**i for i in range(degree))
    return out


def exact_divide(a, b):
    quotient, remainder = P(a).quo_rem(P(b))
    assert not remainder
    return quotient


def independent_e6_polynomials():
    C = [q+q**2+q**7+q**8+q**9+2*q**10+q**11+q**12+q**14,
         q**3+q**4+q**6+q**7+q**9+q**11+q**12+q**13+q**15+q**16,
         q**5+q**8]
    D = P.one()
    for degree in (2,4,5,6,8):
        D *= sum(q**i for i in range(degree))
    A5 = [sum(q**d for d in range(10,16))+2*q**16+2*q**17+sum(q**d for d in range(18,22)),
          q**5+2*q**6+3*sum(q**d for d in range(7,15))+2*q**15+q**16+q**17+q**18,
          q**3+q**4]
    fixture = json.loads((Path(__file__).parent/'data/e6_relative_audited_base.json').read_text())
    A2A2A1 = [P(c) for c in fixture['coefficients']]
    special = {(1,3,4,5,6): A5, (1,2,3,5,6): A2A2A1,
               (1,3,5,6): [(1+q)*c for c in A2A2A1]}
    W = length_polynomial(tuple(range(1,7)))
    for size in range(6):
        for nodes in combinations(range(1,7), size):
            L = length_polynomial(nodes)
            H = exact_divide(W,L)
            counts = special[nodes] if nodes in special else [exact_divide(c*D,L) for c in C]
            yield nodes, H, counts


class E6RelativeImageTests(unittest.TestCase):
    def test_all_63_parabolics_and_all_prime_powers(self):
        for nodes, H, counts in independent_e6_polynomials():
            s = exceptional_root_weight_comparison('E6', retained_nodes=nodes)
            self.assertEqual(s.cohomology_ranks, tuple(int(H[d]) for d in range(H.degree()+1)))
            self.assertEqual(s.cartan_type, 'E6')
            for factor, polynomial in zip((3,9,27), counts):
                self.assertEqual(s.smith_multiplicity_polynomials[factor],
                                 tuple(int(polynomial[d]) for d in range(H.degree()+1)))
            for d in range(s.dimension+1):
                row = s.degree(d)
                expected = {factor: int(c[d]) for factor,c in zip((3,9,27),counts) if c[d]}
                self.assertEqual(row.nonunit_multiplicities, expected)
                self.assertEqual(row.relative_smith.rank, int(H[d]))
                self.assertEqual(s.multiplicity_coefficients[d], sum(expected.values()))
                comparison = row.compare_representations()
                self.assertEqual(comparison.left_scale, 6)
                self.assertEqual(comparison.right_scale, 24)
                self.assertEqual(comparison.sum_mod_left.torsion_multiplicities, expected)
                self.assertEqual(comparison.sum_mod_right.torsion_multiplicities,
                                 {4**d: int(H[d])} if d else {})

    def test_fresh_e6_p1_character_images(self):
        s = exceptional_root_weight_comparison('E6', retained_nodes=(2,3,4,5,6))
        minimal, adjoint = e6_character_source(), e6_character_source(representation='adjoint')
        for d in (0,1,2,3,5,8):
            actual = minimal.degree(d).image().compare(adjoint.degree(d).image())
            expected = s.compare_representations(degree=d)
            self.assertEqual(expected.left_in_right, actual.left_in_right)
            self.assertEqual(expected.right_in_left, actual.right_in_left)
            for group, actual_group in (
                (expected.sum_mod_left, actual.sum.quotient_by(actual.left)),
                (expected.sum_mod_right, actual.sum.quotient_by(actual.right)),
                (expected.left_mod_intersection, actual.left.quotient_by(actual.intersection)),
                (expected.right_mod_intersection, actual.right.quotient_by(actual.intersection))):
                self.assertEqual(group.torsion_factors(), actual_group.torsion_factors)

    def test_shared_three_primary_scales_with_integer_kernels(self):
        borel = BorelMap('A4')
        labels = tuple(borel.simple_reflections[i] for i in range(1,5))
        U = matrix(ZZ, [[1,2,0,3],[0,1,4,1],[0,0,1,2],[0,0,0,1]])
        V = U.transpose()
        row = ExceptionalRelativeDegree((),2,4,0,'test',cartan_type='E6',
                                         prime_power_multiplicities=((3,1),(9,1),(27,1)))
        for a,b in ((1,1),(1,3),(3,1),(3,9),(6,9),(9,6),(6,24),(24,6)):
            expected = row.compare_scaled_sources(root_scale=a, weight_scale=b)
            left = IntegralImageLattice._from_columns(1,labels,
                U*diagonal_matrix(ZZ,[a*a*c for c in (1,3,9,27)])*V,borel)
            right = IntegralImageLattice._from_columns(1,labels,b*b*U*V,borel)
            actual = left.compare(right)
            self.assertEqual(expected.left_in_right, actual.left_in_right)
            self.assertEqual(expected.right_in_left, actual.right_in_left)
            for group, actual_group in (
                (expected.sum_mod_left, actual.sum.quotient_by(left)),
                (expected.sum_mod_right, actual.sum.quotient_by(right)),
                (expected.left_mod_intersection, left.quotient_by(actual.intersection)),
                (expected.right_mod_intersection, right.quotient_by(actual.intersection))):
                self.assertEqual(group.torsion_factors(), actual_group.torsion_factors)

    def test_fresh_a5_full_source_kernels_in_degree_three(self):
        # Rebuild full degree-three root/weight invariants as integer kernels,
        # rather than using a rational invariant generator list or the tables.
        nodes = (1,3,4,5,6)
        borel = BorelMap('E6', omitted_nodes=(2,))
        R = PolynomialRing(ZZ,6,'u')
        u = R.gens()
        exponents = []
        def visit(prefix, remaining):
            if len(prefix) == 5:
                exponents.append(prefix+(remaining,))
            else:
                for a in range(remaining+1):
                    visit(prefix+(a,),remaining-a)
        visit((),3)
        def product(values,powers):
            out = values[0].parent().one()
            for value,a in zip(values,powers):
                out *= value**a
            return out
        monomials = [product(u,p) for p in exponents]
        w = borel.weights
        C = matrix(ZZ,6,6,lambda i,j:borel.simple_roots[i+1].monomial_coefficient(w[j]))
        images = []
        for root_basis in (True,False):
            equations = []
            for node in nodes:
                j = node-1
                if root_basis:
                    reflected = [u[k]-C[k,j]*u[j] for k in range(6)]
                else:
                    alpha = sum(C[j,k]*u[k] for k in range(6))
                    reflected = [u[k]-(alpha if k == j else 0) for k in range(6)]
                differences = [f(*reflected)-f for f in monomials]
                equations.extend([[f.monomial_coefficient(m) for f in differences] for m in monomials])
            K = matrix(ZZ,equations).right_kernel().basis_matrix()
            originals = [borel.simple_roots[i] for i in range(1,7)] if root_basis else w
            substituted = [product(originals,p) for p in exponents]
            labels = tuple(borel.basis(3))
            columns = []
            for row in K:
                expansion = borel.expand(sum(c*f for c,f in zip(row,substituted)))
                columns.append([expansion.get(label,0) for label in labels])
            images.append(IntegralImageLattice._from_columns(3,labels,matrix(ZZ,columns).transpose(),borel))
        M,N = images
        self.assertEqual(N.quotient_by(M).torsion_factors,(27,))
        expected = exceptional_root_weight_comparison('E6',retained_nodes=nodes).degree(3)
        self.assertEqual(expected.relative_smith.torsion_factors(),(27,))

    def test_boundaries_metadata_and_immutable_coefficients(self):
        f = exceptional_root_weight_comparison
        s = f('E6',retained_nodes=(6,5,4,3,1))
        self.assertEqual(s.degree(3).relative_smith.smith_factors(), (1,27))
        self.assertEqual(s.degree(3).as_dict()['cartan_type'], 'E6')
        self.assertIn('e6-',s.degree(3).as_dict()['certificate'])
        self.assertEqual(s.degree(3).multiplicity_2, 0)
        data = s.smith_multiplicity_polynomials
        data[27] = ()
        self.assertEqual(s.smith_multiplicity_polynomials[27][3],1)
        for nodes in ((7,), (1,1), (True,), (2.0,)):
            with self.assertRaises((ValueError,TypeError)):
                f('E6',retained_nodes=nodes)
        point = f('E6',retained_nodes=range(1,7))
        self.assertTrue(point.compare_representations(degree=0).equal)
        self.assertTrue(point.compare_representations(degree=10**12).equal)
        self.assertEqual(point.degree(1).relative_smith.rank,0)
        full = f('E6',retained_nodes=())
        self.assertEqual(full.schubert_rank,51840)
        self.assertEqual(full.dimension,36)
        json.dumps(full.as_dict())
        for rep in ('minimal','adjoint'):
            self.assertTrue(s.compare_representations(degree=3,left=rep,right=rep).equal)
        forward = s.compare_representations(degree=3)
        reverse = s.compare_representations(degree=3,left='adjoint',right='minimal')
        self.assertEqual(forward.sum_mod_left,reverse.sum_mod_right)


if __name__ == '__main__':
    unittest.main()
