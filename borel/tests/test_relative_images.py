"""Independent polynomial and embedded-lattice checks of the E7 summary API.

The 112 formulas and A5 transfers are reconstructed using Sage polynomial
arithmetic, not the catalogue producer. D5+A1/D4+A1 use the coefficients
from the separately audited image certificates, not a fresh proof rerun.
"""
from collections import Counter
from dataclasses import FrozenInstanceError
from itertools import combinations
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from sage.all import CartanType, PolynomialRing, ZZ, WeylGroup, diagonal_matrix, matrix
from sage_borel import (BorelMap, IntegralImageLattice, e7_character_source,
                        exceptional_root_weight_comparison)
from sage_borel.relative_images import ExceptionalRelativeDegree

P = PolynomialRing(ZZ, 'q')
q = P.gen()


def bracket(n, step=1):
    return sum(q**(i*step) for i in range(n))


def poincare(degrees):
    result = P.one()
    for d in degrees:
        result *= bracket(int(d))
    return result


def exact_divide(a, b):
    quotient, remainder = P(a).quo_rem(P(b))
    assert not remainder
    return quotient


def gaussian(blocks, step=1):
    factorial = lambda n: poincare(range(2, n+1))
    denominator = P.one()
    for block in blocks:
        denominator *= factorial(block)
    return exact_divide(factorial(sum(blocks)), denominator)(q**step)


def independent_polynomials():
    """Construct all formulas without reading the runtime catalogue."""
    W = poincare(WeylGroup('E7').degrees())
    denominator = poincare([2,3,4,5,6])*poincare([2,3])
    P0 = exact_divide(W, denominator)
    T0 = P([0,0,0,0,1,0,2,2,2,4,5,4,7,7,8,9,12,8,16,12,15,13,17,
            10,18,12,15,10,15,8,13,8,8,6,8,3,5,3,2,2,2,0,1,0,0,0])
    U0 = exact_divide(q*P0, 1+q)
    T42 = U0*gaussian([4,2])-(U0-T0)*gaussian([2,1], step=4)
    T222 = U0*gaussian([2,2,2])-(U0-T0)*gaussian([1,1,1], step=4)
    fixture = json.loads((Path(__file__).parent/'data/e7_relative_audited_bases.json').read_text())
    special = {
        (2,5,7): T222*bracket(2)*bracket(3),
        (1,2,5,7): T222*bracket(3),
        (2,3,5,7): T222*bracket(3),
        (2,4,5,7): T42*bracket(2)*bracket(3),
        (2,5,6,7): T42*bracket(2)*bracket(3),
        (1,2,3,5,7): T222,
        (1,2,4,5,7): T42*bracket(3),
        (1,2,5,6,7): T42*bracket(3),
        (2,3,5,6,7): T42*bracket(3),
        (2,3,4,5,7): P(fixture['D4+A1']['coefficients']),
        (2,4,5,6,7): T0*bracket(2)*bracket(3),
        (1,2,3,4,5,7): P(fixture['D5+A1']['coefficients']),
        (1,2,3,5,6,7): T42,
        (1,2,4,5,6,7): T0*bracket(3),
        (2,3,4,5,6,7): q**8*(1+q)*(bracket(10,2)+q**4*bracket(6,2)+q**12),
    }
    for size in range(7):
        for nodes in combinations(range(1,8), size):
            # Sage determines the Levi root type independently of the producer's
            # custom graph classification. No E7 group elements are enumerated.
            L = poincare(WeylGroup(CartanType('E7').subtype(nodes)).degrees()) if nodes else P.one()
            H = exact_divide(W, L)
            T = special[nodes] if nodes in special else exact_divide(q*H, 1+q)
            yield nodes, H, T


class E7RelativeImageTests(unittest.TestCase):
    def test_all_parabolic_formulas_and_scaled_groups(self):
        gaps = Counter()
        for nodes, H, T in independent_polynomials():
            with self.subTest(nodes=nodes):
                result = exceptional_root_weight_comparison('E7', retained_nodes=nodes)
                self.assertEqual(result.cohomology_ranks, tuple(int(H[d]) for d in range(H.degree()+1)))
                self.assertEqual(result.multiplicity_coefficients,
                                 tuple(int(T[d]) for d in range(H.degree()+1)))
                gaps[next(d for d in range(H.degree()+1) if T[d])] += 1
                for d in range(result.dimension+1):
                    row = result.degree(d)
                    self.assertEqual(row.relative_smith.rank, int(H[d]))
                    self.assertEqual(row.relative_smith.torsion_multiplicities,
                                     {2: int(T[d])} if T[d] else {})
                    c = row.compare_representations()
                    self.assertEqual(c.sum_mod_left.torsion_multiplicities,
                                     row.relative_smith.torsion_multiplicities)
                    self.assertEqual(c.sum_mod_right.torsion_multiplicities,
                                     {3**d: int(H[d])} if d else {})
        self.assertEqual(gaps, {1:112, 2:12, 4:2, 8:1})

    def test_scaling_against_transformed_integer_lattices(self):
        # Deliberately lose diagonal form in a Schubert coordinate space, then
        # use integer kernels/Smith transforms for sums and intersections.
        borel = BorelMap('A3')
        indices = tuple(borel.simple_reflections[i] for i in range(1,4))
        U = matrix(ZZ, [[1,2,3],[0,1,4],[0,0,1]])
        V = matrix(ZZ, [[1,0,0],[2,1,0],[3,4,1]])
        for m in range(4):
            row = ExceptionalRelativeDegree((), 2, 3, m, 'test')
            for a, b in ((1,1),(2,1),(1,2),(6,10),(9,4),(12,36),(36,12)):
                compressed = row.compare_scaled_sources(root_scale=a, weight_scale=b)
                root = diagonal_matrix(ZZ, [a**2]*(3-m)+[2*a**2]*m)
                weight = b**2*matrix.identity(ZZ, 3)
                left = IntegralImageLattice._from_columns(1, indices, U*root*V, borel)
                right = IntegralImageLattice._from_columns(1, indices, U*weight*V, borel)
                actual = left.compare(right)
                self.assertEqual(compressed.left_in_right, actual.left_in_right)
                self.assertEqual(compressed.right_in_left, actual.right_in_left)
                for expected, actual_group in (
                    (compressed.sum_mod_left, actual.sum.quotient_by(left)),
                    (compressed.sum_mod_right, actual.sum.quotient_by(right)),
                    (compressed.left_mod_intersection, left.quotient_by(actual.intersection)),
                    (compressed.right_mod_intersection, right.quotient_by(actual.intersection))):
                    self.assertEqual(expected.torsion_factors(), actual_group.torsion_factors)

    def test_fresh_e7_p7_character_images(self):
        summary = exceptional_root_weight_comparison('E7', retained_nodes=range(1,7))
        minimal = e7_character_source()
        adjoint = e7_character_source(representation='adjoint')
        for degree in (0,1,2,3,6,10):
            actual = minimal.degree(degree).image().compare(adjoint.degree(degree).image())
            expected = summary.compare_representations(degree=degree)
            self.assertEqual(expected.left_in_right, actual.left_in_right)
            self.assertEqual(expected.right_in_left, actual.right_in_left)
            for group, lattice_group in (
                (expected.sum_mod_left, actual.sum.quotient_by(actual.left)),
                (expected.sum_mod_right, actual.sum.quotient_by(actual.right)),
                (expected.left_mod_intersection, actual.left.quotient_by(actual.intersection)),
                (expected.right_mod_intersection, actual.right.quotient_by(actual.intersection))):
                self.assertEqual(group.torsion_factors(), lattice_group.torsion_factors)

    def test_reversals_equal_sources_and_degrees(self):
        s = exceptional_root_weight_comparison('E7', retained_nodes=range(1,7))
        for d in (0,1,2,10,28,10**12):
            a = s.compare_representations(degree=d)
            b = s.compare_representations(degree=d, left='adjoint', right='minimal')
            self.assertEqual(a.sum_mod_left, b.sum_mod_right)
            self.assertEqual(a.left_mod_intersection, b.right_mod_intersection)
            for rep in ('minimal', 'adjoint'):
                same = s.compare_representations(degree=d, left=rep, right=rep)
                self.assertTrue(same.equal)
                self.assertTrue(same.sum_mod_left.is_trivial)
        self.assertEqual(s.degree(28).relative_smith.smith_factors(), ())
        self.assertEqual(s.degree(0).relative_smith.smith_factors(), (1,))
        point = exceptional_root_weight_comparison('E7', retained_nodes=range(1,8))
        self.assertEqual(point.schubert_rank, 1)
        self.assertTrue(point.compare_representations(degree=0).equal)
        self.assertTrue(point.compare_representations(degree=1).equal)

    def test_compressed_full_flag_and_guards(self):
        with patch('sage_borel.core.WeylGroup', side_effect=AssertionError('no matrix needed')):
            s = exceptional_root_weight_comparison('E7', retained_nodes=())
            self.assertEqual(s.schubert_rank, 2903040)
            middle = s.degree(31).relative_smith
            self.assertGreater(middle.rank, 100000)
            with self.assertRaises(ValueError):
                middle.smith_factors()
            with self.assertRaises(ValueError):
                middle.index(max_bits=10)
            self.assertLess(len(json.dumps(s.as_dict())), 5000)
        two = exceptional_root_weight_comparison('E7', retained_nodes=range(1,7)).degree(1).relative_smith
        self.assertEqual(two.index(), 2)
        self.assertEqual(two.torsion_factors(), (2,))
        self.assertEqual(two.multiplicities, {2:1})

    def test_validation_immutability_and_serialization(self):
        f = exceptional_root_weight_comparison
        for nodes in ((1,1), (0,), (8,), (True,), (1.0,), ('1',)):
            with self.assertRaises((TypeError, ValueError)):
                f('E7', retained_nodes=nodes)
        for kind in ('E5', 'E8', 'e7', None):
            with self.assertRaises(ValueError):
                f(kind, retained_nodes=())
        s = f('E7', retained_nodes=(7,5,2))
        self.assertEqual(s.retained_nodes, (2,5,7))
        for degree in (-1, True, 0.5):
            with self.assertRaises((TypeError, ValueError)):
                s.degree(degree)
        for scale in (0, -1, True, 1.5):
            with self.assertRaises((TypeError, ValueError)):
                s.compare_scaled_sources(degree=1, root_scale=scale)
        for rep in ('Minimal', 56, [], None):
            with self.assertRaises(ValueError):
                s.compare_representations(degree=1, left=rep)
        with self.assertRaises(FrozenInstanceError):
            s.dimension = 1
        data = s.as_dict()
        data['multiplicity_coefficients'][0] = 999
        self.assertEqual(s.multiplicity_coefficients[0], 0)
        json.dumps(s.compare_representations(degree=4).as_dict())


if __name__ == '__main__':
    unittest.main()
