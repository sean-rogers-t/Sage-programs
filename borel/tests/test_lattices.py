"""Independent small lattice checks: Bezout, congruences, minors and geometry."""

from itertools import product
import json
import unittest

from sage.all import ZZ, gcd, matrix

from sage_borel import BorelMap, GeneratorAlgebra, compare_integral_images, two_step_symplectic
from tests.test_integral import determinantal_smith_factors


class LatticeComparisonTests(unittest.TestCase):
    def setUp(self):
        self.borel = BorelMap("A2")
        self.x, self.y = self.borel.weights

    def image(self, *polynomials):
        return GeneratorAlgebra(self.borel, {"g%d" % i: f for i, f in enumerate(polynomials)}).integral_degree(1)

    def test_bezout_and_non_coprime_intersection(self):
        for a, b in ((2, 3), (4, 6), (6, 15), (1, 7)):
            comparison = self.image(a * self.x).compare(self.image(b * self.x))
            self.assertEqual(comparison.sum.smith_factors, (gcd(a, b),))
            self.assertEqual(comparison.intersection.smith_factors, (a * b // gcd(a, b),))
            self.assertTrue(comparison.same_rational_span)
            self.assertEqual(comparison.left_in_right, a % b == 0)
            self.assertEqual(comparison.sum.quotient_by(comparison.left).index, a // gcd(a, b))

    def test_same_smith_distinct_lattices_and_support_alignment(self):
        comparison = self.image(2 * self.x, self.y).compare(self.image(self.x, 2 * self.y))
        self.assertEqual(comparison.left.smith_factors, comparison.right.smith_factors)
        self.assertFalse(comparison.equal)
        self.assertFalse(comparison.left_in_right)
        self.assertEqual(comparison.sum.smith_factors, (1, 1))
        self.assertEqual(comparison.intersection.smith_factors, (2, 2))
        self.assertEqual(comparison.left.quotient_by(comparison.intersection).torsion_factors, (2,))
        disjoint = self.image(self.x).compare(self.image(self.y))
        self.assertEqual(disjoint.intersection.rank, 0)
        self.assertFalse(disjoint.same_rational_span)
        self.assertEqual(disjoint.sum.quotient_by(disjoint.left).free_rank, 1)
        self.assertIsNone(disjoint.sum.quotient_by(disjoint.left).index)

    def test_skew_lattices_against_exhaustive_congruences_and_minors(self):
        # L = {(a*u+c*v,b*v)}. Membership is two explicit divisibility
        # conditions, independent of Smith transformations / rational solving.
        cases = ((2, 3, 1), (3, 2, -1), (2, 2, 1), (1, 4, 3))
        for p, q in product(cases, repeat=2):
            a, b, c = p
            d, e, f = q
            left = self.image(a*self.x, c*self.x+b*self.y)
            right = self.image(d*self.x, f*self.x+e*self.y)
            result = compare_integral_images(left, right)
            expected_sum = matrix(ZZ, [[a, c, d, f], [0, b, 0, e]])
            self.assertEqual(result.sum.smith_factors, determinantal_smith_factors(expected_sum))
            # Full rank determinant identity [Z^2:LcapM][Z^2:L+M]
            # = [Z^2:L][Z^2:M] checks intersection index independently.
            self.assertEqual(result.intersection.saturation_index * result.sum.saturation_index,
                             a*b*d*e)
            for u, v in product(range(-4, 5), repeat=2):
                in_left = v % b == 0 and (u - c*(v//b)) % a == 0
                in_right = v % e == 0 and (u - f*(v//e)) % d == 0
                value = self.borel.expand(u*self.x+v*self.y)
                self.assertEqual(result.intersection.contains(value), in_left and in_right)
            self.assertEqual(result.sum.quotient_by(result.left).as_dict(),
                             result.right.quotient_by(result.intersection).as_dict())
            self.assertEqual(result.sum.quotient_by(result.right).as_dict(),
                             result.left.quotient_by(result.intersection).as_dict())
            json.dumps(result.as_dict())

    def test_free_and_torsion_quotient(self):
        larger = self.image(self.x, self.y).lattice()
        smaller = self.image(6*self.x).lattice()
        quotient = larger.quotient_by(smaller)
        self.assertEqual(quotient.free_rank, 1)
        self.assertEqual(quotient.torsion_factors, (6,))
        self.assertIsNone(quotient.index)
        with self.assertRaises(ValueError):
            smaller.quotient_by(larger)

    def test_proper_subspaces_with_nonprimitive_overlap(self):
        borel = BorelMap("A3")
        x, y, z = borel.weights
        left = GeneratorAlgebra(borel, {"u": 2*x, "v": y}).integral_degree(1)
        right = GeneratorAlgebra(borel, {"u": 3*x, "v": z}).integral_degree(1)
        comparison = left.compare(right)
        self.assertEqual(comparison.sum.smith_factors, (1, 1, 1))
        self.assertEqual(comparison.intersection.smith_factors, (6,))
        self.assertTrue(comparison.intersection.contains(borel.expand(6*x)))
        self.assertFalse(comparison.intersection.contains(borel.expand(3*x)))
        quotient = comparison.sum.quotient_by(left)
        self.assertEqual((quotient.free_rank, quotient.torsion_factors), (1, (2,)))
        self.assertEqual(quotient.as_dict(), comparison.right.quotient_by(comparison.intersection).as_dict())

    def test_zero_lattices_and_degree_zero(self):
        algebra = GeneratorAlgebra(self.borel, {"x": self.x})
        zero = algebra.integral_degree(4)
        result = zero.compare(zero)
        self.assertTrue(result.equal)
        self.assertEqual(result.intersection.rank, 0)
        self.assertEqual(result.sum.quotient_by(zero).index, 1)
        self.assertTrue(result.sum.contains(self.borel.expand(0)))
        constant = algebra.integral_degree(0)
        self.assertEqual(constant.compare(constant).intersection.smith_factors, (1,))
        # A nonzero degree with no source monomials, compared to a nonzero image.
        even = GeneratorAlgebra(self.borel, {"q": self.x**2}).integral_degree(1)
        nonzero = self.image(self.x)
        result = even.compare(nonzero)
        self.assertTrue(result.left_in_right)
        self.assertFalse(result.right_in_left)
        self.assertEqual(result.sum.rank, 1)
        self.assertEqual(result.intersection.rank, 0)
        self.assertFalse(result.left.contains(self.borel.expand(self.x)))

    def test_context_type_degree_and_integrality_checks(self):
        report = self.image(self.x)
        with self.assertRaises(TypeError):
            report.compare(3)
        with self.assertRaises(ValueError):
            report.compare(GeneratorAlgebra(self.borel, {"x": self.x}).integral_degree(2))
        other = BorelMap("A2", omitted_nodes=1)
        with self.assertRaises(ValueError):
            report.compare(GeneratorAlgebra(other, {"x": other.weights[0]}).integral_degree(1))
        self.assertFalse(report.lattice().contains(self.borel.expand(self.x / 2)))
        self.assertFalse(report.lattice().contains(self.borel.expand(self.x + self.y)))
        with self.assertRaises(ValueError):
            report.lattice().contains(self.borel.expand(self.x**2))
        copy = BorelMap("A2")
        self.assertTrue(report.compare(GeneratorAlgebra(copy, {"x": copy.weights[0]}).integral_degree(1)).equal)
        self.assertTrue(report.lattice().basis.is_immutable())

    def test_joint_contains_mixed_products(self):
        flag = two_step_symplectic(3, 1, 1)
        middle = flag.middle_algebra().integral_degree(3)
        tail = flag.tail_algebra().integral_degree(3)
        joint = flag.joint_algebra().integral_degree(3)
        separate = middle.compare(tail).sum
        self.assertTrue(separate.compare(joint).left_in_right)
        self.assertGreater(joint.rank, separate.rank)
        self.assertGreater(joint.lattice().quotient_by(separate).free_rank, 0)


if __name__ == "__main__":
    unittest.main()
