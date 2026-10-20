"""E7 weight-lattice compression and minimal/adjoint image comparison."""
import unittest
from sage.all import ZZ, matrix
from sage_borel import e7_character_source


class E7AdjointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.minimal = e7_character_source()
        cls.adjoint = e7_character_source(representation="adjoint")

    def test_all_primes_adjoint_source(self):
        a, b = self.minimal.springer_map, self.adjoint.springer_map
        self.assertEqual(b.dimension, 133)
        self.assertEqual(b.trace_gram, 3*a.trace_gram)
        self.assertEqual(b.linear_basis_matrix, 36*matrix.identity(ZZ, 7))
        cert = self.adjoint.torus_certificate
        self.assertTrue(cert.certified)
        self.assertEqual(cert.exceptional_primes, ())
        self.assertEqual(cert.cases[0][2], (1,)*7)
        self.assertEqual(cert.cases[0][3], 7)

    def test_divisors_and_direct_laurent_checks(self):
        a, b = self.minimal.degree(1).image(), self.adjoint.degree(1).image()
        self.assertEqual(a.smith_factors, (24,))
        self.assertEqual(b.smith_factors, (36,))
        c = a.compare(b)
        self.assertFalse(c.left_in_right)
        self.assertFalse(c.right_in_left)
        self.assertEqual(c.sum.smith_factors, (12,))
        self.assertEqual(c.intersection.smith_factors, (72,))
        for d in range(3):
            row = self.adjoint.degree(d)
            for f in row.source_basis:
                self.assertTrue(self.adjoint.springer_map.is_integral(
                    self.adjoint.rational_algebra.substitute(f)))
            for m in row.monomials:
                algebra = self.adjoint.rational_algebra
                self.assertEqual(algebra.evaluate(m),
                                 algebra.borel_map.expand(algebra.substitute(m)))

    def test_coprime_image_formulas(self):
        for d in range(6):
            a, b = self.minimal.degree(d).image(), self.adjoint.degree(d).image()
            c = a.compare(b)
            T = c.intersection.basis.inverse()*(3**d*a.basis)
            self.assertEqual(T.denominator(), 1)
            self.assertEqual(abs(T.det()), 1)
            T = c.sum.basis.inverse()*(b.basis/(3**d))
            self.assertEqual(T.denominator(), 1)
            self.assertEqual(abs(T.det()), 1)
            self.assertEqual(c.sum.quotient_by(b).torsion_factors,
                             (3**d,)*a.rank if d else ())

    def test_selection_and_budget_rejections(self):
        explicit = e7_character_source(7, representation="minimal")
        self.assertEqual(explicit.springer_map.theta, self.minimal.springer_map.theta)
        for invalid in ("Adjoint", "133", 133, True, None):
            with self.assertRaises(ValueError):
                e7_character_source(representation=invalid)
        for kw in ({"max_monomials": 1}, {"max_matrix_entries": 1}):
            with self.assertRaises(ValueError):
                self.adjoint.degree(2, **kw)


if __name__ == "__main__":
    unittest.main()
