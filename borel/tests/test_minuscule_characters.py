"""Certified E6/E7 sources, embedded comparisons and compressed coefficients."""
import unittest
from sage.all import QQ, ZZ, CartanMatrix, matrix
from sage_borel import e6_character_source, e7_character_source


class MinusculeCharacterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.e6 = e6_character_source()
        cls.adjoint = e6_character_source(representation="adjoint")
        cls.e7 = e7_character_source()

    def test_torus_sources_and_certificates(self):
        for s, scale, kind in ((self.e6, 6, "E6"), (self.e7, 12, "E7")):
            B = s.springer_map.linear_basis_matrix
            T = (scale*CartanMatrix(kind)).inverse()*B
            self.assertEqual(T.denominator(), 1)
            self.assertEqual(abs(T.det()), 1)
            self.assertTrue(s.torus_certificate.certified)
        s = self.adjoint.springer_map
        self.assertEqual(s.dimension, 78)
        self.assertEqual(s.linear_basis_matrix, 24*matrix.identity(ZZ, 6))
        self.assertEqual(s.trace_gram, 4*self.e6.springer_map.trace_gram)
        self.assertTrue(self.adjoint.torus_certificate.certified)
        self.assertEqual(self.adjoint.torus_certificate.exceptional_primes, ())

    def test_divisor_images_are_incomparable(self):
        a = self.e6.degree(1).image()
        b = self.adjoint.degree(1).image()
        self.assertEqual(a.smith_factors, (18,))
        self.assertEqual(b.smith_factors, (24,))
        c = a.compare(b)
        self.assertFalse(c.left_in_right)
        self.assertFalse(c.right_in_left)
        self.assertEqual(c.sum.smith_factors, (6,))
        self.assertEqual(c.intersection.smith_factors, (72,))
        self.assertEqual(c.sum.quotient_by(a).torsion_factors, (3,))
        self.assertEqual(c.sum.quotient_by(b).torsion_factors, (4,))

    def test_low_images_and_direct_bgg_columns(self):
        for source, factors in ((self.e6, (1, 18, 108, 1944)),
                                (self.adjoint, (1, 24, 576, 13824)),
                                (self.e7, (1, 24, 144, 3456))):
            algebra = source.rational_algebra
            for d, factor in enumerate(factors):
                row = source.degree(d)
                self.assertEqual(row.image().smith_factors, (factor,))
                self.assertTrue(row.basis_matrix.is_immutable())
                for m in row.monomials:
                    self.assertEqual(algebra.evaluate(m),
                                     algebra.borel_map.expand(algebra.substitute(m)))
                # Direct Laurent expansion is deliberately bounded.
                if d <= 2:
                    for f in row.source_basis:
                        self.assertTrue(source.springer_map.is_integral(algebra.substitute(f)))

    def test_source_kernel_and_budgets(self):
        row = self.e7.degree(2)
        self.assertEqual(len(row.source_basis), 2)
        self.assertEqual(row.image().rank, 1)
        for s in (self.e6, self.e7, self.adjoint):
            for kw in ({"max_monomials": 1}, {"max_matrix_entries": 1}):
                with self.assertRaises(ValueError):
                    s.degree(2, **kw)
            for d in (-1, True, QQ(1)/2):
                with self.assertRaises(ValueError):
                    s.degree(d)

    def test_supported_scope(self):
        for node in (0, 2, 6, True, "1", 1.0):
            with self.assertRaises(ValueError):
                e6_character_source(node)
        for node in (0, 1, 6, True, "7", 7.0):
            with self.assertRaises(ValueError):
                e7_character_source(node)
        for representation in ("Adjoint", "78", 78, None, True):
            with self.assertRaises(ValueError):
                e6_character_source(representation=representation)


if __name__ == "__main__":
    unittest.main()
