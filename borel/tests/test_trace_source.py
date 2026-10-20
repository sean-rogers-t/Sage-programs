"""Saturation and exponent-factorization hypotheses of the general theorem."""
from copy import copy
from itertools import combinations_with_replacement
import unittest

from sage.all import LaurentPolynomialRing, QQ, matrix
from sage_borel import (BorelMap, GeneratorAlgebra, IntegralCharacterSource,
                       SpringerCharacterMap, TraceSourceCertificate,
                       springer_character_map)


class TraceSourceTests(unittest.TestCase):
    def test_general_certificate_overcomes_inconclusive_leading_basis(self):
        s = springer_character_map("G2", (0, 1))
        self.assertFalse(s.certify_polynomial_source().certified)
        certificate = s.certify_trace_source()
        self.assertIsInstance(certificate, TraceSourceCertificate)
        self.assertTrue(certificate.certified)
        self.assertTrue(certificate.primitive_coefficients)
        self.assertTrue(certificate.trace_factorization)
        self.assertEqual(certificate.rank, 2)
        self.assertEqual(certificate.support_size, 12)

    def test_higher_representations_and_direct_sums(self):
        cases = [("A1", (5,)), ("A2", (2, 1)), ("G2", (2, 0)),
                 ("F4", (0, 0, 1, 0)), ("E6", (2, 0, 0, 0, 0, 0)),
                 ("E7", (0, 0, 0, 0, 0, 0, 2)),
                 ("E8", (0, 0, 0, 0, 0, 0, 0, 1))]
        for kind, highest in cases:
            self.assertTrue(springer_character_map(kind, highest).certify_trace_source().certified)
        for highest in ([(1, 0), (0, 1)], [(1, 0)]*3, [(2, 0), (0, 0)]):
            self.assertTrue(SpringerCharacterMap("G2", highest).certify_trace_source().certified)

    def test_primitive_coefficients_alone_are_insufficient(self):
        s = copy(springer_character_map("G2", (1, 0)))
        L = LaurentPolynomialRing(QQ, 2, names=("t1", "t2"))
        t, v = L.gens()
        s.linear_pullbacks = (t, t*t+2*v)
        certificate = s.certify_trace_source()
        self.assertTrue(certificate.primitive_coefficients)
        self.assertFalse(certificate.trace_factorization)
        self.assertFalse(certificate.certified)

    def test_trace_factorization_alone_is_insufficient(self):
        s = copy(springer_character_map("G2", (1, 0)))
        s.linear_basis_matrix = 2*s.linear_basis_matrix
        s.linear_pullbacks = tuple(2*f for f in s.linear_pullbacks)
        certificate = s.certify_trace_source()
        self.assertTrue(certificate.trace_factorization)
        self.assertFalse(certificate.primitive_coefficients)
        self.assertFalse(certificate.certified)

    def test_generic_source_uses_the_theorem_and_recovers_symmetric_square(self):
        s = springer_character_map("G2", (2, 0))
        b = BorelMap("G2")
        algebra = GeneratorAlgebra(b, dict(zip(("x", "y"), b.weights)))
        source = IntegralCharacterSource(s, algebra)
        self.assertIsInstance(source.torus_certificate, TraceSourceCertificate)
        row = source.degree(2)
        for f in row.source_basis:
            self.assertTrue(s.is_integral(algebra.substitute(f)))
        monomials = row.monomials
        linears = [sum(s.linear_basis_matrix[i, j]*algebra.domain.gen(i)
                       for i in range(2)) for j in range(2)]
        products = [a*c for a, c in combinations_with_replacement(linears, 2)]
        expected = matrix(QQ, [[f.monomial_coefficient(m) for f in products]
                               for m in monomials])
        transition = row.basis_matrix.inverse()*expected
        self.assertEqual(transition.denominator(), 1)
        self.assertEqual(abs(transition.det()), 1)


if __name__ == "__main__":
    unittest.main()
