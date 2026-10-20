"""Proof-sensitive checks for integral character sources and G2 image algebras."""
import unittest
from sage.all import CartanType, LaurentPolynomialRing, QQ, identity_matrix, matrix
from sage_borel import (SpringerCharacterMap, springer_character_map,
                       polynomial_source_certificate, modular_character_relation,
                       g2_character_algebra)


class CharacterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.g2 = springer_character_map("G2", (1, 0))

    def test_g2_explicit_trace_and_torus_formula(self):
        s = self.g2
        self.assertEqual(s.dimension, 7)
        self.assertEqual(s.trace_gram, matrix(QQ, [[12, -6], [-6, 4]]))
        self.assertEqual(s.linear_basis_matrix, matrix(QQ, [[6, 0], [0, 2]]))
        t, u = s.torus_ring.gens()
        b1 = t**2/u + 2*t - t/u + u/t - 2/t - u/t**2
        b2 = t - t/u + u/t - 1/t
        self.assertEqual(s.linear_pullbacks, (b1, b2))
        self.assertTrue(s.certify_polynomial_source().certified)
        self.assertEqual(s.certify_polynomial_source().exceptional_primes, ())
        self.assertTrue(s.is_integral(s.linear_generators[0]))
        self.assertFalse(s.is_integral(s.ring.gen(0)))
        self.assertEqual(matrix(QQ, 2, 2, lambda i, j:
                         s.theta[i].derivative(s.torus_ring.gen(j))(1, 1)),
                         identity_matrix(QQ, 2))

    def test_repeated_representation_keeps_normalized_map_and_source(self):
        twice = SpringerCharacterMap("G2", [(1, 0), (1, 0)])
        trivial_added = SpringerCharacterMap("G2", [(1, 0), (0, 0)])
        for s in (twice, trivial_added):
            self.assertEqual(s.theta, self.g2.theta)
            self.assertEqual(s.linear_basis_matrix, self.g2.linear_basis_matrix)
        self.assertEqual(twice.trace_gram, 2*self.g2.trace_gram)

    def test_adjoint_requires_different_leading_coordinates_and_prime_case(self):
        s = springer_character_map("G2", (0, 1))
        self.assertEqual(s.dimension, 14)
        self.assertEqual(s.linear_basis_matrix, 4*self.g2.linear_basis_matrix)
        self.assertFalse(s.certify_polynomial_source().certified)
        cert = s.certify_polynomial_source([[1, -3], [0, 1]])
        self.assertTrue(cert.certified)
        self.assertEqual(cert.exceptional_primes, (3,))
        self.assertEqual([c[3] for c in cert.cases], [2, 2])

    def test_rank_one_and_dual_representations(self):
        a1 = springer_character_map("A1", (1,))
        t = a1.torus_ring.gen()
        self.assertEqual(a1.theta, ((t-t**-1)/2,))
        self.assertTrue(a1.certify_polynomial_source().certified)
        a2 = springer_character_map("A2", (1, 0))
        dual = springer_character_map("A2", (0, 1))
        inverses = [t**-1 for t in a2.torus_ring.gens()]
        self.assertEqual(dual.theta, tuple(-f(*inverses) for f in a2.theta))

    def test_exponent_rank_is_over_Q_not_the_prime_field(self):
        L = LaurentPolynomialRing(QQ, 2, names=["t", "v"])
        t, v = L.gens()
        # These remain algebraically independent in characteristic 2,
        # despite their zero Jacobian and even exponent determinant.
        self.assertTrue(polynomial_source_certificate((t**2, v**2)).certified)
        self.assertIsNone(modular_character_relation((t**2, v**2), 2, 3))

    def test_true_higher_degree_modular_obstruction(self):
        L = LaurentPolynomialRing(QQ, 2, names=["t", "v"])
        t, v = L.gens()
        fs = (t, t*t+2*v)
        self.assertFalse(polynomial_source_certificate(fs).certified)
        self.assertIsNone(modular_character_relation(fs, 2, 1))
        relation = modular_character_relation(fs, 2, 2)
        self.assertIsNotNone(relation)
        self.assertEqual((relation(*fs)/2).monomial_coefficient(v), 1)
        # The linear coefficient lattice is primitive (unit minor from
        # t and t^2), but the full polynomial subring is not.
        self.assertTrue(any(c % 2 for c in relation.coefficients()))
        self.assertTrue(all(c.denominator() == 1 for c in (relation(*fs)/2).dict().values()))

    def test_no_obstruction_in_checked_g2_degrees(self):
        for p in (2, 3, 5):
            self.assertIsNone(self.g2.modular_obstruction(p, max_degree=3))

    def test_maximal_parabolic_images_and_integral_ring_relations(self):
        expected = {1: (1, 6, 12, 144, 288, 1728),
                    2: (1, 2, 12, 48, 288, 576)}
        for node in (1, 2):
            a = g2_character_algebra(node)
            for d, factor in enumerate(expected[node]):
                self.assertEqual(a.integral_degree(d).smith_factors, (factor,))
            self.assertEqual(a.integral_degree(6).smith_factors, ())
            linear, quadratic = a.generator_polynomials
            for f in a.generator_polynomials:
                self.assertTrue(self.g2.is_integral(f))
            if node == 1:
                self.assertFalse(a.borel_map.expand(linear**2-3*quadratic))
                self.assertFalse(a.borel_map.expand(quadratic**3))
            else:
                self.assertFalse(a.borel_map.expand(quadratic-3*linear**2))
                self.assertFalse(a.borel_map.expand(linear**6))

    def test_reject_invalid_inputs_and_explicit_resource_bounds(self):
        for weight in ((0, 0), (-1, 0), (1,), (True, 0), (0.5, 0)):
            with self.assertRaises(ValueError):
                springer_character_map("G2", weight)
        with self.assertRaises(ValueError):
            SpringerCharacterMap(["A", 2, 1], [(1, 0, 0)])
        with self.assertRaises(ValueError):
            springer_character_map(CartanType("G2").relabel({1: 2, 2: 1}), (1, 0))
        fs = self.g2.linear_pullbacks
        for p in (1, 4, True, 2.0):
            with self.assertRaises(ValueError):
                modular_character_relation(fs, p)
        for change in ([[2, 0], [0, 1]], [[1, 0]], [[QQ(1)/2, 0], [0, 2]]):
            with self.assertRaises(ValueError):
                polynomial_source_certificate(fs, change)
        for options in ({"max_monomials": 1}, {"max_terms": 1},
                        {"max_matrix_entries": 1}):
            with self.assertRaises(ValueError):
                modular_character_relation(fs, 2, 2, **options)
        with self.assertRaises(ValueError):
            polynomial_source_certificate((fs[0]/2, fs[1]))
        with self.assertRaises(TypeError):
            polynomial_source_certificate((self.g2.ring.gen(0),))
        for node in (0, 3, True):
            with self.assertRaises(ValueError):
                g2_character_algebra(node)


if __name__ == "__main__":
    unittest.main()
