"""Full source lattices, including combinations lost by individual rescaling."""
import unittest
from sage.all import QQ, matrix
from sage_borel import (BorelMap, GeneratorAlgebra, IntegralCharacterSource,
                       f4_character_source, springer_character_map)


class CharacterInvariantTests(unittest.TestCase):
    def test_f4_torus_and_low_degree_images(self):
        for node, linear, quadratic in ((1, 6, (36,)), (2, 6, (36, 72)),
                                         (3, 12, (72, 144)), (4, 12, (144,))):
            s = f4_character_source(node)
            self.assertEqual(s.springer_map.dimension, 26)
            self.assertEqual(s.springer_map.linear_basis_matrix,
                             matrix.diagonal(QQ, [6, 6, 12, 12]))
            self.assertTrue(s.torus_certificate.certified)
            self.assertEqual(s.torus_certificate.exceptional_primes, ())
            self.assertEqual(s.degree(1).image().smith_factors, (linear,))
            self.assertEqual(s.degree(2).image().smith_factors, quadratic)
            for d in range(4):
                row = s.degree(d)
                self.assertTrue(row.basis_matrix.is_immutable())
                for f in row.source_basis:
                    self.assertTrue(s.springer_map.is_integral(
                        s.rational_algebra.substitute(f)))

    def test_source_keeps_zero_image_and_mixed_integral_combinations(self):
        s = f4_character_source(1)
        a = s.rational_algebra
        e1, e2, e3, y = a.domain.gens()
        row = s.degree(2)
        self.assertEqual(len(row.source_basis), 2)
        self.assertEqual(row.image().rank, 1)
        f = 18*(e1-y**2)
        self.assertIn(f, row.source_basis)
        self.assertTrue(s.springer_map.is_integral(a.substitute(f)))
        self.assertTrue(s.springer_map.is_integral(a.substitute(18*(e1+y**2))))
        self.assertFalse(a.evaluate(18*(e1+y**2)))
        old = a.integral_degree(2)
        self.assertTrue(row.image().compare(old).left_in_right)
        self.assertFalse(row.image().compare(old).right_in_left)

    def test_nondiagonal_weight_lattice_and_exact_source_completeness(self):
        s = springer_character_map("A2", (1, 0))
        b = BorelMap("A2")
        a = GeneratorAlgebra(b, dict(zip(("x", "y"), b.weights)))
        source = IntegralCharacterSource(s, a, [[0, -1], [-1, -2]])
        self.assertNotEqual(s.linear_basis_matrix[0, 1], 0)
        # Equality of degree-one lattices tests the transpose convention.
        r = source.degree(1)
        actual = matrix(QQ, [[a.substitute(f).monomial_coefficient(w)
                              for f in r.source_basis] for w in b.weights])
        change = s.linear_basis_matrix.inverse()*actual
        self.assertEqual(change.denominator(), 1)
        self.assertEqual(abs(change.det()), 1)
        for d in (2, 3):
            for f in source.degree(d).source_basis:
                self.assertTrue(s.is_integral(a.substitute(f)))

    def test_subalgebra_scope_empty_degrees_and_rejections(self):
        s = springer_character_map("G2", (1, 0))
        b = BorelMap("G2")
        x, y = b.weights
        a = GeneratorAlgebra(b, {"square": x*x})
        source = IntegralCharacterSource(s, a)
        self.assertEqual(source.degree(1).source_basis, ())
        self.assertEqual(source.degree(1).image().rank, 0)
        self.assertEqual(source.degree(2).source_basis, (36*a.domain.gen(),))
        with self.assertRaises(ValueError):
            IntegralCharacterSource(s, GeneratorAlgebra(b, {"x": x, "xx": x*x}))
        with self.assertRaises(ValueError):
            IntegralCharacterSource(s, GeneratorAlgebra(BorelMap("A1"), {"x": BorelMap("A1").weights[0]}))
        with self.assertRaises(ValueError):
            IntegralCharacterSource(springer_character_map("F4", (0, 0, 0, 1)),
                                    f4_character_source(1).rational_algebra)
        for options in ({"max_monomials": 1}, {"max_matrix_entries": 1}):
            with self.assertRaises(ValueError):
                f4_character_source(2).degree(2, **options)
        for degree in (-1, True, 1.5):
            with self.assertRaises(ValueError):
                source.degree(degree)
        for node in (0, 5, True):
            with self.assertRaises(ValueError):
                f4_character_source(node)


if __name__ == "__main__":
    unittest.main()
