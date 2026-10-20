"""Explicit exceptional source forms and normalization-sensitive images."""
import unittest
from sage.all import QQ, RR
from sage_borel import exceptional_algebra, exceptional_borel


class ExceptionalIntegralTests(unittest.TestCase):
    def test_default_preserves_fractional_catalogue(self):
        algebra = exceptional_algebra("G2", 1)
        with self.assertRaises(ValueError):
            algebra.integral_degree(0)
        self.assertEqual(algebra.generator_polynomials[1],
                         algebra.borel_map.generator_polynomial("quadratic"))

    def test_g2_complete_images_and_explicit_rescaling(self):
        for node, scales in ((1, {"quadratic": 4}),
                             (1, {"quadratic": QQ(4)/3}), (2, {})):
            algebra = exceptional_algebra("G2", node, scales=scales)
            self.assertEqual([algebra.integral_degree(d).smith_factors
                              for d in range(7)],
                             [(1,), (1,), (1,), (2,), (2,), (2,), ()])

    def test_scaling_changes_f4_image(self):
        borel = exceptional_borel("F4", 1)
        original = borel.generator_algebra().integral_degree(4)
        primitive = borel.generator_algebra(scales={"e2": QQ(1)/2,
                                                   "e3": QQ(1)/4}).integral_degree(4)
        self.assertEqual(original.smith_factors, (2, 48))
        self.assertEqual(primitive.smith_factors, (1, 48))
        self.assertTrue(original.compare(primitive).left_in_right)
        self.assertFalse(original.compare(primitive).equal)
        self.assertEqual(borel.generator_image("e2").to_terms()[0]["coefficient"], "-10")

    def test_invalid_scaling(self):
        borel = exceptional_borel("G2", 1)
        for value in (0, -1, True, 0.5, RR(0.5)):
            with self.assertRaises(ValueError):
                borel.generator_algebra(scales={"quadratic": value})
        with self.assertRaises(ValueError):
            borel.generator_algebra(scales={"typo": 4})
        with self.assertRaises(TypeError):
            borel.generator_algebra(scales=[4])


if __name__ == "__main__":
    unittest.main()
