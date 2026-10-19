"""Independent arithmetic checks for selected integral Borel image lattices.

Run ``PYTHONPATH=src sage -python -m unittest tests.test_integral -v`` from
``borel/``.  The Smith invariants are checked by gcds of minors, never by a
second Smith-form call.  Small classical coefficient matrices come from the
independently reflected Cartesian BGG oracle.  Projective-space indices,
Bezout identities, and proper rational subspaces supply explicit geometric
checks.  Set ``BOREL_INTEGRAL_REPORT`` to save the exact check counts.
"""

from dataclasses import FrozenInstanceError
from itertools import combinations, product
import json
import os
from pathlib import Path
import time
import unittest
from unittest.mock import patch

from sage.all import QQ, ZZ, gcd, matrix, prod
from sage.env import SAGE_VERSION

from sage_borel import BorelMap, GeneratorAlgebra
from tests.oracle import OrthonormalBGG
from tests.test_products import ExplicitG2BGG


COUNTS = {}


def count(name, amount=1):
    COUNTS[name] = COUNTS.get(name, 0) + amount


def minor_gcd(coefficient_matrix, size):
    """Return the determinantal divisor, with Delta_0=1 by convention."""
    if size == 0:
        return ZZ.one()
    result = ZZ.zero()
    for rows in combinations(range(coefficient_matrix.nrows()), size):
        for columns in combinations(range(coefficient_matrix.ncols()), size):
            determinant = coefficient_matrix.matrix_from_rows_and_columns(rows, columns).det()
            result = gcd(result, determinant)
    return abs(result)


def determinantal_smith_factors(coefficient_matrix):
    previous = ZZ.one()
    factors = []
    for size in range(1, min(coefficient_matrix.dimensions()) + 1):
        divisor = minor_gcd(coefficient_matrix, size)
        if not divisor:
            break
        factors.append(ZZ(divisor // previous))
        previous = divisor
    return tuple(factors)


def weighted_exponents(weights, degree):
    return {exponents for exponents in product(
        *(range(degree // weight + 1) for weight in weights))
        if sum(e * weight for e, weight in zip(exponents, weights)) == degree}


class IntegralImageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.started = time.time()
        cls.methods_run = 0
        cls.oracles = {}

    def setUp(self):
        type(self).methods_run += 1
        type(self).test_result = self._outcome.result

    @classmethod
    def tearDownClass(cls):
        failed = any(isinstance(test, cls)
                     for test, _ in cls.test_result.failures + cls.test_result.errors)
        report = {
            "suite_result": "failed" if failed else "passed",
            "test_methods_run": cls.methods_run,
            "sage_version": SAGE_VERSION,
            "seconds": round(time.time() - cls.started, 3),
            "command": "PYTHONPATH=src sage -python -m unittest tests.test_integral -v",
            "independent_checks": COUNTS,
            "expected_values_source": (
                "Gcds of all minors; independent Cartesian BGG in A/C/D; "
                "explicit G2 roots; projective-space indices 2^d; Bezout "
                "primitive images; Gr(2,4) Schur classes and rational spans"),
        }
        print("\nIntegral image verification: " + json.dumps(report, sort_keys=True), flush=True)
        if os.environ.get("BOREL_INTEGRAL_REPORT"):
            Path(os.environ["BOREL_INTEGRAL_REPORT"]).write_text(
                json.dumps(report, indent=2) + "\n")

    def oracle(self, kind):
        if kind not in self.oracles:
            self.oracles[kind] = ExplicitG2BGG() if kind == "G2" else OrthonormalBGG(kind)
        return self.oracles[kind]

    def assert_preimage(self, algebra, report, expansion):
        self.assertTrue(report.contains(expansion))
        polynomial = report.preimage(expansion)
        self.assertEqual(polynomial.parent(), algebra.domain)
        self.assertTrue(all(coefficient in ZZ for coefficient in polynomial.coefficients()))
        self.assertTrue(all(sum(e * d for e, d in zip(exponents, algebra.generator_degrees))
                            == report.degree for exponents in polynomial.dict()))
        self.assertEqual(algebra.evaluate(polynomial), expansion)
        count("integer_preimages_reconstructed")
        return polynomial

    def assert_lattice_report(self, report):
        self.assertEqual(report.matrix.base_ring(), ZZ)
        self.assertTrue(report.matrix.is_immutable())
        expected = determinantal_smith_factors(report.matrix)
        self.assertEqual(report.smith_factors, expected)
        self.assertEqual(report.rank, len(expected))
        self.assertEqual(report.image_dimension, len(expected))
        self.assertEqual(report.torsion_factors, tuple(value for value in expected if value > 1))
        self.assertEqual(report.saturation_index, prod(expected))
        self.assertEqual(report.is_saturated, all(value == 1 for value in expected))
        image = report.image_lattice_basis
        saturation = report.saturation_lattice_basis
        self.assertEqual(image.dimensions(), (len(report.schubert_indices), report.rank))
        self.assertEqual(saturation.dimensions(), image.dimensions())
        for lattice in (image, saturation):
            self.assertEqual(lattice.base_ring(), ZZ)
            self.assertTrue(lattice.is_immutable())
            self.assertEqual(lattice.change_ring(QQ).column_space(),
                             report.matrix.change_ring(QQ).column_space())
        # A basis must generate the same lattice as the original columns.
        # Appending its columns cannot change the determinantal divisor.
        self.assertEqual(minor_gcd(image, report.rank), report.saturation_index)
        self.assertEqual(minor_gcd(report.matrix.augment(image), report.rank),
                         report.saturation_index)
        # Primitive columns spanning the same rational space give its full
        # intersection with the ambient integer Schubert lattice.
        self.assertEqual(minor_gcd(saturation, report.rank), 1)
        if report.rank:
            relative = saturation.change_ring(QQ).solve_right(image)
            self.assertTrue(all(value in ZZ for value in relative.list()))
            self.assertEqual(abs(relative.det()), report.saturation_index)
        serialized = report.as_dict()
        json.dumps(serialized)
        self.assertEqual(serialized["saturation_index"], report.saturation_index)
        self.assertEqual(serialized["smith_factors"], list(report.smith_factors))
        for name, lattice in (("matrix", report.matrix), ("image_lattice_basis", image),
                              ("saturation_lattice_basis", saturation)):
            encoded = serialized[name]
            reconstructed = matrix(ZZ, encoded["rows"], encoded["columns"],
                                   {(entry["row"], entry["column"]): ZZ(entry["coefficient"])
                                    for entry in encoded["entries"]})
            self.assertEqual(reconstructed, lattice)
            count("exact_integer_matrix_serializations")
        count("gcd_of_minors_smith_and_lattice_checks")

    def assert_independent_matrix(self, algebra, degree, kind, coordinate_generators):
        report = algebra.integral_degree(degree)
        oracle = self.oracle(kind)
        exponents = tuple(tuple(int(e) for e in monomial.exponents()[0])
                          for monomial in report.monomials)
        self.assertEqual(set(exponents), weighted_exponents(algebra.generator_degrees, degree))
        self.assertEqual(report.domain_dimension, len(exponents))
        columns = []
        expected_support = set()
        for powers in exponents:
            polynomial = prod(f ** e for f, e in zip(coordinate_generators, powers))
            expansion = oracle.expansion(polynomial, degree)
            local = {algebra.borel_map.W.from_reduced_word(w.reduced_word()): ZZ(c)
                     for w, c in expansion.items()}
            expected_support.update(local)
            columns.append(local)
        self.assertEqual(set(report.schubert_indices), expected_support)
        expected_matrix = matrix(ZZ, len(report.schubert_indices), len(columns),
                                 {(row, column): expansion.get(word, 0)
                                  for row, word in enumerate(report.schubert_indices)
                                  for column, expansion in enumerate(columns)})
        self.assertEqual(report.matrix, expected_matrix)
        self.assert_lattice_report(report)
        count("independent_cartesian_image_matrices_" + kind)
        return report

    def test_projective_space_primitive_and_doubled_generator(self):
        borel = BorelMap("A4", omitted_nodes=1)
        h = borel.weights[0]
        for scale in (1, 2):
            algebra = GeneratorAlgebra(borel, {"u": scale * h})
            u, = algebra.domain.gens()
            for degree in range(5):
                report = algebra.integral_degree(degree, include_ambient=True)
                self.assertEqual(report.smith_factors, (ZZ(scale) ** degree,))
                self.assertEqual(report.saturation_index, scale ** degree)
                self.assertEqual(report.ambient_dimension, 1)
                self.assertEqual(report.free_cokernel_rank, 0)
                self.assertTrue(report.full_ambient)
                self.assert_lattice_report(report)
                target = borel.expand((scale * h) ** degree)
                preimage = self.assert_preimage(algebra, report, target)
                self.assertEqual(preimage, u ** degree)
                if scale == 2 and degree:
                    primitive = borel.expand(h ** degree)
                    self.assertFalse(report.contains(primitive))
                    with self.assertRaises(ValueError):
                        report.preimage(primitive)
                    count("projective_nonprimitive_membership_rejections")
                count("projective_known_indices")

    def test_redundant_generators_give_bezout_primitive_images(self):
        borel = BorelMap("A3", omitted_nodes=1)
        h = borel.weights[0]
        algebra = GeneratorAlgebra(borel, {"u": 2 * h, "v": 3 * h})
        for degree in range(4):
            report = algebra.integral_degree(degree)
            self.assertEqual(report.smith_factors, (1,))
            self.assertEqual(report.saturation_index, 1)
            self.assert_lattice_report(report)
            preimage = self.assert_preimage(algebra, report, borel.expand(h ** degree))
            if degree == 1:
                self.assertTrue(any(coefficient < 0 for coefficient in preimage.coefficients()))
            count("bezout_primitive_degree_images")

    def test_classical_matrices_and_smith_factors_by_minors(self):
        for kind in ("A2", "C2"):
            borel, oracle = BorelMap(kind), self.oracle(kind)
            u, v = borel.weights
            algebra = GeneratorAlgebra(borel, {"u": 2 * u, "v": 3 * v})
            inputs = (2 * oracle.omega(1), 3 * oracle.omega(2))
            for degree in range(borel.longest_representative().length() + 1):
                report = self.assert_independent_matrix(algebra, degree, kind, inputs)
                if degree == 1:
                    self.assertEqual(report.smith_factors, (1, 6))
        borel, oracle = BorelMap("A2"), self.oracle("A2")
        algebra = GeneratorAlgebra(borel, {"u": 2 * borel.weights[0],
                                         "v": 6 * borel.weights[1]})
        report = self.assert_independent_matrix(
            algebra, 1, "A2", (2 * oracle.omega(1), 6 * oracle.omega(2)))
        self.assertEqual(report.smith_factors, (2, 6))
        self.assertEqual(report.torsion_factors, (2, 6))

    def test_type_c_tail_and_type_d_euler_lattices(self):
        borel, oracle = BorelMap("C3", omitted_nodes=1), self.oracle("C3")
        algebra = GeneratorAlgebra(borel, {"q1": borel.tail_polynomial(2, 1),
                                         "q2": borel.tail_polynomial(2, 2)})
        x2, x3 = oracle.x[1:]
        for degree in (0, 2, 4, 5):
            report = self.assert_independent_matrix(
                algebra, degree, "C3", (x2 ** 2 + x3 ** 2, x2 ** 2 * x3 ** 2))
            if degree in (2, 4):
                self.assertEqual(report.smith_factors, (1,))
        # IG(1,6)=P5: the corrected tail has q1=-h^2 and q2=+h^4.
        self.assertEqual(borel.tail_elementary(2, 1), -borel.expand(borel.weights[0] ** 2))
        self.assertEqual(borel.tail_elementary(2, 2), borel.expand(borel.weights[0] ** 4))
        count("C3_original_tail_sign_checks", 2)
        borel, oracle = BorelMap("D4", omitted_nodes=2), self.oracle("D4")
        algebra = GeneratorAlgebra(borel, {"q": borel.tail_polynomial(2, 1),
                                         "euler": borel.euler_polynomial(2)})
        x3, x4 = oracle.x[2:]
        for degree in (0, 2, 4):
            report = self.assert_independent_matrix(algebra, degree, "D4",
                                                    (x3 ** 2 + x4 ** 2, x3 * x4))
            for column in report.image_lattice_basis.columns():
                expansion = sum((coefficient * borel.schubert(word)
                                 for word, coefficient in zip(report.schubert_indices, column)),
                                borel.expand(0))
                self.assert_preimage(algebra, report, expansion)
            if degree == 2:
                self.assertEqual(report.smith_factors, (1, 1))
            if degree == 4:
                self.assertEqual(report.smith_factors, (1, 2, 2))
                self.assertEqual(report.saturation_index, 4)
                missing = 0
                for column in report.saturation_lattice_basis.columns():
                    expansion = sum((coefficient * borel.schubert(word)
                                     for word, coefficient in zip(report.schubert_indices, column)),
                                    borel.expand(0))
                    if not report.contains(expansion):
                        with self.assertRaises(ValueError):
                            report.preimage(expansion)
                        missing += 1
                self.assertGreater(missing, 0)
                count("D4_saturation_classes_missing_from_image", missing)
        degree_two = algebra.integral_degree(2)
        self.assert_preimage(algebra, degree_two, borel.euler(2))

    def test_type_c_fixed_character_normalization_changes_integral_index(self):
        borel, oracle = BorelMap("C3", omitted_nodes=2), self.oracle("C3")
        q = borel.tail_polynomial(1, 1)
        # For the one-variable symplectic tail, the literal polynomial
        # character normalization is 4q.  Retaining that normalization must
        # change the integer lattice, even though its rational span agrees.
        primitive = GeneratorAlgebra(borel, {"q": q})
        character_normalized = GeneratorAlgebra(borel, {"z": 4 * q})
        for degree in (2, 4):
            exponent = degree // 2
            small = self.assert_independent_matrix(primitive, degree, "C3", (oracle.x[-1] ** 2,))
            scaled = self.assert_independent_matrix(
                character_normalized, degree, "C3", (4 * oracle.x[-1] ** 2,))
            self.assertEqual(small.smith_factors, (1,))
            self.assertEqual(scaled.smith_factors, (4 ** exponent,))
            self.assertEqual(small.rank, scaled.rank)
            self.assert_preimage(primitive, small, borel.expand(q ** exponent))
            self.assert_preimage(character_normalized, scaled,
                                 borel.expand((4 * q) ** exponent))
            self.assertFalse(scaled.contains(borel.expand(q ** exponent)))
            with self.assertRaises(ValueError):
                scaled.preimage(borel.expand(q ** exponent))
            count("fixed_type_C_character_normalization_indices")

    def test_proper_rational_span_and_missing_support(self):
        borel = BorelMap("A3", omitted_nodes=2)
        w1, h, _ = borel.weights
        c2 = w1 * (h - w1)
        s2, s11 = borel.expand(h ** 2 - c2), borel.expand(c2)
        diagonal = s2 + s11
        algebra = GeneratorAlgebra(borel, {"u": 2 * h ** 2})
        report = algebra.integral_degree(2, include_ambient=True)
        self.assertEqual(report.rank, 1)
        self.assertEqual(report.ambient_dimension, 2)
        self.assertEqual(report.free_cokernel_rank, 1)
        self.assertFalse(report.full_ambient)
        self.assertEqual(report.smith_factors, (2,))
        self.assert_lattice_report(report)
        self.assert_preimage(algebra, report, 2 * diagonal)
        for outside in (diagonal, s2, s11, QQ(1) / 2 * diagonal):
            self.assertFalse(report.contains(outside))
            with self.assertRaises(ValueError):
                report.preimage(outside)
            count("proper_span_or_nonintegral_membership_rejections")
        primitive = GeneratorAlgebra(borel, {"u": h ** 2}).integral_degree(
            2, include_ambient=True)
        self.assertTrue(primitive.is_saturated)
        self.assertEqual(primitive.free_cokernel_rank, 1)
        self.assertFalse(primitive.full_ambient)
        supported = GeneratorAlgebra(borel, {"s2": h ** 2 - c2}).integral_degree(2)
        self.assertEqual(len(supported.schubert_indices), 1)
        self.assertFalse(supported.contains(s11))
        with self.assertRaises(ValueError):
            supported.preimage(s11)
        count("outside_stored_support_rejections")

    def test_rational_polynomials_with_integral_images_are_accepted(self):
        borel = BorelMap("A2")
        u, v = borel.weights
        invariant = u ** 2 - u * v + v ** 2
        polynomial = u ** 2 + invariant / 2
        self.assertTrue(any(coefficient not in ZZ for coefficient in polynomial.coefficients()))
        algebra = GeneratorAlgebra(borel, {"q": polynomial})
        report = algebra.integral_degree(2)
        self.assertEqual(report.smith_factors, (1,))
        self.assert_preimage(algebra, report, borel.expand(u ** 2))
        self.assert_lattice_report(report)
        line = BorelMap("A1")
        zero = GeneratorAlgebra(line, {"z": line.weights[0] ** 2 / 2})
        self.assertEqual(zero.integral_degree(0).smith_factors, (1,))
        self.assertEqual(zero.integral_degree(2).rank, 0)
        count("fractional_input_integral_image_acceptance", 3)

    def test_nonintegral_generator_is_rejected_in_every_degree(self):
        borel = BorelMap("A3", omitted_nodes=1)
        bad = GeneratorAlgebra(borel, {"fraction": borel.weights[0] ** 2 / 2})
        for degree in (0, 1, 2, 3, 10):
            with self.assertRaises(ValueError):
                bad.integral_degree(degree)
            count("nonintegral_generator_rejections_even_at_gaps")
        mixed = GeneratorAlgebra(borel, {"good": 2 * borel.weights[0],
                                        "bad": borel.weights[0] / 2})
        with self.assertRaises(ValueError):
            mixed.integral_degree(0, include_ambient=True)
        count("nonintegral_generator_rejections_even_at_gaps")

    def test_weighted_gaps_above_top_and_zero_images(self):
        borel = BorelMap("A3", omitted_nodes=1)
        h = borel.weights[0]
        algebra = GeneratorAlgebra(borel, {"q": h ** 2})
        for degree, domain_dimension, ambient in ((1, 0, 1), (4, 1, 0), (5, 0, 0)):
            report = algebra.integral_degree(degree, include_ambient=True)
            self.assertEqual(report.domain_dimension, domain_dimension)
            self.assertEqual(report.rank, 0)
            self.assertEqual(report.smith_factors, ())
            self.assertEqual(report.torsion_factors, ())
            self.assertEqual(report.saturation_index, 1)
            self.assertTrue(report.is_saturated)
            self.assertEqual(report.ambient_dimension, ambient)
            self.assertEqual(report.free_cokernel_rank, ambient)
            self.assertEqual(report.full_ambient, ambient == 0)
            self.assert_lattice_report(report)
            self.assert_preimage(algebra, report, borel.expand(0))
            count("empty_weighted_or_above_top_components")
        report = algebra.integral_degree(1)
        self.assertFalse(report.contains(borel.expand(h)))
        line = BorelMap("A1")
        zero = GeneratorAlgebra(line, {"z": line.weights[0] ** 2})
        report = zero.integral_degree(2)
        self.assertEqual(report.matrix.dimensions(), (0, 1))
        self.assert_lattice_report(report)
        self.assert_preimage(zero, report, line.expand(0))

    def test_exceptional_integral_generators_against_explicit_g2(self):
        borel, oracle = BorelMap("G2", omitted_nodes=1), self.oracle("G2")
        algebra = GeneratorAlgebra(borel, {"h": borel.weights[0]})
        for degree in range(6):
            report = self.assert_independent_matrix(algebra, degree, "G2", (oracle.u,))
            self.assertEqual(report.smith_factors, (1 if degree < 3 else 2,))
            self.assert_preimage(algebra, report, borel.expand(borel.weights[0] ** degree))
        # Even a primitive divisor can have nonprimitive powers; this checks
        # actual Schubert coefficients rather than its Cartan denominator.
        self.assertTrue(any(algebra.integral_degree(d).saturation_index > 1
                            for d in range(2, 6)))
        count("G2_primitive_divisor_nonprimitive_power_check")

    def test_validation_membership_context_immutability_and_cache(self):
        borel = BorelMap("A2", omitted_nodes=1)
        algebra = GeneratorAlgebra(borel, {"h": 2 * borel.weights[0]})
        for degree in (-1, True, QQ(3) / 2, "2"):
            with self.assertRaises((TypeError, ValueError)):
                algebra.integral_degree(degree)
            count("invalid_degree_rejections")
        with self.assertRaises((TypeError, ValueError)):
            algebra.integral_degree(1, include_ambient=1)
        report = algebra.integral_degree(1)
        self.assertIs(report, algebra.integral_degree(1))
        self.assertIsNone(report.ambient_dimension)
        self.assertIsNone(report.free_cokernel_rank)
        self.assertIsNone(report.full_ambient)
        for operation in (lambda: report.contains(0), lambda: report.preimage(0)):
            with self.assertRaises(TypeError):
                operation()
            count("non_expansion_membership_rejections")
        for expansion in (borel.schubert([]),
                          borel.schubert([]) + borel.schubert([1]),
                          BorelMap("A2").schubert([1])):
            with self.assertRaises(ValueError):
                report.contains(expansion)
            with self.assertRaises(ValueError):
                report.preimage(expansion)
            count("wrong_degree_or_foreign_context_rejections", 2)
        with self.assertRaises(FrozenInstanceError):
            report.degree = 7
        for lattice in (report.matrix, report.image_lattice_basis,
                        report.saturation_lattice_basis):
            with self.assertRaises((ValueError, TypeError)):
                lattice[0, 0] = 7
            count("immutable_integer_matrices")
        algebra.clear_cache()
        repeated = algebra.integral_degree(1)
        self.assertIsNot(repeated, report)
        self.assertEqual(repeated.matrix, report.matrix)
        self.assertEqual(repeated.smith_factors, report.smith_factors)
        self.assert_preimage(algebra, repeated, borel.expand(2 * borel.weights[0]))
        count("cache_and_frozen_report_checks", 3)
        value = borel.schubert([1])
        for scalar in (ZZ(-2), ZZ(0), ZZ(1), ZZ(2)):
            self.assertEqual(scalar * value, value * scalar)
            self.assertEqual(scalar * value, int(scalar) * value)
            count("sage_integer_scalar_coercion_regressions")

    def test_default_does_not_enumerate_ambient_schubert_basis(self):
        borel = BorelMap("A4", omitted_nodes=2)
        algebra = GeneratorAlgebra(borel, {"u": borel.weights[1] ** 2})
        with patch.object(borel, "basis",
                          side_effect=AssertionError("ambient basis requested")):
            report = algebra.integral_degree(2)
            self.assertIsNone(report.ambient_dimension)
            self.assertIsNone(report.free_cokernel_rank)
            self.assertIsNone(report.full_ambient)
            self.assertEqual(report.rank, 1)
        with patch.object(borel, "basis", wraps=borel.basis) as ambient_basis:
            complete = algebra.integral_degree(2, include_ambient=True)
            self.assertTrue(ambient_basis.called)
        self.assertEqual(complete.ambient_dimension, 2)
        self.assertEqual(complete.free_cokernel_rank, 1)
        self.assertEqual(complete.matrix, report.matrix)
        self.assertIs(algebra.integral_degree(2), report)
        self.assertIs(algebra.integral_degree(2, include_ambient=True), complete)
        count("optional_ambient_basis_and_cache_checks", 2)


if __name__ == "__main__":
    unittest.main()
