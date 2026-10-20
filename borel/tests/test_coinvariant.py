"""Complete-kernel regressions for the direct coinvariant backend.

Small cases compare both production algorithms and independently eliminated
Cartesian invariant ideals.  The rank-nine regression prevents any call to
Schubert expansion; its full kernel follows from the quaternionic projective
space presentation Q[q]/(q^9), with q in polynomial degree two.

Run with ``PYTHONPATH=src sage -python -m unittest tests.test_coinvariant -v``.
Set ``BOREL_COINVARIANT_REPORT`` to write a JSON verification log.
"""

from dataclasses import replace
import json
import os
from pathlib import Path
import time
import unittest
from unittest.mock import patch

from sage.all import CartanType
from sage.env import SAGE_VERSION

from sage_borel import BorelMap, GeneratorAlgebra, TwoStepSymplectic
from tests.test_presentation import EliminationOracle, weighted_degree


CASES = []


class CoinvariantPresentationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.started = time.time()
        cls.methods_run = 0

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
            "command": "PYTHONPATH=src sage -python -m unittest tests.test_coinvariant -v",
            "cases": CASES,
            "expected_values_source": (
                "Independent Cartesian elimination in A/C/D; explicit signed "
                "B2 invariants; BGG crosschecks; Q[q]/(q^9) rank-nine tail"),
        }
        print("\nCoinvariant verification: " + json.dumps(report, sort_keys=True), flush=True)
        if os.environ.get("BOREL_COINVARIANT_REPORT"):
            Path(os.environ["BOREL_COINVARIANT_REPORT"]).write_text(
                json.dumps(report, indent=2) + "\n")

    def assert_backends_agree(self, algebra, label, oracle_family=None, oracle_size=None):
        started = time.time()
        direct = algebra.presentation(backend="coinvariant")
        bgg = algebra.presentation(backend="bgg")
        self.assertEqual(direct.backend, "coinvariant")
        self.assertEqual(bgg.backend, "bgg")
        self.assertEqual(direct.certificate.backend, "coinvariant")
        self.assertEqual(bgg.certificate.backend, "bgg")
        self.assertEqual(direct.relations, bgg.relations)
        self.assertEqual(direct.standard_monomials, bgg.standard_monomials)
        self.assertEqual(direct.hilbert_coefficients, bgg.hilbert_coefficients)
        self.assertEqual(direct.hilbert_series, bgg.hilbert_series)
        self.assertEqual(direct.dimension, bgg.dimension)
        for presentation in (direct, bgg):
            self.assertTrue(presentation.certificate.certified)
            self.assertTrue(presentation.verify(backend="coinvariant").certified)
            self.assertTrue(presentation.verify(backend="bgg").certified)
        if oracle_family is not None:
            oracle = EliminationOracle(algebra, oracle_family, oracle_size)
            monomials = oracle.standard_monomials()
            coefficients = [0] * (max(weighted_degree(m, algebra.generator_degrees)
                                      for m in monomials) + 1)
            for monomial in monomials:
                coefficients[weighted_degree(monomial, algebra.generator_degrees)] += 1
            self.assertEqual(set(direct.relations), set(oracle.expected_basis))
            self.assertEqual(set(direct.standard_monomials), set(monomials))
            self.assertEqual(direct.hilbert_coefficients, tuple(coefficients))
            self.assertEqual(oracle.image_matrix(direct.standard_monomials).rank(),
                             direct.dimension)
        CASES.append({"case": label, "dimension": int(direct.dimension),
                      "relations": len(direct.relations),
                      "independent_elimination": oracle_family is not None,
                      "seconds": round(time.time() - started, 3)})
        return direct

    def test_type_a_full_flag_and_grassmannian(self):
        borel = BorelMap("A2")
        algebra = GeneratorAlgebra(borel, dict(zip(("u", "v"), borel.weights)))
        full = self.assert_backends_agree(algebra, "A2 full flag", "A", 3)
        self.assertEqual(full.hilbert_coefficients, (1, 2, 2, 1))
        borel = BorelMap("A3", omitted_nodes=2)
        w1, w2, _ = borel.weights
        algebra = GeneratorAlgebra(borel, {"c1": w2, "c2": w1 * (w2 - w1)})
        grassmannian = self.assert_backends_agree(algebra, "Gr(2,4)", "A", 4)
        self.assertEqual(grassmannian.hilbert_coefficients, (1, 1, 2, 1, 1))

    def test_type_b_half_weight_convention(self):
        borel = BorelMap("B2")
        algebra = GeneratorAlgebra(borel, dict(zip(("u", "v"), borel.weights)))
        presentation = self.assert_backends_agree(algebra, "B2 full flag")
        u, v = algebra.domain.gens()
        # Fundamental weights are omega1=x1 and omega2=(x1+x2)/2.
        # These signed invariants are specified independently of the production
        # coordinate conversion and normal-form implementation.
        x1, x2 = u, 2 * v - u
        expected = algebra.domain.ideal([x1 ** 2 + x2 ** 2, x1 ** 2 * x2 ** 2])
        self.assertEqual(set(presentation.relations), set(expected.groebner_basis()))
        self.assertEqual(presentation.hilbert_coefficients, (1, 2, 2, 2, 1))
        borel = BorelMap("B2", omitted_nodes=1)
        algebra = GeneratorAlgebra(borel, {"h": borel.weights[0]})
        presentation = self.assert_backends_agree(algebra, "B2/P1 ample class")
        h, = algebra.domain.gens()
        self.assertEqual(presentation.relations, (h ** 4,))
        self.assertEqual(presentation.hilbert_coefficients, (1, 1, 1, 1))

    def test_type_c_full_flag_and_joint_image(self):
        borel = BorelMap("C2")
        algebra = GeneratorAlgebra(borel, dict(zip(("u", "v"), borel.weights)))
        full = self.assert_backends_agree(algebra, "C2 full flag", "C", 2)
        self.assertEqual(full.dimension, 8)
        algebra = TwoStepSymplectic(3, 1, 1).joint_algebra()
        joint = self.assert_backends_agree(algebra, "C3 middle and tail", "C", 3)
        self.assertEqual(joint.dimension, 12)
        self.assertEqual(joint.hilbert_coefficients, (1, 1, 2, 2, 2, 2, 1, 1))

    def test_type_d_euler_generator_is_retained(self):
        borel = BorelMap("D4", omitted_nodes=2)
        _, w2, w3, w4 = borel.weights
        x3, x4 = w3 + w4 - w2, w4 - w3
        algebra = GeneratorAlgebra(borel, {"q": x3 ** 2 + x4 ** 2, "euler": x3 * x4})
        presentation = self.assert_backends_agree(algebra, "D4 tail with Euler", "D", 4)
        q, euler = algebra.domain.gens()
        self.assertNotEqual(presentation.reduce(euler), 0)
        self.assertNotEqual(presentation.reduce(euler ** 2), 0)
        self.assertNotEqual(presentation.reduce(q - euler), 0)
        self.assertEqual(presentation.hilbert_coefficients[2], 2)

    def test_weighted_dependent_and_zero_images(self):
        borel = BorelMap("A4", omitted_nodes=1)
        w = borel.weights[0]
        algebra = GeneratorAlgebra(borel, {"u": w ** 2, "v": w ** 3})
        presentation = self.assert_backends_agree(algebra, "P4 weighted subalgebra", "A", 5)
        self.assertEqual(presentation.hilbert_coefficients, (1, 0, 1, 1, 1))
        algebra = GeneratorAlgebra(borel, {"u": w ** 2, "v": w})
        presentation = self.assert_backends_agree(algebra, "P4 descending weights", "A", 5)
        u, v = algebra.domain.gens()
        self.assertIn(v ** 2 - u, presentation.relations)
        for exponent in (2, 37):
            borel = BorelMap("A1")
            algebra = GeneratorAlgebra(borel, {"z": borel.weights[0] ** exponent})
            presentation = self.assert_backends_agree(algebra, "A1 zero degree %d" % exponent,
                                                       "A", 2)
            z, = algebra.domain.gens()
            self.assertEqual(presentation.relations, (z,))
            self.assertEqual(presentation.standard_monomials, (algebra.domain.one(),))

    def test_rank_nine_tail_full_kernel_without_bgg(self):
        started = time.time()
        borel = BorelMap("C9", omitted_nodes=8)
        algebra = GeneratorAlgebra(borel, {"q": borel.tail_polynomial(1, 1)})
        q, = algebra.domain.gens()
        with patch.object(BorelMap, "expand", side_effect=AssertionError("BGG expansion called")), \
                patch.object(GeneratorAlgebra, "_evaluate_monomial",
                             side_effect=AssertionError("BGG monomial image called")):
            presentation = algebra.presentation(backend="coinvariant")
            self.assertIs(algebra.presentation(), presentation)
            self.assertTrue(presentation.verify().certified)
            self.assertTrue(presentation.verify(backend="coinvariant").certified)
        self.assertEqual(presentation.relations, (q ** 9,))
        self.assertEqual(presentation.standard_monomials, tuple(q ** i for i in range(9)))
        self.assertEqual(presentation.dimension, 9)
        self.assertEqual(presentation.hilbert_coefficients, (1, 0) * 8 + (1,))
        self.assertEqual(presentation.reduce(1 + q + q ** 9 + q ** 25), 1 + q)
        self.assertEqual(presentation.certificate.backend, "coinvariant")
        CASES.append({"case": "C9/P8 one-coordinate tail", "dimension": 9,
                      "relations": ["q^9"], "BGG_calls_permitted": False,
                      "construction_and_verify_seconds": round(time.time() - started, 3)})

    def test_backend_cache_and_serialization(self):
        borel = BorelMap("A2", omitted_nodes=1)
        algebra = GeneratorAlgebra(borel, {"h": borel.weights[0]})
        direct = algebra.presentation()
        bgg = algebra.presentation(backend="bgg")
        self.assertIs(algebra.presentation(backend="auto"), direct)
        self.assertIs(algebra.presentation(backend="coinvariant"), direct)
        self.assertIs(algebra.presentation(backend="bgg"), bgg)
        self.assertIsNot(direct, bgg)
        for presentation in (direct, bgg):
            serialized = presentation.as_dict()
            self.assertEqual(serialized["backend"], presentation.backend)
            self.assertEqual(serialized["certificate"]["backend"], presentation.backend)
            self.assertEqual(presentation.verify().backend, presentation.backend)
            json.dumps(serialized)
        self.assertEqual(direct.verify(backend="bgg").backend, "bgg")
        self.assertEqual(bgg.verify(backend="coinvariant").backend, "coinvariant")
        algebra.clear_cache()
        repeated = algebra.presentation()
        self.assertIsNot(repeated, direct)
        self.assertIsNot(algebra.presentation(backend="bgg"), bgg)
        self.assertEqual(repeated.relations, direct.relations)

    def test_invalid_backend_and_exceptional_fallback(self):
        borel = BorelMap("G2", omitted_nodes=1)
        algebra = GeneratorAlgebra(borel, {"h": borel.weights[0]})
        presentation = algebra.presentation()
        self.assertEqual(presentation.backend, "bgg")
        self.assertIs(algebra.presentation(backend="bgg"), presentation)
        self.assertTrue(presentation.verify().certified)
        h, = algebra.domain.gens()
        self.assertEqual(presentation.relations, (h ** 6,))
        with self.assertRaises(NotImplementedError):
            algebra.presentation(backend="coinvariant")
        for backend in ("invalid", "BGG", "", None, True):
            with self.assertRaises((ValueError, TypeError)):
                algebra.presentation(backend=backend)
        for backend in ("invalid", "BGG", "", True):
            with self.assertRaises((ValueError, TypeError)):
                presentation.verify(backend=backend)

    def test_relabeled_classical_falls_back_without_wrong_coordinates(self):
        kind = CartanType(["B", 2]).relabel({1: 2, 2: 1})
        borel = BorelMap(kind)
        algebra = GeneratorAlgebra(borel, dict(zip(("u", "v"), borel.weights)))
        presentation = algebra.presentation()
        self.assertEqual(presentation.backend, "bgg")
        self.assertIs(algebra.presentation(backend="bgg"), presentation)
        self.assertTrue(presentation.verify().certified)
        self.assertEqual(presentation.dimension, 8)
        self.assertEqual(presentation.hilbert_coefficients, (1, 2, 2, 2, 1))
        with self.assertRaises(NotImplementedError):
            algebra.presentation(backend="coinvariant")
        with self.assertRaises(NotImplementedError):
            presentation.verify(backend="coinvariant")

    def test_tampered_coinvariant_certificate_is_rejected(self):
        borel = BorelMap("A2", omitted_nodes=1)
        algebra = GeneratorAlgebra(borel, {"h": borel.weights[0]})
        presentation = algebra.presentation(backend="coinvariant")
        h, = algebra.domain.gens()
        incomplete = replace(presentation, standard_monomials=(algebra.domain.one(),))
        self.assertFalse(incomplete.verify().certified)
        incorrect = replace(presentation, relations=(h,))
        self.assertFalse(incorrect.verify().relation_images_zero)
        dependent = replace(presentation, standard_monomials=(algebra.domain.one(), h, h))
        self.assertFalse(dependent.verify().basis_images_independent)
        forged = replace(presentation, relations=(h,),
                         standard_monomials=(algebra.domain.one(),),
                         hilbert_coefficients=(1,),
                         hilbert_series=presentation.hilbert_series.parent().one(),
                         certificate=replace(presentation.certificate, ambient_top_degree=0))
        rechecked = forged.verify()
        self.assertEqual(rechecked.ambient_top_degree, 2)
        self.assertFalse(rechecked.relation_images_zero)
        self.assertFalse(rechecked.certified)
        altered = replace(presentation, generator_degrees=(2,))
        self.assertFalse(altered.verify().summary_data_correct)
        self.assertFalse(altered.verify().certified)


if __name__ == "__main__":
    unittest.main()
