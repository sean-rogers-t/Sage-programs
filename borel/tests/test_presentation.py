"""Complete kernel checks independent of the production staircase algorithm.

Expected ideals are obtained by lexicographic elimination from explicit
Cartesian type A/C/D coinvariant ideals and a type G2 weight orbit.  This checks the whole ideal, rather
than only comparing relations up to another chosen degree.  Exceptional
one-divisor examples use the independent ample-class/Hard Lefschetz result.

Run with ``PYTHONPATH=src sage -python -m unittest tests.test_presentation -v``.
Set ``BOREL_PRESENTATION_REPORT`` to write a machine-readable verification log.
"""

from dataclasses import FrozenInstanceError, replace
from itertools import combinations, product
import json
import os
from pathlib import Path
import time
import unittest

from sage.all import PolynomialRing, QQ, matrix, prod
from sage.env import SAGE_VERSION

from sage_borel import BorelMap, GeneratorAlgebra, TwoStepSymplectic


COUNTS = {}
CASES = []


def count(name, amount=1):
    COUNTS[name] = COUNTS.get(name, 0) + amount


def elementary(variables, degree):
    """Subset formula; no library symmetric-polynomial helper is used."""
    return sum((prod(subset) for subset in combinations(variables, degree)), 0)


def weighted_degree(monomial, weights):
    exponent, = monomial.exponents()
    return sum(int(a) * d for a, d in zip(exponent, weights))


class EliminationOracle:
    """A full elimination oracle using explicit independent Weyl invariants."""

    def __init__(self, algebra, family, size):
        self.algebra = algebra
        self.size = size
        self.ring = PolynomialRing(
            QQ, size, names=["x%d" % (i + 1) for i in range(size)], order="degrevlex")
        self.x = self.ring.gens()
        if family == "A":
            self.invariants = tuple(elementary(self.x, i) for i in range(1, size + 1))
            weight_images = tuple(sum(self.x[:i]) - QQ(i) / size * sum(self.x)
                                  for i in range(1, size))
        elif family == "C":
            self.invariants = tuple(elementary(tuple(x ** 2 for x in self.x), i)
                                    for i in range(1, size + 1))
            weight_images = tuple(sum(self.x[:i]) for i in range(1, size + 1))
        elif family == "D":
            self.invariants = tuple(elementary(tuple(x ** 2 for x in self.x), i)
                                    for i in range(1, size)) + (prod(self.x),)
            weight_images = tuple(sum(self.x[:i]) for i in range(1, size - 1)) + (
                (sum(self.x[:-1]) - self.x[-1]) / 2, sum(self.x) / 2)
        elif family == "G2":
            # Nonzero weights of the seven-dimensional representation, in
            # fundamental-weight coordinates.  These power sums generate the
            # rational invariant ring in degrees two and six.  Neither the
            # library's reflection matrices nor its divided differences enter
            # this independent construction.
            u, v = self.x
            orbit = (u, -u, u - v, v - u, 2 * u - v, v - 2 * u)
            self.invariants = tuple(sum(weight ** d for weight in orbit) for d in (2, 6))
            weight_images = self.x
        else:
            raise ValueError("oracle supports A, C, D and G2 only")
        self.generators = tuple(self.ring(f(*weight_images))
                                for f in algebra.generator_polynomials)
        self.coinvariant_basis = self.ring.ideal(self.invariants).groebner_basis()
        self.expected_ideal = self._eliminate()
        self.expected_basis = tuple(self.expected_ideal.groebner_basis())

    def _eliminate(self):
        # Under lex order x1 > ... > xn > formal variables, the x-free
        # members of a Groebner basis generate the entire contraction ideal.
        combined = PolynomialRing(
            QQ, self.size + len(self.generators),
            names=tuple("x%d" % (i + 1) for i in range(self.size)) +
                  self.algebra.generator_names, order="lex")
        variables = combined.gens()
        x, formal = variables[:self.size], variables[self.size:]
        lift = self.ring.hom(x, combined)
        ideal = combined.ideal([lift(f) for f in self.invariants] +
                               [y - lift(f) for y, f in zip(formal, self.generators)])
        basis = ideal.groebner_basis()
        contracted = []
        for relation in basis:
            if all(all(exponent[i] == 0 for i in range(self.size))
                   for exponent in relation.dict()):
                contracted.append(self.algebra.domain({tuple(exponent[self.size:]): coefficient
                                                       for exponent, coefficient in relation.dict().items()}))
        return self.algebra.domain.ideal(contracted)

    def reduce(self, polynomial):
        return self.ring(polynomial).reduce(self.coinvariant_basis)

    def image_matrix(self, monomials):
        images = []
        for monomial in monomials:
            exponent, = monomial.exponents()
            images.append(self.reduce(prod(f ** int(e)
                                           for f, e in zip(self.generators, exponent))))
        support = sorted({tuple(exponent) for image in images for exponent in image.dict()})
        coefficients = [{tuple(exponent): c for exponent, c in image.dict().items()}
                        for image in images]
        return matrix(QQ, len(images), len(support),
                      [row.get(exponent, 0) for row in coefficients for exponent in support])

    def standard_monomials(self):
        leading = [tuple(f.lm().exponents()[0]) for f in self.expected_basis]
        bounds = []
        for position in range(len(self.generators)):
            pure_powers = [exponent[position] for exponent in leading
                           if exponent[position] and all(not exponent[i]
                                                        for i in range(len(self.generators))
                                                        if i != position)]
            if not pure_powers:
                raise AssertionError("independent elimination ideal is not zero-dimensional")
            bounds.append(min(pure_powers))
        exponents = [exponent for exponent in product(*(range(bound) for bound in bounds))
                     if not any(all(a >= b for a, b in zip(exponent, power))
                                for power in leading)]
        return tuple(self.algebra.domain({exponent: 1}) for exponent in exponents)


class KernelPresentationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.started = time.time()
        cls.test_methods_run = 0

    def setUp(self):
        type(self).test_methods_run += 1
        type(self).test_result = self._outcome.result

    @classmethod
    def tearDownClass(cls):
        failed = any(isinstance(test, cls)
                     for test, _ in cls.test_result.failures + cls.test_result.errors)
        report = {"suite_result": "failed" if failed else "passed",
                  "test_methods_run": cls.test_methods_run,
                  "sage_version": SAGE_VERSION,
                  "seconds": round(time.time() - cls.started, 3),
                  "command": "PYTHONPATH=src sage -python -m unittest tests.test_presentation -v",
                  "independent_checks": COUNTS, "cases": CASES,
                  "expected_values_source": "Full lex elimination of Cartesian A/C/D coinvariant ideals and explicit G2 orbit invariants; ample-divisor Hard Lefschetz theorem"}
        print("\nPresentation verification: " + json.dumps(report, sort_keys=True), flush=True)
        if os.environ.get("BOREL_PRESENTATION_REPORT"):
            Path(os.environ["BOREL_PRESENTATION_REPORT"]).write_text(json.dumps(report, indent=2) + "\n")

    def assert_complete_kernel(self, algebra, family, size, label):
        oracle = EliminationOracle(algebra, family, size)
        expected_monomials = oracle.standard_monomials()
        expected_hilbert = [0] * (max(weighted_degree(m, algebra.generator_degrees)
                                     for m in expected_monomials) + 1)
        for monomial in expected_monomials:
            expected_hilbert[weighted_degree(monomial, algebra.generator_degrees)] += 1
        presentation = algebra.presentation()
        self.assertIsInstance(presentation.relations, tuple)
        self.assertIsInstance(presentation.standard_monomials, tuple)
        self.assertIsInstance(presentation.hilbert_coefficients, tuple)
        self.assertEqual(set(presentation.relations), set(oracle.expected_basis))
        self.assertEqual(set(presentation.standard_monomials), set(expected_monomials))
        self.assertEqual(presentation.dimension, len(expected_monomials))
        self.assertEqual(presentation.hilbert_coefficients, tuple(expected_hilbert))
        self.assertEqual(presentation.domain, algebra.domain)
        self.assertEqual(presentation.generator_degrees, algebra.generator_degrees)
        self.assertTrue(presentation.certificate.certified)
        series = presentation.hilbert_series
        t = series.parent().gen()
        self.assertEqual(series, sum(c * t ** d for d, c in enumerate(expected_hilbert)))
        self.assertEqual(series(t=1), presentation.dimension)
        for relation in presentation.relations:
            self.assertEqual(relation.parent(), algebra.domain)
            self.assertEqual(relation.lc(), 1)
            self.assertEqual(relation.reduce(oracle.expected_basis), 0)
            self.assertEqual(oracle.reduce(algebra.substitute(relation)(*
                self._weight_images(oracle, family, size))), 0)
            count("relations_independent_coinvariant_zero")
        self.assertEqual(oracle.image_matrix(presentation.standard_monomials).rank(),
                         presentation.dimension)
        # Include nonhomogeneous inputs and terms far beyond the top degree.
        variables = algebra.domain.gens()
        samples = [algebra.domain(0), algebra.domain(7),
                   sum((i + 2) * y ** (i + 1) for i, y in enumerate(variables)),
                   prod(1 + y + y ** 3 for y in variables),
                   sum(y ** 17 for y in variables)]
        for index, sample in enumerate(samples):
            self.assertEqual(presentation.reduce(sample), sample.reduce(oracle.expected_basis))
            if index < 3:
                self.assertEqual(algebra.evaluate(sample), algebra.evaluate(presentation.reduce(sample)))
            count("normal_forms_independent_groebner")
        self.assertEqual(set(presentation.kernel_ideal().groebner_basis()), set(oracle.expected_basis))
        count("entire_elimination_kernel_and_staircase")
        CASES.append({"case": label, "dimension": int(presentation.dimension),
                      "hilbert_coefficients": list(map(int, presentation.hilbert_coefficients)),
                      "groebner_basis_size": len(presentation.relations)})
        return presentation

    @staticmethod
    def _weight_images(oracle, family, size):
        x = oracle.x
        if family == "A":
            return tuple(sum(x[:i]) - QQ(i) / size * sum(x) for i in range(1, size))
        if family == "C":
            return tuple(sum(x[:i]) for i in range(1, size + 1))
        if family == "G2":
            return x
        return tuple(sum(x[:i]) for i in range(1, size - 1)) + (
            (sum(x[:-1]) - x[-1]) / 2, sum(x) / 2)

    def test_type_a_full_flag_and_grassmannian(self):
        borel = BorelMap("A2")
        algebra = GeneratorAlgebra(borel, dict(zip(("u", "v"), borel.weights)))
        presentation = self.assert_complete_kernel(algebra, "A", 3, "A2 full flag")
        self.assertEqual(presentation.dimension, 6)
        self.assertEqual(presentation.hilbert_coefficients, (1, 2, 2, 1))
        borel = BorelMap("A3", omitted_nodes=2)
        w1, w2, _ = borel.weights
        algebra = GeneratorAlgebra(borel, {"c1": w2, "c2": w1 * (w2 - w1)})
        presentation = self.assert_complete_kernel(algebra, "A", 4, "Gr(2,4)")
        self.assertEqual(presentation.dimension, 6)
        self.assertEqual(presentation.hilbert_coefficients, (1, 1, 2, 1, 1))

    def test_type_c_full_flag_tail_and_joint_images(self):
        borel = BorelMap("C2")
        algebra = GeneratorAlgebra(borel, dict(zip(("u", "v"), borel.weights)))
        presentation = self.assert_complete_kernel(algebra, "C", 2, "C2 full flag")
        self.assertEqual(presentation.dimension, 8)
        family = TwoStepSymplectic(3, 1, 1)
        presentation = self.assert_complete_kernel(family.joint_algebra(), "C", 3,
                                                 "C3 joint middle and tail")
        self.assertEqual(presentation.dimension, 12)
        self.assertEqual(presentation.hilbert_coefficients, (1, 1, 2, 2, 2, 2, 1, 1))
        presentation = self.assert_complete_kernel(family.tail_algebra(), "C", 3, "C3 tail")
        self.assertEqual(presentation.dimension, 3)
        self.assertEqual(presentation.hilbert_coefficients, (1, 0, 1, 0, 1))

    def test_type_d_tail_retains_euler_generator(self):
        borel = BorelMap("D4", omitted_nodes=2)
        _, w2, w3, w4 = borel.weights
        x3, x4 = w3 + w4 - w2, w4 - w3
        algebra = GeneratorAlgebra(borel, {"q": x3 ** 2 + x4 ** 2, "euler": x3 * x4})
        presentation = self.assert_complete_kernel(algebra, "D", 4, "D4 tail with Euler")
        q, euler = algebra.domain.gens()
        self.assertNotEqual(presentation.reduce(euler), 0)
        self.assertNotEqual(presentation.reduce(euler ** 2), 0)
        # The Euler class is independent of q in degree two; its sign cannot
        # be reconstructed from elementary polynomials in squared coordinates.
        self.assertEqual(presentation.hilbert_coefficients[2], 2)
        self.assertNotEqual(presentation.reduce(q - euler), 0)
        count("type_d_euler_independence", 4)

    def test_dependent_generators_and_weighted_holes(self):
        borel = BorelMap("A3", omitted_nodes=1)
        w = borel.weights[0]
        algebra = GeneratorAlgebra(borel, {"u": w, "v": w ** 2})
        presentation = self.assert_complete_kernel(algebra, "A", 4, "P3 dependent generators")
        self.assertEqual(presentation.hilbert_coefficients, (1, 1, 1, 1))
        borel = BorelMap("A4", omitted_nodes=1)
        w = borel.weights[0]
        algebra = GeneratorAlgebra(borel, {"u": w ** 2, "v": w ** 3})
        presentation = self.assert_complete_kernel(algebra, "A", 5, "P4 weighted subalgebra")
        self.assertEqual(presentation.hilbert_coefficients, (1, 0, 1, 1, 1))
        self.assertEqual(presentation.dimension, 4)
        algebra = GeneratorAlgebra(borel, {"u": w ** 2, "v": w})
        presentation = self.assert_complete_kernel(algebra, "A", 5,
                                                 "P4 descending generator weights")
        u, v = algebra.domain.gens()
        self.assertIn(v ** 2 - u, presentation.relations)
        self.assertEqual(presentation.hilbert_coefficients, (1, 1, 1, 1, 1))
        algebra = GeneratorAlgebra(borel, {"u": w, "v": w ** 2, "z": w ** 3})
        presentation = self.assert_complete_kernel(algebra, "A", 5,
                                                 "P4 dependent weights one two three")
        self.assertEqual(presentation.hilbert_coefficients, (1, 1, 1, 1, 1))

    def test_zero_images_and_generators_above_ambient_top_degree(self):
        for exponent in (2, 37):
            borel = BorelMap("A1")
            algebra = GeneratorAlgebra(borel, {"z": borel.weights[0] ** exponent})
            presentation = self.assert_complete_kernel(algebra, "A", 2,
                                                     "A1 zero image degree %d" % exponent)
            z, = algebra.domain.gens()
            self.assertEqual(presentation.relations, (z,))
            self.assertEqual(presentation.standard_monomials, (algebra.domain.one(),))
            self.assertEqual(presentation.dimension, 1)
            self.assertEqual(presentation.hilbert_coefficients, (1,))

    def test_exceptional_ample_divisors(self):
        # G2/P1 and G2/P2 have complex dimension five; F4/P4 has
        # dimension fifteen.  An ample generator has independent nonzero
        # powers through the top degree by Hard Lefschetz, then vanishes.
        for cartan, node, dimension in (("G2", 1, 5), ("G2", 2, 5), ("F4", 4, 15)):
            borel = BorelMap(cartan, omitted_nodes=node)
            algebra = GeneratorAlgebra(borel, {"h": borel.weights[node - 1]})
            h, = algebra.domain.gens()
            presentation = algebra.presentation()
            self.assertEqual(presentation.relations, (h ** (dimension + 1),))
            self.assertEqual(set(presentation.standard_monomials),
                             {h ** i for i in range(dimension + 1)})
            self.assertEqual(presentation.dimension, dimension + 1)
            self.assertEqual(presentation.hilbert_coefficients, (1,) * (dimension + 1))
            self.assertEqual(presentation.reduce(1 + h + h ** (dimension + 7)), 1 + h)
            count("exceptional_ample_class_full_kernel")
            CASES.append({"case": "%s/P%d ample divisor" % (cartan, node),
                          "dimension": dimension + 1,
                          "hilbert_coefficients": [1] * (dimension + 1),
                          "groebner_basis_size": 1})

    def test_g2_full_flag_from_explicit_weight_orbit_invariants(self):
        borel = BorelMap("G2")
        algebra = GeneratorAlgebra(borel, dict(zip(("u", "v"), borel.weights)))
        presentation = self.assert_complete_kernel(algebra, "G2", 2,
                                                 "G2 full flag orbit invariants")
        u, v = algebra.domain.gens()
        self.assertEqual(set(presentation.relations),
                         {v ** 6, u ** 2 - u * v + QQ(1) / 3 * v ** 2})
        self.assertEqual(presentation.dimension, 12)
        self.assertEqual(presentation.hilbert_coefficients, (1, 2, 2, 2, 2, 2, 1))
        count("g2_full_flag_explicit_orbit_invariants")

    def test_immutable_results_fresh_ideals_and_cache_reset(self):
        borel = BorelMap("A2", omitted_nodes=1)
        algebra = GeneratorAlgebra(borel, {"h": borel.weights[0]})
        presentation = algebra.presentation()
        self.assertIs(algebra.presentation(), presentation)
        with self.assertRaises((FrozenInstanceError, AttributeError)):
            presentation.dimension = 500
        first = presentation.kernel_ideal()
        second = presentation.kernel_ideal()
        self.assertIsNot(first, second)
        self.assertEqual(first, second)
        relations = presentation.relations
        monomials = presentation.standard_monomials
        algebra.clear_cache()
        repeated = algebra.presentation()
        self.assertIsNot(repeated, presentation)
        self.assertEqual(repeated.relations, relations)
        self.assertEqual(repeated.standard_monomials, monomials)
        self.assertEqual(repeated.hilbert_coefficients, presentation.hilbert_coefficients)
        certificate = presentation.verify()
        self.assertTrue(certificate.certified)
        for field in ("relation_images_zero", "basis_images_independent",
                      "standard_monomials_reduced", "border_closed", "contains_unit",
                      "groebner_basis_reduced", "summary_data_correct"):
            self.assertTrue(getattr(certificate, field))
        with self.assertRaises((FrozenInstanceError, AttributeError)):
            certificate.contains_unit = False
        # Removing nonconstant basis vectors preserves independence but loses
        # closure under h.  A full-kernel certificate must reject that claim.
        incomplete = replace(presentation, standard_monomials=(algebra.domain.one(),))
        self.assertFalse(incomplete.verify().certified)
        h, = algebra.domain.gens()
        incorrect = replace(presentation, relations=(h,))
        self.assertFalse(incorrect.verify().relation_images_zero)
        dependent = replace(presentation, standard_monomials=(algebra.domain.one(), h, h))
        self.assertFalse(dependent.verify().basis_images_independent)
        no_unit = replace(presentation, standard_monomials=(h, h ** 2))
        self.assertFalse(no_unit.verify().contains_unit)
        # A forged low geometric bound must not hide a nonzero divisor image.
        forged = replace(presentation, relations=(h,),
                         standard_monomials=(algebra.domain.one(),),
                         hilbert_coefficients=(1,),
                         hilbert_series=presentation.hilbert_series.parent().one(),
                         certificate=replace(presentation.certificate, ambient_top_degree=0))
        rechecked = forged.verify()
        self.assertEqual(rechecked.ambient_top_degree, 2)
        self.assertFalse(rechecked.relation_images_zero)
        self.assertFalse(rechecked.certified)
        for altered in (
                replace(presentation, hilbert_coefficients=(1,)),
                replace(presentation, hilbert_series=presentation.hilbert_series + 1),
                replace(presentation, generator_degrees=(2,)),
                replace(presentation, domain=PolynomialRing(QQ, 1, names=("other",)))):
            self.assertFalse(altered.verify().summary_data_correct)
            self.assertFalse(altered.verify().certified)
        serialized = presentation.as_dict()
        self.assertIn("schema_version", serialized)
        json.dumps(serialized)
        count("immutable_fresh_ideal_cache_regressions", 7)
        count("completeness_certificate_and_invalid_claim_rejection", 10)
        count("json_certificate_serialization")


if __name__ == "__main__":
    unittest.main()
