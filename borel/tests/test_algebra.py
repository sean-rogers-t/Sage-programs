"""Independent finite-rank checks for generator images and relations.

Run from ``borel/`` with ``PYTHONPATH=src sage -python -m unittest
tests.test_algebra -v``.  This suite obtains its expected ranks and kernels
from explicit coinvariant ideals in orthonormal x coordinates, rather than
from Schubert multiplication or the library's divided differences.

Set ``BOREL_ALGEBRA_REPORT`` to write the verification counts as JSON.
"""

from dataclasses import FrozenInstanceError
from itertools import product
import json
import os
from pathlib import Path
import platform
import time
import unittest

from sage.all import PolynomialRing, PowerSeriesRing, QQ, matrix, prod
from sage.env import SAGE_VERSION

from sage_borel import BorelMap, GeneratorAlgebra, TwoStepSymplectic, two_step_symplectic


COUNTS = {}
CASE_RESULTS = []


def count(name, amount=1):
    COUNTS[name] = COUNTS.get(name, 0) + amount


def elementary_explicit(variables, degree):
    """Subset definition, separate from the library's coefficient recursion."""
    from itertools import combinations
    return sum((prod(subset) for subset in combinations(variables, degree)), 0)


def weighted_exponents(weights, degree):
    """Brute-force Cartesian product, independent of the production enumerator."""
    return tuple(exponents for exponents in product(
        *(range(degree // weight + 1) for weight in weights))
        if sum(e * weight for e, weight in zip(exponents, weights)) == degree)


def monomial(ring, exponents):
    return ring(prod(x ** exponent for x, exponent in zip(ring.gens(), exponents)))


def coefficient_dict(polynomial):
    # Sage's ETuple monomial keys have a different hash from Python tuples.
    return {tuple(alpha): coefficient for alpha, coefficient in polynomial.dict().items()}


class CoinvariantOracle:
    """Characteristic-zero coinvariant quotient in ordinary x coordinates.

    Type A uses all n elementary symmetric polynomials in x1,...,xn.
    Type C uses all n elementary symmetric polynomials in x1^2,...,xn^2.
    Sage/Singular supplies only exact Groebner reduction and linear algebra.
    No Weyl-group action, root data, BGG operator or library expansion occurs.
    """

    def __init__(self, family, size):
        self.family, self.size = family, size
        self.ring = PolynomialRing(QQ, size, names=["x%d" % i for i in range(1, size + 1)],
                                   order="degrevlex")
        self.x = self.ring.gens()
        invariant_variables = self.x if family == "A" else tuple(x ** 2 for x in self.x)
        self.ideal = self.ring.ideal([
            elementary_explicit(invariant_variables, degree) for degree in range(1, size + 1)])
        # Cache one independently calculated Groebner basis for every case.
        self.groebner = self.ideal.groebner_basis()

    def reduce(self, polynomial):
        return self.ring(polynomial).reduce(self.groebner)

    def omega(self, node):
        if self.family == "A":
            return sum(self.x[:node]) - QQ(node) / self.size * sum(self.x)
        return sum(self.x[:node])

    def input_matrix(self, generators, exponents):
        """Rows are input monomials; this is opposite to the library convention."""
        reductions = [self.reduce(prod(f ** e for f, e in zip(generators, alpha)))
                      for alpha in exponents]
        support = sorted({alpha for f in reductions for alpha in f.dict()})
        return matrix(QQ, len(reductions), len(support),
                      [f.dict().get(alpha, 0) for f in reductions for alpha in support])

    def swap(self, polynomial, left, right):
        images = list(self.x)
        images[left], images[right] = images[right], images[left]
        return self.ring(polynomial)(*images)

    def flip(self, polynomial, position):
        images = list(self.x)
        images[position] = -images[position]
        return self.ring(polynomial)(*images)


class GeneratorAlgebraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.started = time.time()
        cls.oracles = {}
        cls.test_methods_run = 0

    def setUp(self):
        type(self).test_methods_run += 1
        type(self).test_result = self._outcome.result

    @classmethod
    def tearDownClass(cls):
        report = {
            "suite_result": "passed" if not any(isinstance(test, cls)
                                                  for test, _ in cls.test_result.failures + cls.test_result.errors)
                            else "failed",
            "test_methods_run": cls.test_methods_run,
            "sage_version": SAGE_VERSION,
            "python_version": platform.python_version(),
            "command": "PYTHONPATH=src sage -python -m unittest tests.test_algebra -v",
            "seconds": round(time.time() - cls.started, 3),
            "independent_checks": COUNTS,
            "finite_rank_cases": CASE_RESULTS,
            "expected_values_source": "Explicit type A/C invariant ideals; exact independent Groebner reduction",
        }
        summary = {key: value for key, value in report.items() if key != "finite_rank_cases"}
        print("\nAlgebra verification: " + json.dumps(summary, sort_keys=True), flush=True)
        if os.environ.get("BOREL_ALGEBRA_REPORT"):
            Path(os.environ["BOREL_ALGEBRA_REPORT"]).write_text(json.dumps(report, indent=2) + "\n")

    def oracle(self, family, size):
        key = family, size
        if key not in self.oracles:
            self.oracles[key] = CoinvariantOracle(family, size)
        return self.oracles[key]

    def assert_independent_degree(self, algebra, oracle, input_generators, degree,
                                  label, include_ambient=False):
        weights = tuple(int(f.degree()) for f in input_generators)
        exponents = weighted_exponents(weights, degree)
        expected_monomials = tuple(monomial(algebra.domain, alpha) for alpha in exponents)
        report = algebra.degree(degree, include_ambient=include_ambient)
        self.assertEqual(tuple(algebra.monomials(degree)), expected_monomials)
        self.assertEqual(report.monomials, expected_monomials)
        quotient_matrix = oracle.input_matrix(input_generators, exponents)
        self.assertEqual(report.domain_dimension, len(exponents))
        self.assertEqual(report.image_dimension, quotient_matrix.rank())
        self.assertEqual(report.rank, quotient_matrix.rank())
        self.assertEqual(report.kernel_dimension, len(exponents) - quotient_matrix.rank())
        self.assertEqual(report.matrix.ncols(), len(exponents))
        self.assertEqual(report.matrix.nrows(), len(report.schubert_indices))
        self.assertEqual(report.matrix.rank(), report.rank)
        relation_vectors = [
            [coefficient_dict(relation).get(alpha, QQ.zero()) for alpha in exponents]
            for relation in report.relations]
        relation_matrix = matrix(QQ, len(relation_vectors), len(exponents),
                                 [c for row in relation_vectors for c in row])
        self.assertEqual(relation_matrix.row_space(), quotient_matrix.left_kernel())
        # Compare the full kernel space, not just individually vanishing relations.
        self.assertEqual(relation_matrix.row_space(), report.matrix.right_kernel())
        for relation in report.relations:
            substituted = algebra.substitute(relation)
            images = [oracle.omega(i) for i in range(1, algebra.borel_map.rank + 1)]
            self.assertEqual(oracle.reduce(substituted(*images)), 0)
            self.assertEqual(dict(algebra.evaluate(relation)), {})
            count("relations_independently_zero")
        CASE_RESULTS.append({"case": label, "degree": degree,
                             "domain_dimension": len(exponents), "rank": int(report.rank),
                             "kernel_dimension": int(report.kernel_dimension)})
        count("groebner_rank_and_whole_kernel")
        return report

    def test_weighted_monomial_counts_against_generating_functions(self):
        b = BorelMap("A3")
        w = b.weights
        generator_sets = ({"u": w[0]}, {"u": w[0], "v": w[1]},
                          {"u": w[0], "v": w[1] ** 2, "z": w[2] ** 3},
                          {"u": w[0] ** 2, "v": w[1] ** 4})
        for generators in generator_sets:
            algebra = GeneratorAlgebra(b, generators)
            weights = algebra.generator_degrees
            series_ring = PowerSeriesRing(QQ, "t", default_prec=13)
            t = series_ring.gen()
            series = prod(1 / (1 - t ** d) for d in weights)
            for degree in range(13):
                with self.subTest(weights=weights, degree=degree):
                    monomials = tuple(algebra.monomials(degree))
                    self.assertEqual(len(monomials), series[degree])
                    self.assertEqual(monomials, tuple(monomial(algebra.domain, alpha)
                                                     for alpha in weighted_exponents(weights, degree)))
                    self.assertEqual(len(set(monomials)), len(monomials))
                    count("weighted_monomial_generating_function")

    def test_projective_spaces_have_exact_truncation(self):
        for dimension in (1, 2, 3, 4):
            b = BorelMap("A%d" % dimension, omitted_nodes=1)
            algebra = GeneratorAlgebra(b, {"h": b.weights[0]})
            oracle = self.oracle("A", dimension + 1)
            reports = []
            for degree in range(dimension + 3):
                report = self.assert_independent_degree(
                    algebra, oracle, (oracle.omega(1),), degree,
                    "P%d" % dimension, include_ambient=True)
                self.assertEqual(report.rank, int(degree <= dimension))
                self.assertEqual(report.ambient_dimension, int(degree <= dimension))
                reports.append(report)
                count("projective_space_truncation")
            self.assertEqual(algebra.first_relation_degree(dimension), None)
            self.assertEqual(algebra.first_relation_degree(dimension + 1), dimension + 1)
            self.assertEqual(algebra.hilbert_function(dimension + 2),
                             tuple([1] * (dimension + 1) + [0, 0]))
            self.assertEqual(tuple(algebra.through_degree(dimension + 2, include_ambient=True)),
                             tuple(reports))

    def test_full_a2_coinvariant_kernel_and_nonsquare_orientation(self):
        b = BorelMap("A2")
        algebra = GeneratorAlgebra(b, dict(zip(("u", "v"), b.weights)))
        oracle = self.oracle("A", 3)
        for degree in range(8):
            self.assert_independent_degree(algebra, oracle,
                                           (oracle.omega(1), oracle.omega(2)),
                                           degree, "A2 full flag", include_ambient=True)
        report = algebra.degree(2, include_ambient=True)
        self.assertEqual(report.matrix.dimensions(), (2, 3))
        u, v = algebra.domain.gens()
        self.assertEqual(algebra.evaluate(u ** 2 - u * v + v ** 2).as_dict(), {})
        self.assertEqual(algebra.evaluate(u ** 2 * v - u * v ** 2).as_dict(), {})
        self.assertEqual(algebra.hilbert_function(7), (1, 2, 2, 1, 0, 0, 0, 0))
        count("a2_explicit_invariant_relations", 2)

    def test_redundant_weighted_generators_and_inhomogeneous_evaluation(self):
        b = BorelMap("A3", omitted_nodes=1)
        w = b.weights[0]
        algebra = GeneratorAlgebra(b, {"u": w, "v": w ** 2})
        oracle = self.oracle("A", 4)
        u, v = algebra.domain.gens()
        self.assertEqual(algebra.generator_names, ("u", "v"))
        self.assertEqual(algebra.generator_degrees, (1, 2))
        self.assertEqual(algebra.generator_polynomials, (w, w ** 2))
        self.assertEqual(algebra.substitute(7 + 2 * u - 3 * v + u * v),
                         7 + 2 * w - 3 * w ** 2 + w ** 3)
        self.assertEqual(dict(algebra.evaluate(7 + 2 * u - 3 * v + u * v)),
                         dict(b.expand(7 + 2 * w - 3 * w ** 2 + w ** 3)))
        self.assertEqual(dict(algebra.evaluate(v - u ** 2)), {})
        for degree in range(8):
            self.assert_independent_degree(algebra, oracle,
                                           (oracle.omega(1), oracle.omega(1) ** 2),
                                           degree, "redundant P3 generators")
        self.assertEqual(algebra.first_relation_degree(1), None)
        self.assertEqual(algebra.first_relation_degree(2), 2)
        count("inhomogeneous_and_weighted_regressions", 7)

    def test_zero_image_generator_and_empty_degree_component(self):
        b = BorelMap("A1")
        w, = b.weights
        algebra = GeneratorAlgebra(b, {"q": w ** 2})
        oracle = self.oracle("A", 2)
        for degree in range(7):
            report = self.assert_independent_degree(algebra, oracle,
                                                    (oracle.omega(1) ** 2,), degree,
                                                    "A1 invariant zero-image generator")
            self.assertEqual(report.rank, int(degree == 0))
            self.assertEqual(report.domain_dimension, int(degree % 2 == 0))
            count("zero_image_and_empty_domain")
        self.assertEqual(algebra.first_relation_degree(1), None)
        self.assertEqual(algebra.first_relation_degree(2), 2)
        self.assertEqual(algebra.degree(1).matrix.dimensions(), (0, 0))

    def test_validation_strict_degrees_and_immutable_reports(self):
        b = BorelMap("A2", omitted_nodes=1)
        u, v = b.weights
        for generators in ({}, {"z": 0}, {"z": 1}, {"z": u + u ** 2}, {"z": v}):
            with self.subTest(generators=generators):
                with self.assertRaises((ValueError, TypeError)):
                    GeneratorAlgebra(b, generators)
                count("invalid_generators_rejected")
        algebra = GeneratorAlgebra(b, {"u": u, "v": u ** 2})
        for degree in (-1, True, 1.5, QQ(3) / 2, "2", None):
            for method in (algebra.monomials, algebra.degree, algebra.through_degree,
                           algebra.hilbert_function, algebra.first_relation_degree):
                with self.subTest(degree=degree, method=method.__name__):
                    with self.assertRaises((ValueError, TypeError)):
                        method(degree)
                    count("strict_degree_validation")
        report = algebra.degree(2, include_ambient=True)
        self.assertTrue(report.matrix.is_immutable())
        self.assertIs(report.matrix, report.coefficient_matrix)
        with self.assertRaises((ValueError, TypeError)):
            report.matrix[0, 0] = 99
        with self.assertRaises((FrozenInstanceError, AttributeError)):
            report.degree = 19
        self.assertIsInstance(report.relations, tuple)
        self.assertIsInstance(report.monomials, tuple)
        self.assertIsInstance(report.pivot_columns, tuple)
        self.assertIsInstance(report.image_basis_monomials, tuple)
        self.assertEqual(len(report.pivot_columns), report.rank)
        self.assertEqual(report.image_basis_monomials,
                         tuple(report.monomials[i] for i in report.pivot_columns))
        serialized = report.as_dict()
        json.dumps(serialized)
        algebra.clear_cache()
        repeated = algebra.degree(2, include_ambient=True)
        self.assertEqual(repeated.matrix, report.matrix)
        self.assertEqual(repeated.relations, report.relations)
        count("immutable_and_cache_regressions", 12)

    def test_two_step_symplectic_generators_and_independent_kernels(self):
        for n, bsize, k in ((3, 1, 1), (4, 1, 1), (5, 1, 1),
                           (4, 1, 2), (4, 2, 1), (5, 2, 2)):
            family = two_step_symplectic(n, bsize, k)
            self.assertIsInstance(family, TwoStepSymplectic)
            a = n - bsize - k
            self.assertEqual((family.a, family.b, family.k), (a, bsize, k))
            self.assertEqual(family.omitted_nodes, (a, a + bsize))
            oracle = self.oracle("C", n)
            middle = tuple(elementary_explicit(oracle.x[a:a + bsize], j)
                           for j in range(1, bsize + 1))
            tail = tuple(elementary_explicit(tuple(x ** 2 for x in oracle.x[-k:]), j)
                         for j in range(1, k + 1))
            expected_middle = tuple("c%d" % j for j in range(1, bsize + 1))
            expected_tail = tuple("e%d" % j for j in range(1, k + 1))
            self.assertEqual(tuple(family.middle_generators()), expected_middle)
            self.assertEqual(tuple(family.tail_generators()), expected_tail)
            all_actual = tuple(family.middle_generators().values()) + tuple(family.tail_generators().values())
            all_expected = middle + tail
            weight_images = tuple(oracle.omega(j) for j in range(1, n + 1))
            for actual, expected in zip(all_actual, all_expected):
                self.assertEqual(oracle.ring(actual(*weight_images)), expected)
                for node in family.retained_nodes:
                    if node < n:
                        reflected = oracle.swap(expected, node - 1, node)
                    else:
                        reflected = oracle.flip(expected, n - 1)
                    self.assertEqual(reflected, expected)
                    count("two_step_explicit_levi_reflections")
                count("two_step_explicit_coordinate_generators")
            joint = family.joint_algebra()
            self.assertEqual(joint.generator_names, expected_middle + expected_tail)
            self.assertEqual(joint.generator_degrees,
                             tuple(range(1, bsize + 1)) + tuple(2 * j for j in range(1, k + 1)))
            proposed_relations = family.predicted_relations()
            self.assertEqual(len(proposed_relations), bsize + k)
            proposed_degrees = tuple(2 * r for r in range(a + 1, n + 1))
            for relation, expected_degree in zip(proposed_relations, proposed_degrees):
                weighted_degrees = {sum(e * weight for e, weight in zip(alpha, joint.generator_degrees))
                                    for alpha in relation.dict()}
                self.assertEqual(weighted_degrees, {expected_degree})
                substituted = joint.substitute(relation)
                self.assertEqual(oracle.reduce(substituted(*weight_images)), 0)
                count("presentation_relations_independent_groebner")
            # Independent Poincare polynomial of W(C_n)/(S_a x S_b x W(C_k)).
            series_ring = PowerSeriesRing(QQ, "t", default_prec=7)
            t = series_ring.gen()
            c_poincare = lambda rank: prod(sum(t ** j for j in range(2 * i))
                                           for i in range(1, rank + 1))
            symmetric_poincare = lambda size: prod(sum(t ** j for j in range(i))
                                                   for i in range(1, size + 1))
            ambient_series = (c_poincare(n) / (symmetric_poincare(a) *
                                                symmetric_poincare(bsize) * c_poincare(k)))
            for degree in range(7):
                report = self.assert_independent_degree(
                    joint, oracle, all_expected, degree,
                    "C%d two-step b=%d k=%d joint" % (n, bsize, k), include_ambient=True)
                self.assertEqual(report.ambient_dimension, ambient_series[degree])
                count("partial_flag_ambient_poincare_coefficients")
                # Test the entire proposed ideal component against the independent
                # quotient kernel, including all generator multiples in this degree.
                multiples = [relation * multiplier
                             for relation, relation_degree in zip(proposed_relations, proposed_degrees)
                             if relation_degree <= degree
                             for multiplier in joint.monomials(degree - relation_degree)]
                exponents = weighted_exponents(joint.generator_degrees, degree)
                relation_matrix = matrix(QQ, len(multiples), len(exponents), [
                    coefficient_dict(relation).get(alpha, 0) for relation in multiples for alpha in exponents])
                oracle_matrix = oracle.input_matrix(all_expected, exponents)
                self.assertEqual(relation_matrix.row_space(), oracle_matrix.left_kernel())
                self.assertEqual(relation_matrix.rank(), report.kernel_dimension)
                count("whole_presentation_ideal_component")
            self.assertEqual(family.predicted_hilbert_prefix(6), joint.hilbert_function(6))
            count("presentation_hilbert_prefix")
            # The subalgebras are independently checked, not inferred from joint ranks.
            for name, algebra, polynomials in (
                    ("middle", family.middle_algebra(), middle),
                    ("tail", family.tail_algebra(), tail)):
                for degree in (0, 2, 4, 6):
                    report = self.assert_independent_degree(
                        algebra, oracle, polynomials, degree,
                        "C%d two-step b=%d k=%d %s" % (n, bsize, k, name))
                    self.assertLessEqual(report.rank, joint.degree(degree).rank)
                    count("factor_rank_bounded_by_joint")

    def test_two_step_c3_explicit_mixed_relation(self):
        family = TwoStepSymplectic(3, 1, 1)
        algebra = family.joint_algebra()
        c, e = algebra.domain.gens()
        relation = c ** 4 + c ** 2 * e + e ** 2
        self.assertEqual(dict(algebra.evaluate(relation)), {})
        self.assertEqual(dict(algebra.evaluate(e ** 3)), {})
        oracle = self.oracle("C", 3)
        self.assertEqual(oracle.reduce(oracle.x[1] ** 4 + oracle.x[1] ** 2 * oracle.x[2] ** 2 +
                                       oracle.x[2] ** 4), 0)
        self.assertEqual(oracle.reduce(oracle.x[2] ** 6), 0)
        self.assertEqual(algebra.first_relation_degree(3), None)
        self.assertEqual(algebra.first_relation_degree(4), 4)
        count("explicit_two_step_mixed_and_tail_relations", 4)

    def test_two_step_parameter_validation(self):
        for parameters in ((2, 1, 1), (3, 0, 1), (3, 1, 0),
                           (3, -1, 1), (True, 1, 1), (3, True, 1),
                           (3, 1, True), (3.0, 1, 1), (3, QQ(3) / 2, 1), ("3", 1, 1)):
            with self.subTest(parameters=parameters):
                with self.assertRaises((ValueError, TypeError)):
                    TwoStepSymplectic(*parameters)
                count("two_step_parameter_validation")
        family = TwoStepSymplectic(4, 2, 1)
        for name in ("n", "a", "b", "k"):
            with self.assertRaises(AttributeError):
                setattr(family, name, 99)
            count("two_step_readonly_configuration")
        generators = family.joint_generators()
        generators["c1"] = 0
        generators["new"] = 1
        self.assertNotEqual(family.joint_generators()["c1"], 0)
        self.assertNotIn("new", family.joint_generators())
        count("two_step_generator_mapping_copies", 2)


if __name__ == "__main__":
    unittest.main()
