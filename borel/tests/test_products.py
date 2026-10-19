"""Independent checks of Schubert multiplication and inverse Borel sections.

Run from ``borel/`` with
``PYTHONPATH=src sage -python -m unittest tests.test_products -v``.
Expected rank-two products come from independently reflected Cartesian BGG
representatives, including an explicit two-variable G2 realization.  The
Grassmannian checks use the complete Gr(2,4) Littlewood--Richardson table and an
independent Cartesian invariant ideal.  Set ``BOREL_PRODUCTS_REPORT`` to save
the exact check counts and Sage version as JSON.
"""

import json
import os
from pathlib import Path
import time
import unittest
from unittest.mock import patch

from sage.all import CartanType, PolynomialRing, QQ, WeylGroup, prod
from sage.env import SAGE_VERSION

from sage_borel import BorelMap, Expansion
from tests.oracle import OrthonormalBGG


COUNTS = {}


def count(name, amount=1):
    COUNTS[name] = COUNTS.get(name, 0) + amount


def normalized(borel, expansion):
    return {tuple(borel.W.from_reduced_word(list(w.reduced_word())).reduced_word()): QQ(c)
            for w, c in expansion.items() if c}


class ExplicitG2BGG:
    """G2 oracle using six explicit roots and two explicit substitutions.

    The coordinates are fundamental weights u,v.  The short simple root is
    2u-v and the long root is -3u+2v.  This oracle calls no library polynomial,
    reflection, divided-difference, lift, or multiplication methods.
    """

    def __init__(self):
        self.W = WeylGroup("G2", prefix="explicit_g2_s")
        self.w0 = self.W.long_element()
        self.ring = PolynomialRing(QQ, 2, names=("u", "v"))
        self.u, self.v = self.ring.gens()
        u, v = self.u, self.v
        self.roots = (2 * u - v, -3 * u + 2 * v, v - u,
                      u, 3 * u - v, v)
        self.top = self.ring(prod(self.roots) / 12)
        self._schubert = {}
        self._levels = {}
        for w in self.W:
            self._levels.setdefault(w.length(), []).append(w)

    def divided(self, f, node):
        if node == 1:
            reflected, alpha = f(self.v - self.u, self.v), 2 * self.u - self.v
        elif node == 2:
            reflected, alpha = f(self.u, 3 * self.u - self.v), -3 * self.u + 2 * self.v
        else:
            raise ValueError("G2 has nodes 1 and 2")
        quotient, remainder = (f - reflected).quo_rem(alpha)
        if remainder:
            raise ArithmeticError("nonexact G2 oracle divided difference")
        return quotient

    def operate(self, f, w):
        f = self.ring(f)
        for node in reversed(w.reduced_word()):
            f = self.divided(f, node)
        return f

    def schubert(self, w):
        if w not in self._schubert:
            self._schubert[w] = self.operate(self.top, w.inverse() * self.w0)
        return self._schubert[w]

    def expansion(self, f, degree=None):
        f = self.ring(f)
        degrees = [degree] if degree is not None else range(max(0, f.degree()) + 1)
        answer = {}
        for d in degrees:
            for w in self._levels.get(d, ()):
                value = QQ(self.operate(f, w)(0, 0))
                if value:
                    answer[w] = value
        return answer

    def from_fundamental(self, polynomial):
        return self.ring(polynomial(self.u, self.v))


class SchubertProductTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.started = time.time()
        cls.oracles = {}
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
            "command": "PYTHONPATH=src sage -python -m unittest tests.test_products -v",
            "independent_checks": COUNTS,
            "expected_values_source": (
                "Explicit Cartesian BGG discriminants in A/B/C/D; explicit G2 "
                "roots and reflections; complete Gr(2,4) Littlewood--Richardson "
                "table and Cartesian invariant quotient; Poincare duality"),
        }
        print("\nSchubert product verification: " + json.dumps(report, sort_keys=True), flush=True)
        if os.environ.get("BOREL_PRODUCTS_REPORT"):
            Path(os.environ["BOREL_PRODUCTS_REPORT"]).write_text(
                json.dumps(report, indent=2) + "\n")

    def oracle(self, kind):
        if kind not in self.oracles:
            self.oracles[kind] = ExplicitG2BGG() if kind == "G2" else OrthonormalBGG(kind)
        return self.oracles[kind]

    def assertExpansion(self, borel, actual, expected):
        self.assertIsInstance(actual, Expansion)
        self.assertEqual(normalized(borel, actual), normalized(borel, expected))

    def all_basis(self, borel):
        return [w for d in range(borel.longest_representative().length() + 1)
                for w in borel.basis(d)]

    def test_all_rank_two_products_against_independent_discriminants(self):
        for kind in ("A2", "B2", "C2", "G2"):
            borel, oracle = BorelMap(kind), self.oracle(kind)
            self.assertEqual(oracle.operate(oracle.top, oracle.w0), 1)
            for u in oracle.W:
                for v in oracle.W:
                    with self.subTest(kind=kind, left=u.reduced_word(), right=v.reduced_word()):
                        expected = oracle.expansion(oracle.schubert(u) * oracle.schubert(v),
                                                    u.length() + v.length())
                        actual = borel.schubert_product(u.reduced_word(), v.reduced_word())
                        self.assertExpansion(borel, actual, expected)
                        count("all_ordered_rank_two_products_" + kind)

    def test_small_schubert_lifts_and_levi_invariance(self):
        cases = (("A2", None), ("B2", None), ("C2", None), ("G2", None),
                 ("A3", (1, 3)), ("C2", 1), ("D4", 1))
        for kind, omitted in cases:
            borel, oracle = BorelMap(kind, omitted_nodes=omitted), self.oracle(kind)
            for w in self.all_basis(borel):
                with self.subTest(kind=kind, omitted=omitted, word=w.reduced_word()):
                    lift = borel.schubert_polynomial(w)
                    self.assertEqual(lift.parent(), borel.ring)
                    self.assertEqual(lift.total_degree(), w.length())
                    for node in borel.retained_nodes:
                        self.assertEqual(borel.reflect(lift, node), lift)
                        count("literal_levi_invariance")
                    expected = {oracle.W.from_reduced_word(w.reduced_word()): 1}
                    self.assertEqual(normalized(borel, oracle.expansion(
                        oracle.from_fundamental(lift))), normalized(borel, expected))
                    self.assertEqual(borel.expand(lift), borel.schubert(w))
                    self.assertEqual(borel.expand(borel.lift(borel.schubert(w))),
                                     borel.schubert(w))
                    count("independent_schubert_lifts_" + kind)
            for node in borel.omitted_nodes:
                self.assertEqual(borel.schubert_polynomial([node]),
                                 borel.weights[borel.nodes.index(node)])
                count("exact_fundamental_weight_lifts")

    def grassmannian_classes(self):
        borel = BorelMap("A3", omitted_nodes=2)
        w1, w2, _ = borel.weights
        c1, c2 = w2, w1 * (w2 - w1)
        polynomials = {(): borel.ring.one(), (1,): c1, (2,): c1 ** 2 - c2,
                       (1, 1): c2, (2, 1): c1 * c2, (2, 2): c2 ** 2}
        classes = {partition: borel.expand(f) for partition, f in polynomials.items()}
        for value in classes.values():
            self.assertEqual(len(value), 1)
            self.assertEqual(next(iter(value.values())), 1)
        return borel, polynomials, classes

    def test_complete_grassmannian_littlewood_richardson_table(self):
        borel, polynomials, classes = self.grassmannian_classes()
        partitions = list(classes)
        # The complete non-unit table inside the 2-by-2 rectangle.  Every
        # omitted ordered pair is zero; reversed pairs are filled explicitly.
        table = {
            ((1,), (1,)): ((2,), (1, 1)),
            ((1,), (2,)): ((2, 1),),
            ((1,), (1, 1)): ((2, 1),),
            ((1,), (2, 1)): ((2, 2),),
            ((2,), (2,)): ((2, 2),),
            ((1, 1), (1, 1)): ((2, 2),),
        }
        table.update({(right, left): value for (left, right), value in list(table.items())})
        zero = borel.expand(0)
        for left in partitions:
            for right in partitions:
                expected = (classes[right] if not left else classes[left] if not right else
                            sum((classes[p] for p in table.get((left, right), ())), zero))
                with self.subTest(left=left, right=right):
                    self.assertEqual(borel.multiply(classes[left], classes[right]), expected)
                    self.assertEqual(classes[left] * classes[right], expected)
                    self.assertEqual(borel.schubert_product(next(iter(classes[left])),
                                                          next(iter(classes[right]))), expected)
                    count("complete_Gr24_LR_products")

    def test_grassmannian_lifts_in_independent_cartesian_quotient(self):
        borel, polynomials, classes = self.grassmannian_classes()
        oracle = self.oracle("A3")
        x = oracle.x
        coefficients = [oracle.ring.one()] + [oracle.ring.zero()] * 4
        for value in x:
            for i in range(4, 0, -1):
                coefficients[i] += value * coefficients[i - 1]
        invariants = coefficients[1:]
        quotient_basis = oracle.ring.ideal(invariants).groebner_basis()
        for partition, polynomial in polynomials.items():
            actual = oracle.from_fundamental(borel.lift(classes[partition]))
            expected = oracle.from_fundamental(polynomial)
            self.assertEqual((actual - expected).reduce(quotient_basis), 0)
            count("Gr24_cartesian_lift_residues")
        for left in polynomials:
            for right in polynomials:
                result = borel.multiply(classes[left], classes[right])
                actual = oracle.from_fundamental(borel.lift(result))
                expected = oracle.from_fundamental(polynomials[left] * polynomials[right])
                self.assertEqual((actual - expected).reduce(quotient_basis), 0)
                count("Gr24_cartesian_product_residues")

    def test_Gr25_nondivisor_products_below_top_degree(self):
        borel = BorelMap("A4", omitted_nodes=2)
        w1, w2, _, _ = borel.weights
        c1, c2 = w2, w1 * (w2 - w1)
        polynomials = {"2": c1 ** 2 - c2, "11": c2,
                       "31": c1 ** 2 * c2 - c2 ** 2, "22": c2 ** 2}
        classes = {partition: borel.expand(f) for partition, f in polynomials.items()}
        for value in classes.values():
            self.assertEqual(len(value), 1)
            self.assertEqual(next(iter(value.values())), 1)
        self.assertEqual(borel.longest_representative().length(), 6)
        self.assertEqual(next(iter(classes["31"])).length(), 4)
        # All three results lie strictly below top degree, so complementary
        # Poincare duality cannot compute these nondivisor products for us.
        cases = (("2", "2", classes["31"] + classes["22"]),
                 ("2", "11", classes["31"]),
                 ("11", "11", classes["22"]))
        for left, right, expected in cases:
            with self.subTest(left=left, right=right):
                self.assertEqual(classes[left] * classes[right], expected)
                self.assertEqual(borel.schubert_product(next(iter(classes[left])),
                                                      next(iter(classes[right]))), expected)
                count("Gr25_below_top_nondivisor_LR_products")
        divisor = borel.expand(c1)
        self.assertEqual(borel.integrate(divisor ** 6), 5)
        self.assertEqual(borel.intersection_number(*([divisor] * 6)), 5)
        count("Gr25_degree_six_intersections", 2)

    def test_arbitrary_two_step_flag_all_products(self):
        borel, oracle = BorelMap("A3", omitted_nodes=(1, 3)), self.oracle("A3")
        basis = self.all_basis(borel)
        self.assertEqual(len(basis), 12)
        for u in basis:
            for v in basis:
                left = oracle.W.from_reduced_word(u.reduced_word())
                right = oracle.W.from_reduced_word(v.reduced_word())
                expected = oracle.expansion(oracle.schubert(left) * oracle.schubert(right),
                                            u.length() + v.length())
                self.assertExpansion(borel, borel.schubert_product(u, v), expected)
                count("A3_two_step_all_ordered_products")

    def test_type_d_quadric_middle_classes_and_euler(self):
        borel, oracle = BorelMap("D4", omitted_nodes=1), self.oracle("D4")
        basis = self.all_basis(borel)
        self.assertEqual(len(basis), 8)
        middle = borel.basis(3)
        self.assertEqual(len(middle), 2)
        for left in middle:
            for right in basis:
                u = oracle.W.from_reduced_word(left.reduced_word())
                v = oracle.W.from_reduced_word(right.reduced_word())
                expected = oracle.expansion(oracle.schubert(u) * oracle.schubert(v),
                                            left.length() + right.length())
                self.assertExpansion(borel, borel.schubert_product(left, right), expected)
                count("D4_quadric_middle_products")
        euler = borel.euler(3)
        h3 = borel.expand(borel.weights[0] ** 3)
        point = borel.schubert(borel.longest_representative())
        self.assertEqual(h3 * h3, 2 * point)
        self.assertEqual(euler * euler, -2 * point)
        self.assertEqual(h3 * euler, borel.expand(0))
        self.assertEqual(borel.integrate(euler * euler), -2)
        count("D4_quadric_geometric_intersections", 4)

    def test_poincare_pairings_and_multiple_intersections(self):
        for kind, omitted in (("A2", None), ("B2", None), ("C2", 1),
                              ("G2", 1), ("A3", (1, 3)), ("D4", 1)):
            borel = BorelMap(kind, omitted_nodes=omitted)
            top_degree = borel.longest_representative().length()
            for w in self.all_basis(borel):
                for v in borel.basis(top_degree - w.length()):
                    expected = 1 if v == borel.poincare_dual(w) else 0
                    self.assertEqual(borel.intersection_number(borel.schubert(w),
                                                              borel.schubert(v)), expected)
                    count("complete_complementary_poincare_pairing_entries")
        borel, _, classes = self.grassmannian_classes()
        self.assertEqual(borel.intersection_number(*([classes[(1,)]] * 4)), 2)
        self.assertEqual(borel.intersection_number(classes[(2,)], classes[(1, 1)]), 0)
        self.assertEqual(borel.integrate(3 * classes[(2, 2)] + 7 * classes[()]), 3)
        self.assertEqual(borel.intersection_number(), 0)
        count("Gr24_multiple_and_inhomogeneous_integrals", 4)

    def test_inhomogeneous_linearity_scalars_units_and_associativity(self):
        borel, _, classes = self.grassmannian_classes()
        unit, h, s2, s11 = (classes[p] for p in ((), (1,), (2,), (1, 1)))
        zero = borel.expand(0)
        left, right = 2 * unit - h + QQ(3) / 2 * s2, 3 * h + s11
        expected = 6 * h + 2 * s11 - 3 * (s2 + s11) + QQ(7) / 2 * classes[(2, 1)]
        self.assertEqual(left * right, expected)
        self.assertEqual(borel.multiply(left, right), expected)
        self.assertEqual(borel.expand(borel.lift(left)), left)
        self.assertEqual(borel.lift(2 * left - right), 2 * borel.lift(left) - borel.lift(right))
        for value in (zero, unit, left, right, classes[(2, 2)]):
            self.assertEqual(value * unit, value)
            self.assertEqual(unit * value, value)
            self.assertEqual(value * zero, zero)
            self.assertEqual(zero * value, zero)
            count("unit_and_zero_products", 4)
        for first, second, third in ((h, h, h), (left, h, right), (s2, h, h)):
            self.assertEqual((first * second) * third, first * (second * third))
            count("nontrivial_associativity_checks")
        self.assertEqual(borel.lift(zero), 0)
        self.assertEqual(borel.lift(unit), 1)
        self.assertEqual(QQ(3) / 2 * h, h * (QQ(3) / 2))
        count("linearity_mixed_degree_and_scalar_checks", 7)

    def test_relabeled_cartan_type_lifts_and_products(self):
        relabel = {1: 3, 2: 5}
        borel = BorelMap(CartanType(["B", 2]).relabel(relabel))
        oracle = self.oracle("B2")
        def relabeled_expected(expansion):
            return {borel.W.from_reduced_word([relabel[i] for i in w.reduced_word()]): c
                    for w, c in expansion.items()}
        for w in self.all_basis(borel):
            lift = borel.schubert_polynomial(w)
            self.assertEqual(normalized(borel, relabeled_expected(oracle.expansion(
                oracle.from_fundamental(lift)))), normalized(borel, {w: 1}))
            self.assertEqual(borel.expand(lift), borel.schubert(w))
            count("relabeled_B2_independent_lifts")
        for u in oracle.W:
            for v in oracle.W:
                expected = relabeled_expected(oracle.expansion(
                    oracle.schubert(u) * oracle.schubert(v), u.length() + v.length()))
                actual = borel.schubert_product([relabel[i] for i in u.reduced_word()],
                                                [relabel[i] for i in v.reduced_word()])
                self.assertExpansion(borel, actual, expected)
                count("relabeled_B2_all_ordered_products")

    def test_exceptional_divisors_agree_with_root_coroot_chevalley(self):
        for kind, node in (("F4", 4), ("E6", 1), ("E7", 7), ("E8", 8)):
            borel = BorelMap(kind, omitted_nodes=node)
            divisor = borel.schubert([node])
            self.assertEqual(borel.schubert_polynomial([node]),
                             borel.weights[borel.nodes.index(node)])
            self.assertEqual(divisor * divisor, borel.chevalley(node, [node]))
            self.assertEqual(borel.schubert_product([node], [node]),
                             borel.chevalley(node, [node]))
            count("exceptional_divisor_lifts_and_products", 3)

    def test_low_degree_E7_inverse_avoids_global_discriminant(self):
        for omitted in (None, 7):
            borel = BorelMap("E7", omitted_nodes=omitted)
            basis = borel.basis(2)
            self.assertEqual(len(basis), 27 if omitted is None else 1)
            # Constructing the 63-root E7 discriminant is unnecessary here.
            # Forbid that computational path, rather than imposing an unstable
            # wall-clock threshold on exact rational linear algebra.
            with patch.object(BorelMap, "_bgg_inverse_polynomial",
                              side_effect=AssertionError("global discriminant path used")):
                for word in basis:
                    polynomial = borel.schubert_polynomial(word)
                    self.assertEqual(polynomial.total_degree(), 2)
                    self.assertEqual(borel.expand(polynomial), borel.schubert(word))
                    for node in borel.retained_nodes:
                        self.assertEqual(borel.reflect(polynomial, node), polynomial)
                    count("E7_low_degree_lifts_without_discriminant")
                mixed = borel.schubert([]) + sum(
                    ((i + 1) * borel.schubert(word) for i, word in enumerate(basis)),
                    borel.expand(0))
                self.assertEqual(borel.expand(borel.lift(mixed)), mixed)
                count("E7_full_and_partial_mixed_lift_roundtrips")
            self.assertIsNone(borel._inverse_discriminant)
            count("E7_global_discriminant_not_constructed")

    def test_context_validation_and_invalid_schubert_indices(self):
        full, partial = BorelMap("A2"), BorelMap("A2", omitted_nodes=1)
        foreign = partial.schubert([])
        for operation in (lambda: full.multiply(full.schubert([]), foreign),
                          lambda: full.schubert([]) * foreign,
                          lambda: full.lift(foreign), lambda: full.integrate(foreign),
                          lambda: full.intersection_number(full.schubert([]), foreign)):
            with self.assertRaises(ValueError):
                operation()
            count("foreign_context_rejections")
        for operation in (lambda: full.multiply(1, full.schubert([])),
                          lambda: full.multiply(full.schubert([]), 1),
                          lambda: full.lift(full.weights[0]),
                          lambda: full.integrate(1)):
            with self.assertRaises(TypeError):
                operation()
            count("non_expansion_input_rejections")
        for operation in (lambda: partial.schubert([2]),
                          lambda: partial.schubert_polynomial([2]),
                          lambda: partial.schubert_product([1], [2]),
                          lambda: full.schubert([7]),
                          lambda: full.schubert_polynomial([True])):
            with self.assertRaises((TypeError, ValueError)):
                operation()
            count("invalid_schubert_index_rejections")
        with self.assertRaises(TypeError):
            full.schubert(self.oracle("A2").W.one())
        count("foreign_weyl_parent_rejections")

    def test_point_rank_one_and_cache_reset(self):
        point = BorelMap("A2", omitted_nodes=[])
        unit = point.schubert([])
        self.assertEqual(point.schubert_polynomial([]), 1)
        self.assertEqual(unit * unit, unit)
        self.assertEqual(point.integrate(7 * unit), 7)
        self.assertEqual(point.intersection_number(), 1)
        self.assertEqual(point.intersection_number(2 * unit, 3 * unit), 6)
        line = BorelMap("A1")
        h = line.schubert([1])
        self.assertEqual(line.schubert_polynomial([1]), line.weights[0])
        self.assertEqual(h * h, line.expand(0))
        self.assertEqual(line.intersection_number(h), 1)
        borel, _, classes = self.grassmannian_classes()
        value = classes[(1,)] + 2 * classes[(2,)]
        lift, product = borel.lift(value), value * value
        borel.clear_cache()
        self.assertEqual(borel.lift(value), lift)
        self.assertEqual(value * value, product)
        self.assertEqual(borel.expand(borel.lift(value)), value)
        count("point_rank_one_and_cache_boundary_checks", 11)

    def test_powers_and_exponent_validation(self):
        borel, _, classes = self.grassmannian_classes()
        unit, h, point = (classes[p] for p in ((), (1,), (2, 2)))
        zero = borel.expand(0)
        self.assertEqual(h ** 0, unit)
        self.assertEqual(zero ** 0, unit)
        self.assertEqual(h ** 1, h)
        self.assertEqual(h ** 2, classes[(2,)] + classes[(1, 1)])
        self.assertEqual(h ** 3, 2 * classes[(2, 1)])
        self.assertEqual(h ** 4, 2 * point)
        self.assertEqual(h ** 5, zero)
        self.assertEqual(h ** 64, zero)
        self.assertEqual((unit + h) ** 3, unit + 3 * h + 3 * h ** 2 + h ** 3)
        for exponent in (-1, True, QQ(3) / 2, "2"):
            with self.assertRaises((TypeError, ValueError)):
                h ** exponent
        with self.assertRaises((TypeError, ValueError)):
            pow(h, 2, 5)
        count("power_boundary_and_geometric_checks", 14)

    def test_explicit_inverse_methods_and_validation(self):
        for kind, omitted in (("A2", 1), ("G2", 1)):
            borel, oracle = BorelMap(kind, omitted_nodes=omitted), self.oracle(kind)
            for w in self.all_basis(borel):
                for method in ("linear", "bgg", "auto"):
                    lift = borel.schubert_polynomial(w, method=method)
                    expected = {oracle.W.from_reduced_word(w.reduced_word()): 1}
                    self.assertEqual(normalized(borel, oracle.expansion(
                        oracle.from_fundamental(lift))), normalized(borel, expected))
                    self.assertEqual(borel.expand(borel.lift(borel.schubert(w), method=method)),
                                     borel.schubert(w))
                    count("explicit_inverse_method_roundtrips_" + kind)
            mixed = borel.schubert([]) + 2 * borel.schubert([omitted])
            for method in ("linear", "bgg", "auto"):
                self.assertEqual(borel.expand(borel.lift(mixed, method=method)), mixed)
                count("explicit_inverse_method_mixed_roundtrips")
            with self.assertRaises((TypeError, ValueError)):
                borel.schubert_polynomial([], method="unknown")
            with self.assertRaises((TypeError, ValueError)):
                borel.lift(borel.expand(0), method="unknown")
            count("invalid_inverse_method_rejections", 2)


if __name__ == "__main__":
    unittest.main()
