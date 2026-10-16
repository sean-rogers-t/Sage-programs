"""From borel/, run ``PYTHONPATH=src sage -python -m unittest discover -s tests -v``.

Set BOREL_EXHAUSTIVE=1 to compare every divisor product in B2, C2, B3,
C3, and D4 (1,088 independent products).  The normal suite includes all
products in the first four groups and the low-length D4 products.  Optional
BOREL_TEST_REPORT names a JSON file in which to record exact check counts.
"""

import json
import os
from pathlib import Path
import platform
import time
import unittest

from sage.all import CartanType, QQ, WeylGroup, prod
from sage.env import SAGE_VERSION

from sage_borel import BorelMap, Expansion, classical_borel, to_type_a
from tests.oracle import OrthonormalBGG, elementary


COUNTS = {}


def count(name, amount=1):
    COUNTS[name] = COUNTS.get(name, 0) + amount


def normalized(engine, expansion):
    answer = {}
    for w, coefficient in expansion.items():
        if coefficient:
            local = engine.W.from_reduced_word(list(w.reduced_word()))
            answer[tuple(local.reduced_word())] = QQ(coefficient)
    return answer


class BorelLibraryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.oracles = {}
        cls.started = time.time()

    @classmethod
    def tearDownClass(cls):
        report = {"sage_version": SAGE_VERSION, "python_version": platform.python_version(),
                  "independent_checks": COUNTS, "seconds": round(time.time() - cls.started, 3),
                  "exhaustive_divisor_products": os.environ.get("BOREL_EXHAUSTIVE") == "1"}
        print("\nExact verification counts: " + json.dumps(report, sort_keys=True), flush=True)
        if os.environ.get("BOREL_TEST_REPORT"):
            Path(os.environ["BOREL_TEST_REPORT"]).write_text(json.dumps(report, indent=2) + "\n")

    def oracle(self, kind):
        if kind not in self.oracles:
            self.oracles[kind] = OrthonormalBGG(kind)
        return self.oracles[kind]

    def assertExpansion(self, engine, actual, expected):
        self.assertEqual(normalized(engine, actual), normalized(engine, expected))

    def test_a1_constants_and_inhomogeneous_polynomials(self):
        b = BorelMap("A1")
        w, = b.weights
        identity, reflection = b.W.one(), b.W.simple_reflections()[1]
        self.assertEqual(dict(b.expand(0)), {})
        self.assertEqual(dict(b.expand(7)), {identity: 7})
        self.assertEqual(dict(b.expand(w)), {reflection: 1})
        self.assertEqual(dict(b.expand(7 - 3 * w + 5 * w ** 2)),
                         {identity: 7, reflection: -3})
        self.assertEqual(b.coefficient(w, [1]), 1)
        self.assertEqual(b.coefficient(w ** 2, [1]), 0)
        self.assertEqual(b.divided_difference(w ** 3, 1), w ** 2)
        self.assertEqual(b.reflect(w, 1), -w)
        count("a1_boundary_cases", 8)

    def test_reflection_and_divided_difference_in_orthonormal_coordinates(self):
        for kind in ("A3", "B3", "C3", "D4"):
            b, oracle = BorelMap(kind), self.oracle(kind)
            f = sum((i + 1) * w ** 3 for i, w in enumerate(b.weights))
            f += prod(b.weights[:2]) + 5
            expected_f = oracle.from_fundamental(f)
            for node in range(1, b.rank + 1):
                with self.subTest(kind=kind, node=node):
                    self.assertEqual(oracle.from_fundamental(b.reflect(f, node)),
                                     oracle.reflect(expected_f, node))
                    self.assertEqual(oracle.from_fundamental(b.divided_difference(f, node)),
                                     oracle.divided(expected_f, node))
                    count("coordinate_reflection_and_difference", 2)
            self.assertExpansion(b, b.expand(f), oracle.expansion(expected_f))
            count("inhomogeneous_bgg_expansions")

    def test_independent_divisor_products(self):
        exhaustive = os.environ.get("BOREL_EXHAUSTIVE") == "1"
        for kind in ("B2", "C2", "B3", "C3", "D4"):
            b, oracle = BorelMap(kind), self.oracle(kind)
            self.assertEqual(oracle.operate(oracle.top, oracle.w0), 1)
            self.assertEqual(oracle.schubert(oracle.W.one()), 1)
            for w in oracle.W:
                if kind == "D4" and not exhaustive and w.length() > 3:
                    continue
                for node in range(1, b.rank + 1):
                    with self.subTest(kind=kind, word=w.reduced_word(), node=node):
                        expected = oracle.expansion(oracle.omega(node) * oracle.schubert(w),
                                                    w.length() + 1)
                        actual = b.chevalley(node, list(w.reduced_word()))
                        self.assertExpansion(b, actual, expected)
                        count("discriminant_chevalley_" + kind)

    def test_repeated_products_against_independent_bgg(self):
        for kind in ("A3", "B3", "C3", "D4"):
            b, oracle = BorelMap(kind), self.oracle(kind)
            examples = ([], [1], [1, 1], [1, 2, 1], [b.rank, b.rank, 1],
                        [b.rank - 1, b.rank, 2, 1])
            for nodes in examples:
                with self.subTest(kind=kind, nodes=nodes):
                    expected = oracle.expansion(prod(oracle.omega(i) for i in nodes), len(nodes))
                    self.assertExpansion(b, b.divisor_product(nodes), expected)
                    count("iterated_bgg_products")

    def test_classical_elementary_polynomials_and_euler(self):
        for kind in ("A3", "B3", "C3", "D4"):
            b, oracle = BorelMap(kind), self.oracle(kind)
            x = oracle.x
            if kind[0] == "A":
                x = tuple(a - sum(x) / len(x) for a in x)
            for k in range(len(x) + 1):
                tail = x[len(x) - k:] if k else ()
                for degree in range(k + 2):
                    squared = kind[0] != "A"
                    expected_poly = elementary([a ** 2 if squared else a for a in tail], degree)
                    with self.subTest(kind=kind, k=k, degree=degree):
                        self.assertEqual(oracle.from_fundamental(b.tail_polynomial(k, degree)),
                                         expected_poly)
                        expected = oracle.expansion(expected_poly)
                        self.assertExpansion(b, b.tail_elementary(k, degree), expected)
                        count("tail_polynomial_and_expansion", 2)
                for degree in (0, 1, k, k + 1):
                    self.assertExpansion(b, b.head_elementary(k, degree),
                                         oracle.expansion(elementary(x[:k], degree)))
                    count("head_bgg_expansions")
            if kind[0] == "D":
                for k in range(len(x) + 1):
                    tail = x[len(x) - k:] if k else ()
                    euler_poly = prod(tail)
                    self.assertExpansion(b, b.euler(k), oracle.expansion(euler_poly))
                    self.assertExpansion(b, b.expand(b.tail_polynomial(k, k)),
                                         oracle.expansion(euler_poly ** 2))
                    count("d_euler_and_square", 2)

    def test_corrected_type_d_historical_counterexample(self):
        b = BorelMap("D4", omitted_nodes=1)
        expected = {b.W.from_reduced_word([2, 1]): -1}
        self.assertExpansion(b, b.tail_elementary(3, 1), expected)
        # The fork product stays in D: the historical repeated routine called B.
        full, oracle = BorelMap("D4"), self.oracle("D4")
        self.assertExpansion(full, full.divisor_product([4, 4, 3, 1]),
                             oracle.expansion(oracle.omega(4) ** 2 * oracle.omega(3) * oracle.omega(1)))
        self.assertExpansion(full, full.divisor_product([1, 1, 1]),
                             oracle.expansion(oracle.omega(1) ** 3))
        count("legacy_d_regressions", 3)

    def test_type_d_split_middle_and_euler_signs(self):
        b = BorelMap("D4", omitted_nodes=2)
        expected = {b.W.from_reduced_word([1, 2]): 1,
                    b.W.from_reduced_word([3, 2]): -1,
                    b.W.from_reduced_word([4, 2]): -1}
        self.assertExpansion(b, b.tail_elementary(2, 1), expected)
        # These are the two distinct middle classes, with their individual signs.
        self.assertExpansion(b, b.euler(2),
                             {b.W.from_reduced_word([3, 2]): 1,
                              b.W.from_reduced_word([4, 2]): -1})
        d5 = BorelMap("D5", omitted_nodes=2)
        self.assertExpansion(d5, d5.euler(3),
                             {d5.W.from_reduced_word([4, 3, 2]): -1,
                              d5.W.from_reduced_word([5, 3, 2]): 1})
        count("d_typed_middle_and_euler_formulas", 3)

    def test_type_c_one_dimensional_head_sign_regression(self):
        # IG(1,6)=P5: the corrected i=2 expression gives one point in degree4,
        # whereas the literal dissertation sign typo gives coefficient5.
        b = BorelMap("C3", omitted_nodes=1)
        expected = {b.W.from_reduced_word([2, 3, 2, 1]): 1}
        self.assertExpansion(b, b.tail_elementary(2, 2), expected)
        count("c_dissertation_sign_regression")

    def test_arbitrary_parabolic_and_empty_quotient(self):
        b, oracle = BorelMap("A3", omitted_nodes=(1, 3)), self.oracle("A3")
        w1, w2, w3 = b.weights
        self.assertExpansion(b, b.expand((w1 + 2 * w3) ** 3),
                             oracle.expansion((oracle.omega(1) + 2 * oracle.omega(3)) ** 3))
        with self.assertRaises(ValueError):
            b.expand(w2)
        for degree in range(7):
            expected = {tuple(b.W.from_reduced_word(w.reduced_word()).reduced_word())
                        for w in oracle.levels().get(degree, [])
                        if (w * oracle.W.simple_reflections()[2]).length() > w.length()}
            actual = {tuple(w.reduced_word()) for w in b.basis(degree)}
            self.assertEqual(actual, expected)
            count("arbitrary_parabolic_degree_sets")
        point = BorelMap("A2", omitted_nodes=[])
        self.assertEqual(dict(point.expand(4)), {point.W.one(): 4})
        with self.assertRaises(ValueError):
            point.expand(point.weights[0])
        self.assertEqual(list(point.basis(0)), [point.W.one()])
        self.assertEqual(list(point.basis(1)), [])
        count("parabolic_image_and_rejections", 6)

    def test_poincare_duality(self):
        projective = BorelMap("A2", omitted_nodes=2)
        top = projective.W.from_reduced_word([1, 2])
        self.assertEqual(projective.poincare_dual([]), top)
        self.assertEqual(projective.poincare_dual(top), projective.W.one())
        self.assertEqual(projective.poincare_dual([2]), projective.W.simple_reflections()[2])
        self.assertEqual(projective.poincare_dual([2]).parent(), projective.W)
        with self.assertRaises(ValueError):
            projective.poincare_dual([1])
        count("a_poincare_regressions", 5)
        self.assertEqual(projective.longest_representative(), top)
        spin = BorelMap("D5", omitted_nodes=5)
        classes = [w for degree in range(11) for w in spin.basis(degree)]
        self.assertEqual(len(classes), 16)
        for w in classes:
            dual = spin.poincare_dual(w)
            self.assertEqual(dual.parent(), spin.W)
            self.assertEqual(dual.length() + w.length(), 10)
            self.assertEqual(spin.poincare_dual(dual), w)
            self.assertFalse(any(dual.has_descent(i, side="right") for i in (1, 2, 3, 4)))
            count("d_spin_poincare_classes")
        # Confirm the intersection pairing itself, independently of the
        # group-theoretic formula: integrate discriminant-derived representatives.
        oracle = self.oracle("D5")
        top = oracle.W.from_reduced_word(spin.longest_representative().reduced_word())
        for w in classes:
            left = oracle.W.from_reduced_word(w.reduced_word())
            right = oracle.W.from_reduced_word(spin.poincare_dual(w).reduced_word())
            self.assertEqual(oracle.operate(oracle.schubert(left) * oracle.schubert(right), top), 1)
            count("d_spin_discriminant_intersection_pairings")
        middle = [w for w in classes if w.length() == 5]
        for left in middle:
            for right in middle:
                a = oracle.W.from_reduced_word(left.reduced_word())
                z = oracle.W.from_reduced_word(right.reduced_word())
                expected = 1 if right == spin.poincare_dual(left) else 0
                self.assertEqual(oracle.operate(oracle.schubert(a) * oracle.schubert(z), top), expected)
                count("d_spin_middle_pairing_matrix_entries")

    def test_classical_factory_conventions(self):
        for family, n, k, expected_type, omitted in (
                ("A", 4, 2, "A3", 2), ("B", 4, 0, "B4", 4),
                ("C", 4, 2, "C4", 2), ("D", 5, 2, "D5", 3)):
            b = classical_borel(family, n, k)
            reference = BorelMap(expected_type, omitted_nodes=omitted)
            self.assertEqual(b.W.cartan_type(), reference.W.cartan_type())
            self.assertEqual(dict(b.expand(b.weights[omitted - 1])),
                             {b.W.simple_reflections()[omitted]: 1})
            with self.assertRaises(ValueError):
                b.expand(b.weights[(omitted % b.rank)])
            count("classical_factory_conventions", 3)
        for args in (("A", 4, 0), ("A", 4, 4), ("B", 4, 4),
                     ("C", 4, -1), ("D", 5, 1), ("D", 3, 2)):
            with self.subTest(args=args), self.assertRaises((TypeError, ValueError)):
                classical_borel(*args)
            count("classical_factory_invalid_inputs")

    def test_signed_embedding_uses_matrix_action(self):
        # Signed order is (e1,...,en,-en,...,-e1).  Obtain the expected
        # permutation directly from its coordinate matrix, without any words.
        for family, n in (("B", 3), ("C", 3), ("D", 4), ("C", 10)):
            W = WeylGroup([family, n], prefix="sigma")
            samples = ([n], [1, n, n - 1, 2], [n, n - 1, n, 1])
            for word in samples:
                with self.subTest(family=family, n=n, word=word):
                    w = W.from_reduced_word(word)
                    matrix = w.matrix()
                    first_half = []
                    for i in range(n):
                        entries = [(j, matrix[j, i]) for j in range(n) if matrix[j, i]]
                        self.assertEqual(len(entries), 1)
                        j, sign = entries[0]
                        first_half.append(j + 1 if sign == 1 else 2 * n - j)
                    expected = first_half + [2 * n + 1 - a for a in reversed(first_half)]
                    embedded = to_type_a(w)
                    actual = [next(j + 1 for j in range(2 * n) if embedded.matrix()[j, i] == 1)
                              for i in range(2 * n)]
                    self.assertEqual(actual, expected)
                    self.assertEqual(embedded.parent().cartan_type(),
                                     WeylGroup(["A", 2 * n - 1]).cartan_type())
                    count("signed_matrix_embedding_cases")

    def test_relabeled_classical_coordinate_helpers_reject(self):
        kind = CartanType(["B", 2]).relabel({1: 2, 2: 1})
        b = BorelMap(kind)
        # The root-data engine still supports this finite Cartan type.
        self.assertEqual(dict(b.expand(b.weights[0])), {b.W.simple_reflections()[1]: 1})
        # Standard classical x_i formulas have a specified Dynkin numbering.
        with self.assertRaises(ValueError):
            b.classical_coordinates()
        with self.assertRaises(ValueError):
            to_type_a(b.W.simple_reflections()[1])
        count("relabeled_classical_helpers", 3)

    def test_exceptional_types_without_weyl_group_enumeration(self):
        for kind in ("G2", "F4", "E6", "E7", "E8"):
            b = BorelMap(kind)
            for node, weight in enumerate(b.weights, 1):
                self.assertEqual(dict(b.expand(weight)), {b.W.simple_reflections()[node]: 1})
                self.assertEqual(b.coefficient(weight, [node]), 1)
                self.assertEqual(dict(b.divisor_product([node])), dict(b.expand(weight)))
                count("exceptional_linear_classes", 3)
            # A W-invariant quadratic constructed by summing over roots must vanish.
            from sage.all import RootSystem
            lattice = RootSystem(kind).weight_lattice()
            invariant = b.ring.zero()
            for root in lattice.positive_roots():
                linear = sum(QQ(root[i]) * b.weights[i - 1] for i in range(1, b.rank + 1))
                invariant += linear ** 2
            self.assertEqual(dict(b.expand(invariant)), {})
            count("exceptional_quadratic_invariants")

    def test_serialization_and_immutability(self):
        b = BorelMap("C3")
        answer = b.expand(5 + b.weights[0] + b.weights[2] ** 2)
        self.assertEqual(answer.to_terms(), b.expand(5 + b.weights[2] ** 2 + b.weights[0]).to_terms())
        self.assertEqual(json.loads(json.dumps(answer.to_terms())), answer.to_terms())
        with self.assertRaises(TypeError):
            answer[b.W.one()] = 6
        with self.assertRaises(AttributeError):
            answer._borel = None
        absent = b.W.simple_reflections()[2]
        self.assertNotIn(absent, answer)
        self.assertEqual(answer.get(absent, 99), 99)
        self.assertEqual(answer[absent], 0)
        for w in answer:
            self.assertEqual(w.parent(), b.W)
        other = BorelMap("B3")
        with self.assertRaises((TypeError, ValueError)):
            b.chevalley(1, other.W.one())
        partial = BorelMap("C3", omitted_nodes=1)
        self.assertNotEqual(b.expand(1), partial.expand(1))
        with self.assertRaises(ValueError):
            Expansion(partial, {partial.W.simple_reflections()[2]: 1})
        count("typed_serialization_regressions", 11)

    def test_expansion_arithmetic_context_and_cache_reset(self):
        b = BorelMap("C3")
        first, second = b.expand(b.weights[0]), b.expand(b.weights[1])
        self.assertEqual(dict(2 * first - first + second),
                         dict(b.expand(b.weights[0] + b.weights[1])))
        self.assertEqual(dict(first - first), {})
        self.assertEqual(dict(b.multiply_divisor(2 * first + second, 3)),
                         dict(b.expand((2 * b.weights[0] + b.weights[1]) * b.weights[2])))
        partial = BorelMap("C3", omitted_nodes=1).expand(1)
        with self.assertRaises(ValueError):
            first + partial
        product = b.divisor_product([1, 2, 3, 1]).to_terms()
        b.clear_cache()
        self.assertEqual(b.last_stats, {})
        self.assertEqual(b.divisor_product([1, 2, 3, 1]).to_terms(), product)
        count("expansion_arithmetic_and_cache_cases", 6)

    def test_invalid_inputs(self):
        for kind in ("A0", "A2~", [["A", 1], ["A", 1]]):
            with self.subTest(kind=kind), self.assertRaises((TypeError, ValueError, NotImplementedError)):
                BorelMap(kind)
        for nodes in (0, 4, (1, 4), (1.5,), True, False, (QQ(3) / 2,), (1.0,)):
            with self.subTest(nodes=nodes), self.assertRaises((TypeError, ValueError)):
                BorelMap("C3", omitted_nodes=nodes)
        b = BorelMap("C3")
        for node in (0, 4, 1.5, True, False, QQ(3) / 2, 1.0):
            with self.subTest(node=node), self.assertRaises((TypeError, ValueError)):
                b.divided_difference(b.weights[0], node)
        for args in ((-1, 1), (4, 1), (1, -1), (1.5, 1), (1, 1.5)):
            with self.subTest(args=args), self.assertRaises((TypeError, ValueError)):
                b.tail_elementary(*args)
        for degree in (-1, 1.5, QQ(3) / 2, True, False, 1.0):
            with self.subTest(degree=degree), self.assertRaises((TypeError, ValueError)):
                b.basis(degree)
        count("input_validation_cases", 29)


if __name__ == "__main__":
    unittest.main(verbosity=2)
