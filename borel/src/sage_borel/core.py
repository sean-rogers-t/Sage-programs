"""Exact Borel characteristic maps using sparse divided-difference frontiers.

Convention: beta(omega_i) = sigma_(s_i).  A polynomial of degree d has its
coefficient at sigma_w given by the constant term of partial_w(f), with the
rightmost simple reflection of w applied first.  Polynomial degree is half the
cohomological degree.
"""

from collections import defaultdict
from operator import index as integer_index

from sage.all import CartanType, PolynomialRing, QQ, RootSystem, WeylGroup, ZZ

from .expansion import Expansion
from .operations import OperationsMixin


class BorelMap(OperationsMixin):
    """The Borel map for a finite irreducible crystallographic root system.

    ``omitted_nodes=None`` means the full flag variety.  An integer means a
    maximal parabolic; an iterable permits any parabolic.  Retained nodes
    generate its Levi Weyl group, and inputs must be invariant under them.
    An empty iterable therefore gives G/G.  Numbering is Sage's Dynkin numbering.

    This engine never enumerates the whole Weyl group to evaluate a polynomial.
    It follows only nonzero divided differences, merging equal Weyl elements
    reached by different reduced words without adding their values twice.
    """

    def __init__(self, cartan_type, omitted_nodes=None):
        try:
            kind = CartanType(cartan_type)
        except (TypeError, ValueError, KeyError, IndexError) as error:
            raise ValueError("a finite irreducible Cartan type is required") from error
        if not kind.is_finite() or not kind.is_irreducible() or not kind.is_crystallographic():
            raise ValueError("a finite irreducible crystallographic Cartan type is required")
        self.cartan_type = kind
        self.rank = int(kind.rank())
        if self.rank < 1:
            raise ValueError("Cartan type must have positive rank")
        self.nodes = tuple(int(i) for i in kind.index_set())
        self._node_positions = {i: position for position, i in enumerate(self.nodes)}
        if omitted_nodes is None:
            omitted = self.nodes
        elif isinstance(omitted_nodes, (int, type(ZZ.one()))):
            omitted = (self._validate_node(omitted_nodes),)
        else:
            try:
                omitted = tuple(self._validate_node(i) for i in omitted_nodes)
            except TypeError as error:
                raise ValueError("omitted_nodes must be an integer or an iterable of nodes") from error
        if len(set(omitted)) != len(omitted):
            raise ValueError("omitted_nodes contains a repeated Dynkin node")
        self.omitted_nodes = tuple(sorted(omitted))
        self.retained_nodes = tuple(i for i in self.nodes if i not in self.omitted_nodes)
        self.W = WeylGroup(kind, prefix="s")
        self.simple_reflections = self.W.simple_reflections()
        # Explicit number of variables retains the multivariate interface at
        # rank one, so exponent dictionaries have the same shape in every type.
        self.ring = PolynomialRing(QQ, self.rank,
                                   names=["w{}".format(i) for i in self.nodes])
        self.weights = tuple(self.ring.gens())
        lattice = RootSystem(kind).weight_lattice()
        roots = lattice.simple_roots()
        self.simple_roots = {
            i: self.ring(sum(QQ(coefficient) * self.weights[self._node_positions[int(j)]]
                             for j, coefficient in roots[i].monomial_coefficients().items()))
            for i in self.nodes
        }
        self._reflection_images = {}
        for i in self.nodes:
            images = list(self.weights)
            images[self._node_positions[i]] -= self.simple_roots[i]
            self._reflection_images[i] = tuple(images)
        self._basis_levels = [(self.W.one(),)]
        self._word_cache = {self.W.one(): ()}
        self._parabolic_longest = None
        self._longest = None
        self.last_stats = {}

    def _validate_node(self, node):
        if isinstance(node, bool):
            raise ValueError("Dynkin node must be an integer, not a boolean")
        try:
            value = integer_index(node)
        except (TypeError, ValueError) as error:
            raise ValueError("Dynkin node must be an integer") from error
        if value not in self.nodes:
            raise ValueError("Dynkin node {} is not in {}".format(value, self.nodes))
        return int(value)

    def _polynomial(self, polynomial):
        try:
            return self.ring(polynomial)
        except (TypeError, ValueError) as error:
            raise TypeError("input must be a polynomial over QQ in {}".format(self.weights)) from error

    def _coerce_weyl(self, word):
        if hasattr(word, "parent") and hasattr(word, "reduced_word"):
            if word.parent() is not self.W:
                raise TypeError("Weyl element has a different parent; supply its Dynkin word explicitly")
            return word
        try:
            indices = tuple(self._validate_node(i) for i in word)
        except TypeError as error:
            raise TypeError("supply a Weyl element or an iterable of Dynkin nodes") from error
        element = self.W.one()
        for i in indices:
            element = element * self.simple_reflections[i]
        return element

    def _word(self, element):
        if element not in self._word_cache:
            self._word_cache[element] = tuple(int(i) for i in element.reduced_word())
        return self._word_cache[element]

    def _expansion(self, coefficients=()):
        return Expansion(self, coefficients)

    def reflect(self, polynomial, node):
        """Apply a simple reflection to a fundamental-weight polynomial."""
        i = self._validate_node(node)
        f = self._polynomial(polynomial)
        return self.ring(f(*self._reflection_images[i]))

    def divided_difference(self, polynomial, node):
        """Return (f - s_i f)/alpha_i, checked by exact polynomial division."""
        i = self._validate_node(node)
        f = self._polynomial(polynomial)
        if not f:
            return self.ring.zero()
        # s_i changes omega_i only.  Most nodes therefore cost no polynomial
        # substitution when a frontier value has small weight support.
        if f.degree(self.weights[self._node_positions[i]]) == 0:
            return self.ring.zero()
        numerator = f - self.reflect(f, i)
        if not numerator:
            return self.ring.zero()
        quotient, remainder = numerator.quo_rem(self.simple_roots[i])
        if remainder:
            raise ArithmeticError("divided difference failed exact polynomial divisibility")
        return self.ring(quotient)

    def _validate_invariance(self, f):
        noninvariant = [i for i in self.retained_nodes if self.reflect(f, i) != f]
        if noninvariant:
            raise ValueError("input is not invariant under retained Levi nodes {}".format(noninvariant))

    def _is_minimal(self, element):
        return not any(element.has_descent(i, side="right") for i in self.retained_nodes)

    def basis(self, degree):
        """Return minimal right-coset representatives in one polynomial degree.

        This optional basis routine uses cached breadth-first generation; expand
        does not call it.  Asking for a large complete basis can still be costly.
        """
        if isinstance(degree, bool):
            raise ValueError("degree must be a nonnegative integer, not a boolean")
        try:
            d = integer_index(degree)
        except (TypeError, ValueError) as error:
            raise ValueError("degree must be a nonnegative integer") from error
        if d < 0:
            raise ValueError("degree must be a nonnegative integer")
        d = int(d)
        while len(self._basis_levels) <= d:
            new_degree = len(self._basis_levels)
            candidates = set()
            for w in self._basis_levels[-1]:
                for i in self.nodes:
                    if w.has_descent(i, side="left"):
                        continue
                    v = self.simple_reflections[i] * w
                    if self._is_minimal(v):
                        candidates.add(v)
            self._basis_levels.append(tuple(sorted(candidates, key=self._word)))
        return self._basis_levels[d]

    def coefficient(self, polynomial, word):
        """Compute one Schubert coefficient without generating other classes."""
        f = self._polynomial(polynomial)
        self._validate_invariance(f)
        element = self._coerce_weyl(word)
        if not self._is_minimal(element):
            raise ValueError("Weyl element is not a minimal right-coset representative")
        for i in reversed(self._word(element)):
            f = self.divided_difference(f, i)
            if not f:
                return QQ.zero()
        return QQ(f(*([0] * self.rank)))

    def expand(self, polynomial):
        """Evaluate the Borel map and return an exact ``Expansion``.

        Homogeneous components are computed separately.  The frontier at depth
        j contains partial_w(f) for the distinct reachable length-j elements.
        Braid-equivalent paths are checked for equality and counted only once.
        ``last_stats`` reports algorithmic work, without relying on wall time.
        """
        f = self._polynomial(polynomial)
        self._validate_invariance(f)
        stats = {"input_monomials": len(f.monomials()), "homogeneous_degrees": [],
                 "divided_differences": 0, "zero_differences": 0,
                 "duplicate_paths": 0, "frontier_sizes": {}, "max_frontier": 0,
                 "output_terms": 0}
        components = defaultdict(dict)
        for exponent, value in f.dict().items():
            components[int(sum(exponent))][exponent] = value
        output = {}
        for degree in sorted(components):
            stats["homogeneous_degrees"].append(degree)
            component = self.ring(components[degree])
            frontier = {self.W.one(): component}
            sizes = [1]
            stats["max_frontier"] = max(stats["max_frontier"], 1)
            for depth in range(degree):
                following = {}
                for w, value in frontier.items():
                    for i in self.nodes:
                        if w.has_descent(i, side="left"):
                            continue
                        difference = self.divided_difference(value, i)
                        stats["divided_differences"] += 1
                        if not difference:
                            stats["zero_differences"] += 1
                            continue
                        v = self.simple_reflections[i] * w
                        if not self._is_minimal(v):
                            continue
                        if v in following:
                            stats["duplicate_paths"] += 1
                            if following[v] != difference:
                                raise ArithmeticError("braid-equivalent divided differences disagree")
                        else:
                            following[v] = difference
                frontier = following
                sizes.append(len(frontier))
                stats["max_frontier"] = max(stats["max_frontier"], len(frontier))
                if not frontier:
                    # All longer divided differences of this component vanish.
                    break
            if len(sizes) == degree + 1:
                for w, value in frontier.items():
                    if value and value.total_degree() != 0:
                        raise ArithmeticError("a degree-matching coefficient is not constant")
                    if value:
                        output[w] = output.get(w, QQ.zero()) + QQ(value)
            stats["frontier_sizes"][str(degree)] = sizes
        result = self._expansion(output)
        stats["output_terms"] = len(result)
        self.last_stats = stats
        return result

    def poincare_dual(self, word):
        """Return the Schubert index w_0*w*w_(0,P) in the same G/P.

        The longest Levi element is found by successive length-increasing
        simple reflections in retained nodes, without enumerating W_P.
        """
        element = self._coerce_weyl(word)
        if not self._is_minimal(element):
            raise ValueError("Poincare duality requires a minimal right-coset representative")
        if self._parabolic_longest is None:
            longest = self.W.one()
            while True:
                for i in self.retained_nodes:
                    if not longest.has_descent(i, side="right"):
                        longest = longest * self.simple_reflections[i]
                        break
                else:
                    break
            self._parabolic_longest = longest
        if self._longest is None:
            self._longest = self.W.long_element()
        result = self._longest * element * self._parabolic_longest
        if not self._is_minimal(result):
            raise ArithmeticError("Poincare-dual index failed the minimal representative check")
        return result

    def longest_representative(self):
        """Return the top Schubert index w_0*w_(0,P), without enumerating W."""
        return self.poincare_dual(self.W.one())

    def clear_cache(self):
        """Release generated bases, reduced words, and per-class product caches."""
        self._basis_levels = [(self.W.one(),)]
        self._word_cache = {self.W.one(): ()}
        self._clear_operation_caches()
        self.last_stats = {}
