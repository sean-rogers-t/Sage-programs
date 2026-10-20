"""Sparse divisor operations and classical-coordinate generators.

This module never enumerates a Weyl group.  Chevalley multiplication scans its
positive-root reflections once per newly encountered Schubert class and caches
the covers.  Classical elementary polynomials are constructed by coefficient
recurrences rather than enumerating subsets.
"""

from math import comb
from operator import index as integer_index

from sage.all import CartanType, QQ


def _integer(value, name, minimum=0):
    """Accept Python/Sage integers without silently truncating rationals."""
    if isinstance(value, bool):
        raise TypeError("%s must be an integer, not a boolean" % name)
    try:
        result = integer_index(value)
    except TypeError as error:
        raise TypeError("%s must be an integer" % name) from error
    if minimum is not None and result < minimum:
        raise ValueError("%s must be at least %s" % (name, minimum))
    return result


def elementary_polynomial(ring, values, degree):
    """Return ``e_degree(values)`` using O(len(values) * degree) operations.

    The empty elementary polynomial is one and degrees above the alphabet size
    give zero.  ``values`` are coerced into ``ring``; no subsets are formed.
    """
    degree = _integer(degree, "degree")
    values = tuple(ring(value) for value in values)
    if degree > len(values):
        return ring.zero()
    coefficients = [ring.one()] + [ring.zero() for _ in range(degree)]
    for position, value in enumerate(values, 1):
        for j in range(min(position, degree), 0, -1):
            coefficients[j] += value * coefficients[j - 1]
    return coefficients[degree]


def complete_polynomial(ring, values, degree):
    """Return ``h_degree(values)`` by the complete-function recurrence."""
    degree = _integer(degree, "degree")
    coefficients = [ring.one()] + [ring.zero() for _ in range(degree)]
    for value in values:
        value = ring(value)
        # Ascending order permits repeated occurrences of the current variable.
        for j in range(1, degree + 1):
            coefficients[j] += value * coefficients[j - 1]
    return coefficients[degree]


class OperationsMixin:
    """Additional operations for :class:`sage_borel.core.BorelMap`."""

    def _operation_node(self, node):
        node = _integer(node, "simple-root index", minimum=None)
        if node not in self.cartan_type.index_set():
            raise ValueError("simple-root index %s is not in %s" %
                             (node, self.cartan_type.index_set()))
        return node

    def _operation_divisor_node(self, node):
        node = self._operation_node(node)
        omitted = self.omitted_nodes
        if omitted is not None and node not in omitted:
            raise ValueError(
                "omega_%s is not a divisor on this parabolic quotient; "
                "use an omitted simple-root index" % node)
        return node

    def _chevalley_reflections(self):
        if not hasattr(self, "_operation_reflections"):
            reflections = self.W.reflections()
            fundamental = self.W.domain().fundamental_weights()
            nodes = tuple(self.cartan_type.index_set())
            self._operation_node_positions = {
                node: position for position, node in enumerate(nodes)}
            # In all finite types, <omega_i, alpha^vee> is computed in Sage's
            # ambient root realization.  In particular, root lengths and the
            # type D fork are handled by associated_coroot(), not guessed.
            data = []
            for root in reflections.keys():
                coroot = root.associated_coroot()
                coefficients = tuple(QQ(fundamental[node].dot_product(coroot))
                                     for node in nodes)
                data.append((reflections[root], coefficients))
            self._operation_reflections = tuple(data)
            self._operation_covers = {}
            self._operation_chevalley = {}
        return self._operation_reflections

    def _chevalley_covers(self, element):
        reflections = self._chevalley_reflections()
        if element not in self._operation_covers:
            target_length = element.length() + 1
            self._operation_covers[element] = tuple(
                (target, coefficients)
                for reflection, coefficients in reflections
                for target in (element * reflection,)
                if target.length() == target_length)
        return self._operation_covers[element]

    def chevalley(self, node, word_or_element):
        r"""Multiply ``epsilon_w`` by ``beta(omega_node)``.

        ``word_or_element`` is a Weyl-group element belonging to this map, or a
        sequence of simple-root indices in multiplication order.  The result
        is a sparse :class:`Expansion`.  Only omitted nodes define divisors on
        a parabolic quotient.
        """
        node = self._operation_divisor_node(node)
        element = self._coerce_weyl(word_or_element)
        if not self._is_minimal(element):
            raise ValueError("Weyl element is not a minimal right-coset representative")
        self._chevalley_reflections()
        key = (node, element)
        if key not in self._operation_chevalley:
            position = self._operation_node_positions[node]
            terms = {
                target: coefficients[position]
                for target, coefficients in self._chevalley_covers(element)
                if coefficients[position]}
            self._operation_chevalley[key] = self._expansion(terms)
        return self._operation_chevalley[key]

    def multiply_divisor(self, expansion, node):
        """Multiply an expansion by a fundamental-weight divisor sparsely."""
        from .expansion import Expansion

        node = self._operation_divisor_node(node)
        if not isinstance(expansion, Expansion):
            raise TypeError("expansion must be a Schubert Expansion")
        self._expansion()._check_context(expansion)
        terms = {}
        for element, coefficient in expansion.items():
            for target, factor in self.chevalley(node, element).items():
                terms[target] = terms.get(target, QQ.zero()) + coefficient * factor
                if not terms[target]:
                    del terms[target]
        return self._expansion(terms)

    def divisor_product(self, indices):
        """Return ``beta(prod(omega_i for i in indices))``.

        Empty products give the unit; one-factor products give the divisor.
        Coefficients are combined after each multiplication, avoiding the
        historical routines' duplicate-term lists and repeated root scans.
        """
        if isinstance(indices, (str, bytes)):
            raise TypeError("indices must be an iterable of simple-root indices")
        try:
            nodes = tuple(self._operation_divisor_node(node) for node in indices)
        except TypeError as error:
            raise TypeError("indices must be an iterable of simple-root indices") from error
        result = self._expansion({self.W.one(): QQ.one()})
        for node in nodes:
            result = self.multiply_divisor(result, node)
            if not result:
                break
        return result

    def classical_coordinates(self):
        r"""Return the standard coordinate weights ``(x_1,...,x_n)``.

        For ``A_r`` there are ``r+1`` coordinates and their sum is zero; this
        is also the characteristic-map image of the GL(r+1) coordinates, whose
        central character maps to zero.  In B and D the half-weight/fork
        conventions use Bourbaki/Sage numbering exactly.
        """
        if not hasattr(self, "_operation_coordinates"):
            if self.cartan_type != CartanType([self.cartan_type.type(), self.rank]):
                raise ValueError("classical coordinate generators require standard Dynkin numbering")
            kind = self.cartan_type.type()
            weights = self.weights
            zero = self.ring.zero()
            if kind == "A":
                coordinates = tuple(weights[j] - (weights[j - 1] if j else zero)
                                    for j in range(self.rank)) + (-weights[-1],)
            elif kind == "C":
                coordinates = tuple(weights[j] - (weights[j - 1] if j else zero)
                                    for j in range(self.rank))
            elif kind == "B":
                coordinates = tuple(weights[j] - (weights[j - 1] if j else zero)
                                    for j in range(self.rank - 1)) + (
                    2 * weights[-1] - (weights[-2] if self.rank > 1 else zero),)
            elif kind == "D":
                coordinates = tuple(weights[j] - (weights[j - 1] if j else zero)
                                    for j in range(self.rank - 2)) + (
                    weights[-2] + weights[-1] - weights[-3],
                    weights[-1] - weights[-2])
            else:
                raise ValueError("standard coordinate generators require type A, B, C or D")
            self._operation_coordinates = coordinates
        return self._operation_coordinates

    def _classical_count(self, count, name):
        count = _integer(count, name)
        size = len(self.classical_coordinates())
        if count > size:
            raise ValueError("%s must be between 0 and %s" % (name, size))
        return count

    def _classical_squared(self, squared):
        if squared is None:
            return self.cartan_type.type() != "A"
        if not isinstance(squared, bool):
            raise TypeError("squared must be a boolean or None")
        return squared

    def tail_polynomial(self, k, degree, squared=None):
        """Return an elementary polynomial in the last ``k`` coordinates.

        ``degree`` is the elementary index.  Type A defaults to ``e_degree(x)``;
        B/C/D default to ``e_degree(x^2)``.  Thus the latter has polynomial
        degree ``2*degree`` and cohomological degree ``4*degree``.  An explicit
        boolean ``squared`` overrides this default.  This returns the literal
        polynomial, before passage to the Borel quotient.
        """
        k = self._classical_count(k, "k")
        degree = _integer(degree, "degree")
        squared = self._classical_squared(squared)
        coordinates = self.classical_coordinates()
        values = coordinates[len(coordinates) - k:] if k else ()
        if squared:
            values = tuple(value * value for value in values)
        return elementary_polynomial(self.ring, values, degree)

    def _tail_shortcut_permitted(self, head_count, squared):
        """Whether Levi generators preserve the chosen coordinate alphabets."""
        if self.omitted_nodes is None:
            return True
        kind = self.cartan_type.type()
        if not ((kind == "A" and not squared) or
                (kind in ("B", "C", "D") and squared)):
            return False
        forbidden = {head_count} if 0 < head_count < len(self.classical_coordinates()) else set()
        if kind == "D" and head_count == self.rank - 1:
            forbidden.add(self.rank)
        return forbidden.issubset(set(self.omitted_nodes))

    def tail_elementary(self, k, degree, squared=None):
        r"""Expand a classical tail generator, using a quotient shortcut.

        When cheaper, use ``e_d(tail) = (-1)^d h_d(head)`` (type A), or
        ``e_d(tail^2) = (-1)^d h_d(head^2)`` (B/C/D), in the coinvariant ring.
        These equalities follow by inverting the total elementary generating
        function of all coordinates.  They are equalities after applying beta,
        not polynomial identities.  ``tail_polynomial`` always retains the
        unmodified input.  Nonmatching parabolic contexts use the literal input
        so that the engine's invariance validation remains effective.
        """
        k = self._classical_count(k, "k")
        degree = _integer(degree, "degree")
        squared = self._classical_squared(squared)
        if degree > k:
            return self.expand(self.ring.zero())
        if degree == 0:
            return self.expand(self.ring.one())
        coordinates = self.classical_coordinates()
        head_count = len(coordinates) - k
        kind = self.cartan_type.type()
        valid_identity = ((kind == "A" and not squared) or
                          (kind in ("B", "C", "D") and squared))
        complete_cost = comb(head_count + degree - 1, degree) if head_count else 0
        if (valid_identity and complete_cost < comb(k, degree) and
                self._tail_shortcut_permitted(head_count, squared)):
            values = coordinates[:head_count]
            if squared:
                values = tuple(value * value for value in values)
            polynomial = (-1) ** degree * complete_polynomial(self.ring, values, degree)
        else:
            polynomial = self.tail_polynomial(k, degree, squared=squared)
        return self.expand(polynomial)

    def head_polynomial(self, count, degree):
        """Return ``e_degree(x_1,...,x_count)`` before applying beta."""
        count = self._classical_count(count, "count")
        degree = _integer(degree, "degree")
        return elementary_polynomial(self.ring, self.classical_coordinates()[:count], degree)

    def head_elementary(self, count, degree):
        """Expand ``e_degree(x_1,...,x_count)`` in the Schubert basis."""
        return self.expand(self.head_polynomial(count, degree))

    def euler_polynomial(self, k):
        """Return the type D Euler product of the last ``k`` coordinates."""
        if self.cartan_type.type() != "D":
            raise ValueError("the Euler tail generator is provided only for type D")
        k = self._classical_count(k, "k")
        result = self.ring.one()
        coordinates = self.classical_coordinates()
        for value in coordinates[len(coordinates) - k:] if k else ():
            result *= value
        return result

    def euler(self, k):
        """Expand the type D Euler tail product in the Schubert basis."""
        return self.expand(self.euler_polynomial(k))

    def operation_cache_info(self):
        """Return cover/divisor cache sizes without initializing root data."""
        return {
            "positive_root_reflections": len(getattr(self, "_operation_reflections", ())),
            "covered_classes": len(getattr(self, "_operation_covers", {})),
            "divisor_products": len(getattr(self, "_operation_chevalley", {})),
        }

    def _clear_operation_caches(self):
        """Drop potentially large per-class caches; retain fixed root data."""
        if hasattr(self, "_operation_covers"):
            self._operation_covers.clear()
            self._operation_chevalley.clear()
