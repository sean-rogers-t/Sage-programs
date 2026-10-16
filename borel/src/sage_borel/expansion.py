"""Immutable, exact Schubert expansions.

The keys are Sage Weyl-group elements in the parent used by a ``BorelMap``.
Coefficients are rational numbers.  No polynomial representatives are hidden in
this class: it records a linear combination of Schubert classes only.
"""

from collections.abc import Mapping
from types import MappingProxyType

from sage.all import QQ


class Expansion(Mapping):
    """An immutable Schubert-basis expansion with its flag-variety context.

    Usually constructed by ``BorelMap.expand`` or ``BorelMap.chevalley``.  Missing
    valid Weyl-group keys have coefficient zero.  ``to_terms`` is suitable for
    JSON serialization and uses a deterministic choice of reduced word.
    """

    __slots__ = ("_borel", "_data", "_keys")

    def __setattr__(self, name, value):
        if hasattr(self, name):
            raise AttributeError("Schubert expansions are immutable")
        object.__setattr__(self, name, value)

    def __init__(self, borel_map, coefficients=()):
        data = {}
        items = coefficients.items() if isinstance(coefficients, Mapping) else coefficients
        for word, coefficient in items:
            element = borel_map._coerce_weyl(word)
            if not borel_map._is_minimal(element):
                raise ValueError("Schubert index is not a minimal right-coset representative")
            value = QQ(coefficient)
            if value:
                data[element] = data.get(element, QQ.zero()) + value
        data = {element: value for element, value in data.items() if value}
        self._borel = borel_map
        self._data = MappingProxyType(data)
        self._keys = tuple(sorted(data, key=lambda w: (w.length(), tuple(w.reduced_word()))))

    @property
    def cartan_type(self):
        return self._borel.cartan_type

    @property
    def omitted_nodes(self):
        return self._borel.omitted_nodes

    @property
    def borel_map(self):
        return self._borel

    def __getitem__(self, word):
        element = self._borel._coerce_weyl(word)
        return self._data.get(element, QQ.zero())

    def __iter__(self):
        return iter(self._keys)

    def __len__(self):
        return len(self._keys)

    def __contains__(self, word):
        # Mapping's default implementation detects missing keys via KeyError;
        # our coefficient indexing intentionally returns zero instead.
        element = self._borel._coerce_weyl(word)
        return element in self._data

    def get(self, word, default=None):
        """Use ordinary Mapping.get semantics; indexing itself defaults to zero."""
        element = self._borel._coerce_weyl(word)
        return self._data.get(element, default)

    def coefficient(self, word):
        """Return a coefficient; accept a Weyl element or a word of nodes."""
        return self[word]

    def as_dict(self):
        """Return a mutable copy, leaving this expansion unchanged."""
        return dict(self._data)

    def to_terms(self):
        """Return exact coefficients and reduced words as JSON-safe values."""
        return [
            {"word": [int(i) for i in w.reduced_word()], "coefficient": str(self._data[w])}
            for w in self._keys
        ]

    def _check_context(self, other):
        if not isinstance(other, Expansion):
            return False
        if (self.cartan_type != other.cartan_type or
                self.omitted_nodes != other.omitted_nodes or
                self._borel.W is not other._borel.W):
            raise ValueError("Schubert expansions belong to different flag-variety contexts")
        return True

    def __eq__(self, other):
        if isinstance(other, Expansion):
            if (self.cartan_type != other.cartan_type or
                    self.omitted_nodes != other.omitted_nodes or
                    self._borel.W is not other._borel.W):
                return False
        return Mapping.__eq__(self, other)

    def __add__(self, other):
        if other == 0:
            return self
        if not self._check_context(other):
            return NotImplemented
        data = self.as_dict()
        for w, value in other.items():
            data[w] = data.get(w, QQ.zero()) + value
        return Expansion(self._borel, data)

    def __radd__(self, other):
        return self.__add__(other)

    def __neg__(self):
        return Expansion(self._borel, ((w, -value) for w, value in self.items()))

    def __sub__(self, other):
        if not isinstance(other, Expansion):
            return NotImplemented
        return self + (-other)

    def __mul__(self, scalar):
        if isinstance(scalar, Expansion):
            return NotImplemented
        try:
            value = QQ(scalar)
        except (TypeError, ValueError):
            return NotImplemented
        return Expansion(self._borel, ((w, value * c) for w, c in self.items()))

    __rmul__ = __mul__

    def __repr__(self):
        if not self:
            return "0"
        return " + ".join(
            "{}*sigma{}".format(self._data[w], tuple(int(i) for i in w.reduced_word()))
            for w in self._keys
        )
