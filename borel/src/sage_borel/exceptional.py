"""Named exceptional Levi invariants and their Borel characteristic images.

The catalogue stores polynomial inputs and provenance, never precomputed
Schubert answers.  All images are evaluated by :class:`BorelMap` using exact
divided differences.  Coordinates use Sage's standard Dynkin numbering and
the normalized Springer convention ``Theta_i -> omega_i``.
"""

from collections import OrderedDict
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
import json
from operator import index as integer_index
from types import MappingProxyType

from sage.all import CartanType, QQ

from .core import BorelMap


_EXCEPTIONAL_TYPES = ("G2", "F4", "E6", "E7", "E8")


def _freeze(value):
    """Copy JSON metadata into recursively read-only Python values."""
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    return value


def _thaw(value):
    if hasattr(value, "items"):
        return {key: _thaw(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_thaw(item) for item in value]
    return value


@dataclass(frozen=True)
class GeneratorInfo:
    """An immutable catalogue entry; degree is polynomial degree.

    ``source`` identifies the theorem and row; ``normalization`` describes the
    coordinate convention; ``corrections`` records audited changes to printed
    formulas.  ``as_dict()`` returns an independent JSON-compatible copy.
    """

    kind: str
    parabolic: int
    name: str
    aliases: tuple
    polynomial_degree: int
    spec: object
    coordinates: object
    source: object
    normalization: object
    corrections: tuple

    @property
    def cohomological_degree(self):
        return 2 * self.polynomial_degree

    def as_dict(self):
        return {field: _thaw(getattr(self, field))
                for field in self.__dataclass_fields__}


@lru_cache(maxsize=1)
def _catalogue():
    path = files("sage_borel").joinpath("data", "exceptional_generators.json")
    with path.open("r", encoding="utf-8") as stream:
        data = json.load(stream)
    if data.get("schema_version") != 1:
        raise RuntimeError("unsupported exceptional generator catalogue schema")
    entries = {}
    for record in data["generators"]:
        source = dict(data.get("sources", {}).get(record.get("source", {}).get("key"), {}))
        source.update(record.get("source", {}))
        info = GeneratorInfo(
            kind=record["kind"], parabolic=int(record["parabolic"]),
            name=record["name"], aliases=tuple(record.get("aliases", ())),
            polynomial_degree=int(record["polynomial_degree"]),
            spec=_freeze(record["spec"]),
            coordinates=_freeze(record.get("coordinates", {})),
            source=_freeze(source),
            normalization=_freeze(record.get("normalization", {})),
            corrections=_freeze(record.get("corrections", ())))
        key = (info.kind, info.parabolic)
        entries.setdefault(key, []).append(info)
    entries = {key: tuple(value) for key, value in entries.items()}
    for key, infos in entries.items():
        seen = set()
        for info in infos:
            names = (info.name,) + info.aliases
            if seen.intersection(names) or len(set(names)) != len(names):
                raise RuntimeError("duplicate exceptional generator name in {}".format(key))
            seen.update(names)
    return (MappingProxyType(entries), _freeze(data.get("conventions", {})),
            _freeze(data.get("coverage", {})))


def _exceptional_type(cartan_type):
    try:
        kind = CartanType(cartan_type)
        label = "{}{}".format(kind.type(), kind.rank())
    except (TypeError, ValueError, KeyError, IndexError, AttributeError) as error:
        raise ValueError("exceptional type must be G2, F4, E6, E7 or E8") from error
    if label not in _EXCEPTIONAL_TYPES or kind != CartanType(label):
        raise ValueError("exceptional type must be G2, F4, E6, E7 or E8 with standard Dynkin numbering")
    return kind, label


def _parabolic_node(parabolic, rank):
    if isinstance(parabolic, bool):
        raise ValueError("parabolic must be an integer Dynkin node")
    try:
        node = integer_index(parabolic)
    except (TypeError, ValueError) as error:
        raise ValueError("parabolic must be an integer Dynkin node") from error
    if not 1 <= node <= rank:
        raise ValueError("parabolic must lie between 1 and {}".format(rank))
    return int(node)


def _entries(cartan_type, parabolic):
    kind, label = _exceptional_type(cartan_type)
    node = _parabolic_node(parabolic, int(kind.rank()))
    entries, _, _ = _catalogue()
    try:
        return kind, label, node, entries[label, node]
    except KeyError as error:
        available = sorted(p for (t, p) in entries if t == label)
        text = ", ".join("P{}".format(p) for p in available) or "none"
        raise NotImplementedError(
            "the generator catalogue does not yet cover {}/P{}; available {} parabolics: {}. "
            "Use BorelMap directly for other invariant polynomials.".format(label, node, label, text)
        ) from error


def _elementary(values, degree, ring):
    """Dynamic programming avoids enumerating subsets of coordinate values."""
    if degree < 0:
        raise RuntimeError("negative degree in exceptional catalogue")
    if degree > len(values):
        return ring.zero()
    coefficients = [ring.one()] + [ring.zero() for _ in range(degree)]
    for count, value in enumerate(values, 1):
        for index in range(min(count, degree), 0, -1):
            coefficients[index] += value * coefficients[index - 1]
    return coefficients[degree]


class ExceptionalBorel(BorelMap):
    """A maximal-parabolic Borel map with named exceptional generators.

    ``generators()`` discovers canonical names and their precise polynomial
    specifications.  ``generator_polynomial(name)`` returns the invariant in
    fundamental weights, and ``generator_image(name)`` computes its Schubert
    expansion.  Source metadata is available through ``generator_info(name)``.
    """

    def __init__(self, cartan_type, parabolic):
        kind, label, node, entries = _entries(cartan_type, parabolic)
        super().__init__(kind, omitted_nodes=node)
        self.exceptional_type = label
        self.parabolic = node
        self._generator_infos = entries
        self._generator_lookup = {
            name: info for info in entries for name in (info.name,) + info.aliases
        }
        self.catalogue_coverage = _catalogue()[2].get(label, MappingProxyType({}))
        self._generator_polynomials = OrderedDict()

    def generators(self):
        """Return immutable descriptors in published/source order."""
        return self._generator_infos

    def generator_info(self, name):
        """Look up a canonical name or an explicitly recorded alias."""
        if not isinstance(name, str):
            raise TypeError("generator name must be a string")
        try:
            return self._generator_lookup[name]
        except KeyError as error:
            available = ", ".join(info.name for info in self._generator_infos)
            partial = (" This catalogue currently contains only the displayed source generators."
                       if not self.catalogue_coverage.get("complete_generator_lists", False) else "")
            raise KeyError("unknown generator {!r} for {}/P{}; choose from: {}.{}".format(
                name, self.exceptional_type, self.parabolic, available, partial)) from error

    def _coordinate(self, coefficients):
        if len(coefficients) != self.rank:
            raise RuntimeError("coordinate rank differs from exceptional Cartan rank")
        return self.ring(sum(QQ(value) * weight
                             for value, weight in zip(coefficients, self.weights)))

    def _psi_polynomial(self, degree):
        if self.exceptional_type != "E7" or self.parabolic != 7:
            raise RuntimeError("the psi convention requires E7/P7")
        _, conventions, _ = _catalogue()
        definition = conventions["psi_e7"]
        # These are the six A5 coordinates and the orthogonal A1 coordinate
        # used for the 27-weight E6 power sums in Part III, Theorem 8.2(g).
        x_vectors = definition["x_simple_root_coefficients"]
        z_vector = definition["z_simple_root_coefficients"]
        if len(x_vectors) != 6 or any(len(vector) != self.rank for vector in x_vectors) or len(z_vector) != self.rank:
            raise RuntimeError("invalid E7 psi root-coordinate specification")
        roots = tuple(self.simple_roots[i] for i in self.nodes)
        x = [self.ring(sum(QQ(value) * root for value, root in zip(vector, roots)))
             for vector in x_vectors]
        z = self.ring(sum(QQ(value) * root for value, root in zip(z_vector, roots)))
        values = [value + z for value in x] + [value - z for value in x]
        values += [-x[i] - x[j] for i in range(len(x)) for j in range(i + 1, len(x))]
        return self.ring(sum(value ** degree for value in values))

    def _build_generator(self, info):
        spec = info.spec
        operation = spec["operation"]
        coordinates = {int(index): self._coordinate(vector)
                       for index, vector in info.coordinates.items()}
        if operation == "omega":
            result = self.weights[int(spec["index"]) - 1]
        elif operation == "linear":
            result = coordinates[int(spec["variable"])]
        elif operation == "elementary":
            exponent = 2 if spec.get("squared", False) else 1
            values = [coordinates[int(index)] ** exponent for index in spec["variables"]]
            result = _elementary(values, int(spec["degree"]), self.ring)
        elif operation == "psi":
            result = self._psi_polynomial(int(spec["degree"]))
        elif operation == "g2quadratic":
            result = (self.weights[1] - QQ(3) / 2 * self.weights[0]) ** 2
        elif operation == "g2product":
            result = self.weights[0] * self.weights[1] - self.weights[0] ** 2
        elif operation in ("product", "euler"):
            result = self.ring.one()
            for index in spec["variables"]:
                result *= coordinates[int(index)]
        else:
            raise RuntimeError("unsupported exceptional generator operation {!r}".format(operation))
        result = self.ring(result)
        if result and result.total_degree() != info.polynomial_degree:
            raise RuntimeError("generator degree differs from its catalogue specification")
        if any(self.reflect(result, node) != result for node in self.retained_nodes):
            raise ArithmeticError("catalogued generator {!r} is not Levi invariant".format(info.name))
        return result

    def generator_polynomial(self, name):
        """Construct the exact named invariant in ``self.ring``.

        At most 16 expanded polynomials are cached per object.  Schubert
        coefficients are never read from the catalogue or cached here.
        """
        info = self.generator_info(name)
        if info.name in self._generator_polynomials:
            polynomial = self._generator_polynomials.pop(info.name)
        else:
            polynomial = self._build_generator(info)
        self._generator_polynomials[info.name] = polynomial
        if len(self._generator_polynomials) > 16:
            self._generator_polynomials.popitem(last=False)
        return polynomial

    def generator_image(self, name):
        """Compute a named generator's Schubert expansion using the core."""
        return self.expand(self.generator_polynomial(name))

    def clear_cache(self):
        """Clear operation/frontier caches and expanded generator inputs."""
        super().clear_cache()
        self._generator_polynomials.clear()


def exceptional_borel(cartan_type, parabolic):
    """Create a reusable exceptional maximal-parabolic generator context."""
    return ExceptionalBorel(cartan_type, parabolic)


def exceptional_generators(cartan_type, parabolic):
    """Discover the named generators without constructing a Weyl group."""
    return _entries(cartan_type, parabolic)[3]


@lru_cache(maxsize=16)
def _cached_borel(kind, parabolic):
    return ExceptionalBorel(kind, parabolic)


def _convenience_context(cartan_type, parabolic):
    kind, label = _exceptional_type(cartan_type)
    node = _parabolic_node(parabolic, int(kind.rank()))
    return _cached_borel(label, node)


def exceptional_polynomial(cartan_type, parabolic, generator):
    """Return a named invariant polynomial in fundamental-weight variables."""
    return _convenience_context(cartan_type, parabolic).generator_polynomial(generator)


def exceptional_image(cartan_type, parabolic, generator):
    """Compute a named invariant's exact Schubert expansion.

    Example: ``exceptional_image('F4', parabolic=1, generator='e2')``.
    Use ``exceptional_borel`` to reuse a context and access its weights/ring.
    """
    return _convenience_context(cartan_type, parabolic).generator_image(generator)
