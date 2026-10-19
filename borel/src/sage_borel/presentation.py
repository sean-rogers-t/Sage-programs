"""Complete rational image presentations by exact monomial staircases.

This is a finite-algebra construction, not an inference from a Hilbert prefix.
All degrees are polynomial degrees (half the cohomological degrees).
"""

from dataclasses import dataclass, field, replace
from heapq import heappop, heappush

from sage.all import PolynomialRing, QQ

from .algebra import _polynomial_terms


def _exponents(monomial):
    return tuple(int(e) for e in monomial.exponents()[0])


def _divides(left, right):
    return all(a <= b for a, b in zip(left, right))


def _weight(exponents, degrees):
    return sum(e * d for e, d in zip(exponents, degrees))


def _order_key(exponents, degrees):
    # Refine positive weighted degree by the domain's ordinary degrevlex.
    # Every relation is weighted homogeneous, so its leading monomial agrees
    # with ordinary degrevlex even when the generator degrees differ.
    return (_weight(exponents, degrees), sum(exponents),
            tuple(-e for e in reversed(exponents)))


class _ImageEchelon:
    """Incremental sparse elimination, optionally carrying formal lifts."""

    def __init__(self, key):
        self.rows = {}
        self._key = key

    def insert(self, expansion, lift=None):
        vector = dict(expansion)
        while vector:
            pivot = min(vector, key=self._key)
            coefficient = vector[pivot]
            if pivot not in self.rows:
                vector = {index: value / coefficient for index, value in vector.items()}
                if lift is not None:
                    lift /= coefficient
                self.rows[pivot] = (vector, lift)
                return True, lift
            row, row_lift = self.rows[pivot]
            for index, value in row.items():
                residual = vector.get(index, QQ.zero()) - coefficient * value
                if residual:
                    vector[index] = residual
                else:
                    vector.pop(index, None)
            if lift is not None:
                lift -= coefficient * row_lift
        return False, lift


class _BggEvaluator:
    """The original exact Schubert-coordinate evaluator."""

    backend = "bgg"

    def __init__(self, algebra, top_degree):
        self.algebra = algebra
        self.top_degree = top_degree
        self.key = algebra.borel_map._word

    def image(self, exponents):
        if _weight(exponents, self.algebra.generator_degrees) > self.top_degree:
            return {}
        return self.algebra._evaluate_monomial(exponents)


def resolve_backend(algebra, backend):
    """Validate a backend name and resolve auto without constructing a quotient."""
    if not isinstance(backend, str) or backend not in ("auto", "coinvariant", "bgg"):
        raise ValueError("backend must be 'auto', 'coinvariant', or 'bgg'")
    if backend == "bgg":
        return backend
    from .coinvariant import CoinvariantEvaluator
    if CoinvariantEvaluator.supports(algebra):
        return "coinvariant"
    if backend == "auto":
        return "bgg"
    raise NotImplementedError(
        "the coinvariant backend requires standard-numbered types A, B, C or D; "
        "use backend='bgg' for this Cartan type")


def _evaluator(algebra, top_degree, backend):
    resolved = resolve_backend(algebra, backend)
    if resolved == "coinvariant":
        from .coinvariant import CoinvariantEvaluator
        return CoinvariantEvaluator(algebra, top_degree)
    return _BggEvaluator(algebra, top_degree)


def _polynomial_image(evaluator, polynomial):
    result = {}
    for exponents, coefficient in polynomial.dict().items():
        for index, value in evaluator.image(tuple(exponents)).items():
            residual = result.get(index, QQ.zero()) + coefficient * value
            if residual:
                result[index] = residual
            else:
                result.pop(index, None)
    return result


@dataclass(frozen=True)
class KernelCertificate:
    """Exact checks proving that a presentation describes the entire image.

    Border closure proves spanning; independent images prove injectivity of
    the proposed quotient. Together with vanishing relations they prove the
    full kernel equality. This certificate makes no novelty claim.
    """

    relation_images_zero: bool
    basis_images_independent: bool
    standard_monomials_reduced: bool
    border_closed: bool
    contains_unit: bool
    groebner_basis_reduced: bool
    summary_data_correct: bool
    ambient_top_degree: int
    basis_degree_ranks: tuple
    checked_relations: int
    checked_basis_monomials: int
    checked_border_monomials: int
    backend: str

    @property
    def certified(self):
        return all((self.relation_images_zero, self.basis_images_independent,
                    self.standard_monomials_reduced, self.border_closed,
                    self.contains_unit, self.groebner_basis_reduced,
                    self.summary_data_correct))

    def as_dict(self):
        return {
            "certified": self.certified,
            "backend": self.backend,
            "method": ("exact coinvariant normal forms and finite monomial border closure"
                       if self.backend == "coinvariant" else
                       "exact Schubert images and finite monomial border closure"),
            "checks": {
                name: getattr(self, name) for name in (
                    "relation_images_zero", "basis_images_independent",
                    "standard_monomials_reduced", "border_closed",
                    "contains_unit", "groebner_basis_reduced", "summary_data_correct")
            },
            "ambient_top_degree": self.ambient_top_degree,
            "basis_degree_ranks": [
                {"degree": degree, "rank": rank}
                for degree, rank in self.basis_degree_ranks
            ],
            "checked_relations": self.checked_relations,
            "checked_basis_monomials": self.checked_basis_monomials,
            "checked_border_monomials": self.checked_border_monomials,
        }


@dataclass(frozen=True)
class KernelPresentation:
    """A complete kernel ideal and a monomial basis of the rational image.

    ``relations`` is the reduced monic Groebner basis in ``domain``'s ordinary
    degrevlex order. ``hilbert_series`` is a finite QQ[t] polynomial, using the
    actual positive generator degrees. ``kernel_ideal()`` returns a fresh Sage
    ideal so mutating it cannot alter this immutable report.
    """

    domain: object
    generator_degrees: tuple
    relations: tuple
    standard_monomials: tuple
    hilbert_coefficients: tuple
    hilbert_series: object
    certificate: KernelCertificate
    backend: str
    _algebra: object = field(repr=False, compare=False)

    @property
    def dimension(self):
        return len(self.standard_monomials)

    def kernel_ideal(self):
        return self.domain.ideal(list(self.relations))

    def reduce(self, polynomial):
        """Return the unique representative in the standard-monomial basis."""
        return self._algebra._formal_polynomial(polynomial).reduce(list(self.relations))

    def verify(self, backend=None):
        """Recheck completeness, optionally with a different exact backend.

        The default uses this presentation's backend. ``backend='bgg'`` gives
        an independent Schubert-coordinate cross-check of a coinvariant result,
        which may be expensive at higher rank.
        """
        # Recompute the geometric bound instead of trusting recorded metadata.
        # An artificially lowered bound could otherwise hide a nonzero image.
        top_degree = int(self._algebra.borel_map.longest_representative().length())
        evaluator = _evaluator(self._algebra, top_degree,
                               self.backend if backend is None else backend)
        certificate = _verify(self._algebra, self.relations, self.standard_monomials,
                              top_degree, evaluator)
        coefficients, series = _hilbert_data(self.standard_monomials,
                                             self._algebra.generator_degrees)
        correct = (self.domain is self._algebra.domain
                   and self.generator_degrees == self._algebra.generator_degrees
                   and self.hilbert_coefficients == coefficients
                   and self.hilbert_series == series)
        return replace(certificate, summary_data_correct=bool(correct))

    def as_dict(self):
        """Return source data, relations, basis and certificate for exact replay."""
        from . import __version__
        algebra = self._algebra
        return {
            "schema_version": 2,
            "library_version": __version__,
            "backend": self.backend,
            "cartan_type": str(algebra.borel_map.cartan_type),
            "omitted_nodes": list(algebra.borel_map.omitted_nodes),
            "generator_names": list(algebra.generator_names),
            "generator_degrees": list(self.generator_degrees),
            "cartan_variables": [str(x) for x in algebra.borel_map.weights],
            "generator_polynomials": [
                {"polynomial": str(f), "terms": _polynomial_terms(f)}
                for f in algebra.generator_polynomials
            ],
            "term_order": str(self.domain.term_order()),
            "relations": [
                {"polynomial": str(f), "terms": _polynomial_terms(f)}
                for f in self.relations
            ],
            "standard_monomials": [
                {"polynomial": str(m), "exponents": list(_exponents(m))}
                for m in self.standard_monomials
            ],
            "dimension": self.dimension,
            "hilbert_coefficients": list(self.hilbert_coefficients),
            "hilbert_series": str(self.hilbert_series),
            "certificate": self.certificate.as_dict(),
        }


def _verify(algebra, relations, monomials, top_degree, evaluator):
    degrees = algebra.generator_degrees
    exponents = tuple(_exponents(m) for m in monomials)
    standard = set(exponents)
    leading = tuple(_exponents(f.lm()) for f in relations)
    reduced = len(standard) == len(exponents) and all(
        not any(_divides(lm, alpha) for lm in leading) for alpha in exponents)
    unit = (0,) * len(degrees) in standard
    groebner_reduced = len(set(leading)) == len(leading)
    for position, relation in enumerate(relations):
        terms = relation.dict()
        lm = leading[position]
        groebner_reduced = groebner_reduced and relation.lc() == 1 and all(
            _weight(tuple(alpha), degrees) == _weight(lm, degrees)
            for alpha in terms)
        groebner_reduced = groebner_reduced and not any(
            _divides(other, lm) for j, other in enumerate(leading) if j != position)
        groebner_reduced = groebner_reduced and all(
            tuple(alpha) in standard for alpha in terms if tuple(alpha) != lm)

    border = set()
    for alpha in standard:
        for position in range(len(degrees)):
            child = list(alpha)
            child[position] += 1
            child = tuple(child)
            if child not in standard:
                border.add(child)
    closed = all(any(_divides(lm, alpha) for lm in leading) for alpha in border)
    zero_relations = all(not _polynomial_image(evaluator, f) for f in relations)
    echelon = {}
    independent = True
    for alpha in exponents:
        degree = _weight(alpha, degrees)
        if degree not in echelon:
            echelon[degree] = _ImageEchelon(evaluator.key)
        inserted, _ = echelon[degree].insert(evaluator.image(alpha))
        independent = independent and inserted
    return KernelCertificate(
        relation_images_zero=bool(zero_relations),
        basis_images_independent=bool(independent),
        standard_monomials_reduced=bool(reduced),
        border_closed=bool(closed),
        contains_unit=bool(unit),
        groebner_basis_reduced=bool(groebner_reduced),
        summary_data_correct=True,
        ambient_top_degree=top_degree,
        basis_degree_ranks=tuple((degree, len(rows.rows))
                                for degree, rows in sorted(echelon.items())),
        checked_relations=len(relations),
        checked_basis_monomials=len(monomials),
        checked_border_monomials=len(border),
        backend=evaluator.backend,
    )


def _hilbert_data(monomials, degrees):
    highest_degree = max((_weight(_exponents(m), degrees) for m in monomials),
                         default=0)
    coefficients = [0] * (highest_degree + 1)
    for monomial in monomials:
        coefficients[_weight(_exponents(monomial), degrees)] += 1
    return tuple(coefficients), PolynomialRing(QQ, "t")(coefficients)


def build_presentation(algebra, backend="auto"):
    """Construct and certify a complete presentation without a degree cutoff."""
    ring = algebra.domain
    degrees = algebra.generator_degrees
    count = len(degrees)
    top_degree = int(algebra.borel_map.longest_representative().length())
    evaluator = _evaluator(algebra, top_degree, backend)
    queue, queued = [], set()

    def enqueue(alpha):
        if alpha not in queued:
            queued.add(alpha)
            heappush(queue, (_order_key(alpha, degrees), alpha))

    enqueue((0,) * count)
    leading, relations, monomials = [], [], []
    echelon = {}
    while queue:
        _, alpha = heappop(queue)
        if any(_divides(lm, alpha) for lm in leading):
            continue
        degree = _weight(alpha, degrees)
        monomial = ring({alpha: QQ.one()})
        if degree not in echelon:
            echelon[degree] = _ImageEchelon(evaluator.key)
        independent, lift = echelon[degree].insert(
            evaluator.image(alpha), monomial)
        if independent:
            monomials.append(monomial)
            for position in range(count):
                child = list(alpha)
                child[position] += 1
                enqueue(tuple(child))
        else:
            relation = ring(lift)
            if relation.lc() != 1 or _exponents(relation.lm()) != alpha:
                raise ArithmeticError("staircase relation failed the leading-term check")
            relations.append(relation)
            leading.append(alpha)

    relations = tuple(sorted(relations, key=lambda f: f.lm(), reverse=True))
    monomials = tuple(monomials)
    certificate = _verify(algebra, relations, monomials, top_degree, evaluator)
    if not certificate.certified:
        raise ArithmeticError("full kernel presentation failed its completeness certificate")
    coefficients, series = _hilbert_data(monomials, degrees)
    return KernelPresentation(
        domain=ring, generator_degrees=degrees, relations=relations,
        standard_monomials=monomials, hilbert_coefficients=coefficients,
        hilbert_series=series, certificate=certificate,
        backend=evaluator.backend,
        _algebra=algebra,
    )
