"""Closed Smith multiplicities for the fixed orthogonal Springer sources.

The formulas are the unpublished, internally audited deductions in
TYPE_B_GENERAL_SMITH.md and TYPE_D_EVEN_GLUING_THEOREM.md.  This module
evaluates them with exact integer arithmetic; it constructs neither a
Weyl group nor a Schubert matrix.  It does not compute Smith basis vectors.
"""

from dataclasses import dataclass

from .operations import _integer


def _gaussian_prefix(n, k, limit, cache):
    """Coefficients of [n choose k]_q through limit, via its product formula."""
    k = min(k, n - k)
    limit = min(limit, k * (n - k))
    key = (n, k, limit)
    if key not in cache:
        coefficients = [1] + [0] * limit
        # Factors with i > limit are 1 modulo q^(limit+1).  The order
        # below divides by 1-q^i, then multiplies by 1-q^(n-k+i).
        for i in range(1, min(k, limit) + 1):
            for d in range(i, limit + 1):
                coefficients[d] += coefficients[d - i]
            shift = n - k + i
            for d in range(limit, shift - 1, -1):
                coefficients[d] -= coefficients[d - shift]
        cache[key] = tuple(coefficients)
    return cache[key]


def _gaussian_product_coefficient(parameters, degree, cache):
    """One coefficient of a product, truncating and using palindromicity."""
    if any(n < 0 or k < 0 or k > n for n, k in parameters):
        return 0
    top = sum(k * (n - k) for n, k in parameters)
    if not 0 <= degree <= top:
        return 0
    degree = min(degree, top - degree)
    if degree == 0:
        return 1
    product = (1,)
    for n, k in parameters:
        right = _gaussian_prefix(n, k, degree, cache)
        if len(right) == 1:
            continue
        result = [0] * min(degree + 1, len(product) + len(right) - 1)
        for i, first in enumerate(product):
            if first:
                for j in range(min(len(right), len(result) - i)):
                    result[i + j] += first * right[j]
        product = tuple(result)
    return product[degree] if degree < len(product) else 0


def _type_b_counts(n, k, degree, cache):
    """Exponent counts, including the point n=k used by the D Euler term."""
    if degree < 0 or degree % 2 or degree > 2 * k * (n - k):
        return {}
    if n == k:
        return {0: 1} if degree == 0 else {}
    length, m = min(k, n - k), degree // 2
    r, c = n // 2, (n + 1) // 2
    counts = {}
    # Expand q^(e*c) [r choose e]_q U_(n-2e,length-e)(q).
    for exponent in range(min(length, m // c) + 1):
        count = 0
        for j in range((length - exponent) // 2 + 1):
            remaining = m - exponent * c - j * (c - length + j)
            if remaining < 0:
                break
            count += _gaussian_product_coefficient(
                ((r, exponent), (r - exponent, j),
                 (c - exponent - j, length - exponent - 2 * j)),
                remaining, cache)
        if count:
            counts[exponent] = count
    return counts


def _unit_gluing(n, k, degree, cache):
    if k % 2 or degree % 2:
        return 0
    a, h = n - k, k // 2
    index = degree // 2 - h * (a // 2 + 1)
    return _gaussian_product_coefficient(
        ((h + (a + 1) // 2 - 1, h),), index, cache)


@dataclass(frozen=True)
class OrthogonalSmithSummary:
    """Compressed Smith factors in polynomial degree ``degree``.

    ``exponent_multiplicities`` is a tuple of (e, count) pairs for the
    factors 2^e, including units.  ``multiplicities`` returns a fresh
    factor-to-count dictionary.  The index is in the saturation of the
    selected rational image, not necessarily in all cohomology.

    Use ``smith_factors()`` or ``saturation_index()`` to explicitly
    materialize potentially large outputs.  Neither is expanded by
    ``as_dict()``.  Membership, preimages, and lattice basis vectors
    require ``GeneratorAlgebra.integral_degree()`` instead.
    """

    family: str
    n: int
    k: int
    degree: int
    sector: str
    normalization: str
    exponent_multiplicities: tuple
    unit_gluing_dimension: int

    @property
    def cohomological_degree(self):
        return 2 * self.degree

    @property
    def rank(self):
        return sum(count for _, count in self.exponent_multiplicities)

    @property
    def multiplicities(self):
        return {2 ** exponent: count
                for exponent, count in self.exponent_multiplicities}

    @property
    def saturation_exponent(self):
        """The exact integer E such that [saturation : image] = 2^E."""
        return sum(exponent * count
                   for exponent, count in self.exponent_multiplicities)

    @property
    def is_saturated(self):
        return self.saturation_exponent == 0

    def smith_factors(self, max_rank=100000):
        """Expand the sorted tuple, rejecting rank above the explicit limit."""
        max_rank = _integer(max_rank, "max_rank")
        if self.rank > max_rank:
            raise ValueError("rank exceeds max_rank; use multiplicities or "
                             "increase max_rank explicitly")
        return tuple(2 ** exponent
                     for exponent, count in self.exponent_multiplicities
                     for _ in range(count))

    def saturation_index(self, max_bits=1000000):
        """Materialize 2^E, subject to its bit-length limit; zero image has 1."""
        max_bits = _integer(max_bits, "max_bits", minimum=1)
        exponent = self.saturation_exponent
        if exponent + 1 > max_bits:
            raise ValueError("index exceeds max_bits; use saturation_exponent "
                             "or increase max_bits explicitly")
        return 1 << exponent

    def as_dict(self):
        """JSON-safe compressed report; no expanded factors or enormous index."""
        return {
            "family": self.family, "n": self.n, "k": self.k,
            "degree": self.degree,
            "cohomological_degree": self.cohomological_degree,
            "sector": self.sector, "normalization": self.normalization,
            "rank": self.rank,
            "exponent_multiplicities": [
                {"exponent": exponent, "multiplicity": count}
                for exponent, count in self.exponent_multiplicities],
            "saturation_exponent": self.saturation_exponent,
            "is_saturated": self.is_saturated,
            "normalized_unit_gluing_dimension": self.unit_gluing_dimension,
            "scope": "image in its rational saturation; no ambient basis computed",
            "method": "closed orthogonal Smith formulas",
        }


def orthogonal_smith(family, n, k, degree, *, sector="full",
                     normalization="normalized"):
    r"""Compute the closed B/D Springer-image Smith formula without matrices.

    Parameters use the Kumar papers: ``OG(n-k, 2*n+1)`` in B and
    ``OG(n-k, 2*n)`` in D.  Ranges: B n>=2, 1<=k<n;
    D n>=4, 2<=k<n.  ``degree=d`` means H^(2d), as in integral_degree.

    The normalized inputs are q_i=e_i(x_tail^2), and in D also
    p=prod(x_tail), with p^2=q_k.  ``sector`` is ``full``, ``square``,
    or (D only) ``euler``; the last means the module p*ZZ[q], not the
    algebra generated by p alone.  ``character`` normalization uses
    4^i*q_i and 2^k*p, multiplying every degree-d factor by 2^d.
    This does not change the character lattice to Spin.

    ``unit_gluing_dimension`` counts normalized unit gluing for the
    full D image; it is zero for sector-only and B reports.  It is
    unchanged by normalization and is not a count of unit factors
    in the scaled character form.

    The formula gives multiplicities, not explicit Smith directions.
    It is the internally audited unpublished result documented in
    TYPE_B_GENERAL_SMITH.md and TYPE_D_EVEN_GLUING_THEOREM.md.
    """
    if not isinstance(family, str) or family.upper() not in ("B", "D"):
        raise ValueError("family must be B or D")
    family = family.upper()
    n = _integer(n, "n", minimum=4 if family == "D" else 2)
    k = _integer(k, "k", minimum=2 if family == "D" else 1)
    degree = _integer(degree, "degree")
    if k >= n:
        raise ValueError("k must be smaller than n")
    if sector not in ("full", "square", "euler"):
        raise ValueError("sector must be full, square or euler")
    if family == "B" and sector == "euler":
        raise ValueError("the Euler sector is only defined in type D")
    if normalization not in ("normalized", "character"):
        raise ValueError("normalization must be normalized or character")

    # Cache only within this call, so a sequence of high-rank requests
    # cannot leave an unbounded global polynomial cache behind.
    cache, counts, gluing = {}, {}, 0
    if family == "B" or sector in ("full", "square"):
        counts = _type_b_counts(n, k, degree, cache)
    if family == "D" and sector in ("full", "euler"):
        for exponent, count in _type_b_counts(n - 1, k, degree - k, cache).items():
            counts[exponent] = counts.get(exponent, 0) + count
        if sector == "full":
            gluing = _unit_gluing(n, k, degree, cache)
            if gluing:
                counts[0] = counts.get(0, 0) - gluing
                counts[1] = counts.get(1, 0) + gluing
    if any(count < 0 for count in counts.values()):
        raise ArithmeticError("negative multiplicity in the orthogonal formula")
    shift = degree if normalization == "character" else 0
    pairs = tuple((exponent + shift, count)
                  for exponent, count in sorted(counts.items()) if count)
    return OrthogonalSmithSummary(family, n, k, degree, sector,
                                  normalization, pairs, gluing)
