"""Exact images and relations of chosen Borel generators.

The degreewise reports describe vector spaces in specified degrees. The separate
``presentation()`` method certifies the complete kernel ideal and image basis.
All degrees here are polynomial degrees, half the cohomological degrees.
"""

from collections.abc import Mapping
from dataclasses import dataclass, replace
from keyword import iskeyword
from operator import index as integer_index

from sage.all import PolynomialRing, QQ, matrix

from .core import BorelMap


def _degree(value):
    """Require an actual nonnegative integer, without rounding or truncation."""
    if isinstance(value, bool):
        raise ValueError("degree must be a nonnegative integer, not a boolean")
    try:
        result = integer_index(value)
    except (TypeError, ValueError) as error:
        raise ValueError("degree must be a nonnegative integer") from error
    if result < 0:
        raise ValueError("degree must be a nonnegative integer")
    return int(result)


def _polynomial_terms(polynomial):
    return [
        {"exponents": [int(e) for e in exponents], "coefficient": str(coefficient)}
        for exponents, coefficient in sorted(polynomial.dict().items())
    ]


@dataclass(frozen=True)
class DegreeAnalysis:
    """One exact homogeneous component of a chosen generator-algebra map.

    Columns of ``coefficient_matrix`` are ``monomials``; its rows are
    ``schubert_indices``.  Only Schubert classes that occur in those columns are
    stored.  ``relations`` is a vector-space basis of this degree of the kernel,
    rather than a claimed list of ideal generators.  The matrix is immutable.
    """

    degree: int
    domain_dimension: int
    image_dimension: int
    kernel_dimension: int
    monomials: tuple
    schubert_indices: tuple
    coefficient_matrix: object
    relations: tuple
    pivot_columns: tuple
    image_basis_monomials: tuple
    ambient_dimension: object = None

    @property
    def rank(self):
        return self.image_dimension

    @property
    def matrix(self):
        return self.coefficient_matrix

    def as_dict(self):
        """Return a JSON-safe report with exact rational coefficients."""
        return {
            "degree": self.degree,
            "cohomological_degree": 2 * self.degree,
            "domain_dimension": self.domain_dimension,
            "image_dimension": self.image_dimension,
            "rank": self.image_dimension,
            "kernel_dimension": self.kernel_dimension,
            "ambient_dimension": self.ambient_dimension,
            "monomials": [str(monomial) for monomial in self.monomials],
            "schubert_indices": [
                [int(node) for node in element.reduced_word()]
                for element in self.schubert_indices
            ],
            "relations": [
                {"polynomial": str(relation), "terms": _polynomial_terms(relation)}
                for relation in self.relations
            ],
            "pivot_columns": list(self.pivot_columns),
            "image_basis_monomials": [str(monomial) for monomial in self.image_basis_monomials],
            "matrix": {
                "rows": int(self.coefficient_matrix.nrows()),
                "columns": int(self.coefficient_matrix.ncols()),
                "entries": [
                    {"row": int(row), "column": int(column), "coefficient": str(value)}
                    for (row, column), value in sorted(self.coefficient_matrix.dict().items())
                    if value
                ],
            },
        }


class GeneratorAlgebra:
    """Map a weighted polynomial algebra to the cohomology of a flag variety.

    ``generators`` is a nonempty ordered mapping from distinct variable names
    to nonzero, homogeneous, positive-degree polynomials in ``borel_map.ring``.
    The polynomials must be invariant under its retained Levi nodes.  Algebraic
    independence is not assumed, and a generator is allowed to have zero Borel
    image.  Insertion order determines the formal variables and matrix columns.

    ``degree(d)`` computes the image and kernel exactly over QQ without asking
    for the full ambient Schubert basis.  Set ``include_ambient=True`` only when
    the dimension of that basis is also wanted.
    """

    def __init__(self, borel_map, generators):
        if not isinstance(borel_map, BorelMap):
            raise TypeError("borel_map must be a BorelMap")
        if not isinstance(generators, Mapping):
            raise TypeError("generators must be an ordered mapping of names to polynomials")
        entries = tuple(generators.items())
        if not entries:
            raise ValueError("at least one generator is required")
        names = tuple(name for name, _ in entries)
        if any(not isinstance(name, str) or not name.isidentifier() or iskeyword(name)
               for name in names):
            raise ValueError("generator names must be valid, non-keyword identifiers")
        if len(set(names)) != len(names):
            raise ValueError("generator names must be distinct")
        polynomials = tuple(borel_map._polynomial(polynomial) for _, polynomial in entries)
        degrees = []
        for name, polynomial in zip(names, polynomials):
            if not polynomial or not polynomial.is_homogeneous() or polynomial.total_degree() <= 0:
                raise ValueError("generator {} must be nonzero, homogeneous, and of positive degree".format(name))
            borel_map._validate_invariance(polynomial)
            degrees.append(int(polynomial.total_degree()))
        self._borel_map = borel_map
        self._generator_names = names
        self._generator_degrees = tuple(degrees)
        self._generator_polynomials = polynomials
        # The explicit variable count keeps the same multivariate interface
        # when only one generator was supplied.
        self._domain = PolynomialRing(QQ, len(names), names=names)
        self.clear_cache()

    @property
    def borel_map(self):
        return self._borel_map

    @property
    def generator_names(self):
        return self._generator_names

    @property
    def generator_degrees(self):
        return self._generator_degrees

    @property
    def generator_polynomials(self):
        return self._generator_polynomials

    @property
    def domain(self):
        return self._domain

    def clear_cache(self):
        """Clear this algebra's reports and monomial caches, retaining Borel data."""
        self._monomial_cache = {}
        self._substitution_cache = {}
        self._image_cache = {}
        self._degree_cache = {}
        self._ambient_degree_cache = {}
        self._presentation_cache = None
        self._generator_power_cache = [
            {0: self.borel_map.ring.one(), 1: polynomial}
            for polynomial in self.generator_polynomials
        ]

    def _formal_polynomial(self, polynomial):
        try:
            return self.domain(polynomial)
        except (TypeError, ValueError) as error:
            raise TypeError("input must be a polynomial over QQ in {}".format(self.generator_names)) from error

    def _substitute_monomial(self, exponents):
        if exponents not in self._substitution_cache:
            value = self.borel_map.ring.one()
            for position, exponent in enumerate(exponents):
                powers = self._generator_power_cache[position]
                if exponent not in powers:
                    powers[exponent] = self.generator_polynomials[position] ** exponent
                value *= powers[exponent]
            self._substitution_cache[exponents] = value
        return self._substitution_cache[exponents]

    def substitute(self, polynomial):
        """Substitute the chosen Cartan polynomials into a formal polynomial."""
        formal = self._formal_polynomial(polynomial)
        answer = self.borel_map.ring.zero()
        for exponents, coefficient in formal.dict().items():
            answer += coefficient * self._substitute_monomial(tuple(int(e) for e in exponents))
        return answer

    def _evaluate_monomial(self, exponents):
        if exponents not in self._image_cache:
            self._image_cache[exponents] = self.borel_map.expand(self._substitute_monomial(exponents))
        return self._image_cache[exponents]

    def evaluate(self, polynomial):
        """Return the exact Schubert expansion of a formal polynomial's image."""
        formal = self._formal_polynomial(polynomial)
        coefficients = {}
        for exponents, coefficient in formal.dict().items():
            expansion = self._evaluate_monomial(tuple(int(e) for e in exponents))
            for element, value in expansion.items():
                coefficients[element] = coefficients.get(element, QQ.zero()) + coefficient * value
        return self.borel_map._expansion(coefficients)

    def monomials(self, degree):
        """Return formal monomials of exactly the requested weighted degree.

        Enumeration follows only solutions of sum(e_i * generator_degree_i) =
        degree; it does not construct and discard all lower total-degree terms.
        Exponent tuples are ordered lexicographically.
        """
        d = _degree(degree)
        if d not in self._monomial_cache:
            exponents = []
            weights = self.generator_degrees

            def visit(position, remaining, prefix):
                weight = weights[position]
                if position == len(weights) - 1:
                    if remaining % weight == 0:
                        exponents.append(prefix + (remaining // weight,))
                    return
                for exponent in range(remaining // weight + 1):
                    visit(position + 1, remaining - exponent * weight, prefix + (exponent,))

            visit(0, d, ())
            self._monomial_cache[d] = tuple(
                self.domain({powers: QQ.one()}) for powers in exponents
            )
        return self._monomial_cache[d]

    def degree(self, degree, include_ambient=False):
        """Compute a homogeneous image matrix and a basis of its right kernel.

        The default stores only nonzero Schubert support.  Ambient dimensions
        are optional because obtaining a complete ambient basis can cost much
        more than computing the selected generator images.
        """
        d = _degree(degree)
        if not isinstance(include_ambient, bool):
            raise ValueError("include_ambient must be a boolean")
        if d not in self._degree_cache:
            monomials = self.monomials(d)
            expansions = tuple(
                self._evaluate_monomial(tuple(int(e) for e in monomial.exponents()[0]))
                for monomial in monomials
            )
            support = set()
            for expansion in expansions:
                support.update(expansion)
            indices = tuple(sorted(support, key=self.borel_map._word))
            positions = {element: row for row, element in enumerate(indices)}
            entries = {
                (positions[element], column): value
                for column, expansion in enumerate(expansions)
                for element, value in expansion.items()
            }
            coefficients = matrix(QQ, len(indices), len(monomials), entries, sparse=True)
            if not monomials:
                pivots = ()
                relations = ()
            elif not indices:
                pivots = ()
                relations = monomials
            else:
                pivots = tuple(int(column) for column in coefficients.pivots())
                relations = tuple(
                    self.domain(sum(coefficient * monomial
                                    for coefficient, monomial in zip(vector, monomials)))
                    for vector in coefficients.right_kernel().basis()
                )
            coefficients.set_immutable()
            rank = len(pivots)
            self._degree_cache[d] = DegreeAnalysis(
                degree=d,
                domain_dimension=len(monomials),
                image_dimension=rank,
                kernel_dimension=len(monomials) - rank,
                monomials=monomials,
                schubert_indices=indices,
                coefficient_matrix=coefficients,
                relations=relations,
                pivot_columns=pivots,
                image_basis_monomials=tuple(monomials[column] for column in pivots),
            )
        if include_ambient:
            if d not in self._ambient_degree_cache:
                self._ambient_degree_cache[d] = replace(
                    self._degree_cache[d], ambient_dimension=len(self.borel_map.basis(d))
                )
            return self._ambient_degree_cache[d]
        return self._degree_cache[d]

    def through_degree(self, maximum_degree, include_ambient=False):
        """Return immutable reports for degrees zero through a finite bound."""
        bound = _degree(maximum_degree)
        return tuple(self.degree(d, include_ambient=include_ambient) for d in range(bound + 1))

    def hilbert_function(self, maximum_degree):
        """Return image dimensions in degrees zero through a finite bound."""
        return tuple(report.image_dimension for report in self.through_degree(maximum_degree))

    def first_relation_degree(self, maximum_degree):
        """Find the first positive kernel degree up to a bound, or return None.

        None means no relation was found through that finite bound; it is not
        a claim that the complete map is injective.
        """
        bound = _degree(maximum_degree)
        for d in range(1, bound + 1):
            if self.degree(d).kernel_dimension:
                return d
        return None

    def presentation(self):
        """Return a certified full kernel ideal, image basis and Hilbert series.

        This finite staircase construction takes no guessed degree cutoff.
        It uses exact Schubert images and checks border closure and basis
        independence. Runtime depends on the selected image's dimension;
        a full large flag ring can still be expensive. The immutable result
        is cached until ``clear_cache()``.
        """
        if self._presentation_cache is None:
            from .presentation import build_presentation
            self._presentation_cache = build_presentation(self)
        return self._presentation_cache
