"""Integral image lattices of fixed rational Borel generators.

The Schubert classes form an integral basis of the cohomology of a complex
flag variety.  Once each selected generator has integral Schubert
coefficients, its products define a map from the weighted polynomial ring
over ZZ.  In degree d, its image is the column lattice of the integer
Schubert coefficient matrix M.  A Smith decomposition U*M*V=D records both
the image's saturation and the finite quotient of that saturation by the
image.  No denominators of the fixed input generators are cleared.
"""

from dataclasses import dataclass, field, replace

from sage.all import QQ, ZZ, identity_matrix, matrix, prod, vector

from .expansion import Expansion


def _matrix_dict(value):
    """Serialize a possibly empty integer matrix without losing its shape."""
    return {
        "rows": int(value.nrows()),
        "columns": int(value.ncols()),
        "entries": [
            {"row": int(row), "column": int(column), "coefficient": str(coefficient)}
            for (row, column), coefficient in sorted(value.dict().items())
            if coefficient
        ],
    }


@dataclass(frozen=True)
class IntegralDegreeAnalysis:
    r"""One exact integral image lattice in polynomial degree ``degree``.

    Columns of ``coefficient_matrix`` are the images of ``monomials``; rows
    are the stored ``schubert_indices``.  Only Schubert classes occurring in
    these columns are stored.  Both lattice-basis matrices also use those
    rows and have ``rank`` columns.  All matrices are immutable over ZZ.

    ``smith_factors`` includes every positive nonzero Smith diagonal entry,
    including 1.  The factors greater than 1 describe
    ``saturation(image)/image``.  Its order is ``saturation_index``; a
    saturated image can still have smaller rank than the full cohomology.
    Ambient dimensions are obtained only when explicitly requested.
    """

    degree: int
    domain_dimension: int
    image_dimension: int
    monomials: tuple
    schubert_indices: tuple
    coefficient_matrix: object
    smith_factors: tuple
    torsion_factors: tuple
    image_lattice_basis: object
    saturation_lattice_basis: object
    _smith_left: object = field(repr=False, compare=False)
    _smith_right: object = field(repr=False, compare=False)
    _borel_map: object = field(repr=False, compare=False)
    _domain: object = field(repr=False, compare=False)
    ambient_dimension: object = None

    @property
    def rank(self):
        return self.image_dimension

    @property
    def kernel_dimension(self):
        """Rank of the kernel of the integer coefficient map."""
        return self.domain_dimension - self.image_dimension

    @property
    def matrix(self):
        return self.coefficient_matrix

    @property
    def saturation_index(self):
        """Index in the saturation, with value 1 for the zero lattice."""
        return ZZ(prod(self.smith_factors))

    @property
    def is_saturated(self):
        return not self.torsion_factors

    @property
    def free_cokernel_rank(self):
        """Free rank of full cohomology/image, if its dimension was requested."""
        if self.ambient_dimension is None:
            return None
        return self.ambient_dimension - self.image_dimension

    @property
    def full_ambient(self):
        """Whether the image spans the full rational degree, if known.

        This does not assert integral surjectivity: nontrivial Smith factors
        can still obstruct it when the ranks agree.
        """
        if self.ambient_dimension is None:
            return None
        return self.image_dimension == self.ambient_dimension

    def lattice(self):
        """Return an embedded lattice for exact sums, intersections and quotients."""
        from .lattices import IntegralImageLattice
        # integral_degree already computed and froze these exact bases.
        return IntegralImageLattice(
            self.degree, self.schubert_indices, self.image_lattice_basis,
            self.smith_factors, self.saturation_lattice_basis, self._borel_map)

    def compare(self, other):
        """Compare with another degree report or an IntegralImageLattice."""
        from .lattices import compare_integral_images
        return compare_integral_images(self, other)

    def _integer_preimage_coefficients(self, expansion):
        if not isinstance(expansion, Expansion):
            raise TypeError("cohomology class must be a Schubert Expansion")
        self._borel_map._expansion()._check_context(expansion)
        if any(int(word.length()) != self.degree for word in expansion):
            raise ValueError("cohomology class must be homogeneous of polynomial degree {} "
                             "(or zero)".format(self.degree))

        positions = {word: row for row, word in enumerate(self.schubert_indices)}
        coefficients = [ZZ.zero()] * len(self.schubert_indices)
        for word, coefficient in expansion.items():
            if word not in positions or QQ(coefficient).denominator() != 1:
                return None
            coefficients[positions[word]] = ZZ(coefficient)
        transformed = self._smith_left * vector(ZZ, coefficients)
        if any(transformed[self.rank:]):
            return None
        solution = [ZZ.zero()] * self.domain_dimension
        for position, factor in enumerate(self.smith_factors):
            if transformed[position] % factor:
                return None
            solution[position] = transformed[position] // factor
        # U*M*V=D, so M*(V*y)=target whenever D*y=U*target.
        return tuple(self._smith_right * vector(ZZ, solution))

    def contains(self, expansion):
        """Test membership in the fixed generators' integer image lattice.

        Zero is allowed in every degree.  A wrong flag-variety context or a
        nonzero term in a different degree raises ValueError.  A rational
        class with noninteger coefficients, or a same-degree class outside
        the column lattice, returns False.
        """
        return self._integer_preimage_coefficients(expansion) is not None

    def preimage(self, expansion):
        """Return a formal polynomial with integer coefficients mapping to a class.

        The polynomial lies in the original algebra's rational formal ring,
        but all its coefficients are integers.  It has this report's weighted
        degree, unless it is zero.  Integer preimages are not unique when
        there are relations.  A class outside the lattice raises ValueError.
        """
        coefficients = self._integer_preimage_coefficients(expansion)
        if coefficients is None:
            raise ValueError("cohomology class is not in the integral image lattice "
                             "in polynomial degree {}".format(self.degree))
        return self._domain(sum(coefficient * monomial
                                for coefficient, monomial in zip(coefficients, self.monomials)))

    def as_dict(self):
        """Return a JSON-safe report with exact integer matrix coefficients."""
        return {
            "degree": self.degree,
            "cohomological_degree": 2 * self.degree,
            "domain_dimension": self.domain_dimension,
            "image_dimension": self.image_dimension,
            "rank": self.rank,
            "kernel_dimension": self.kernel_dimension,
            "ambient_dimension": self.ambient_dimension,
            "free_cokernel_rank": self.free_cokernel_rank,
            "full_ambient": self.full_ambient,
            "smith_factors": [int(value) for value in self.smith_factors],
            "torsion_factors": [int(value) for value in self.torsion_factors],
            "saturation_index": int(self.saturation_index),
            "is_saturated": self.is_saturated,
            "monomials": [str(monomial) for monomial in self.monomials],
            "schubert_indices": [
                [int(node) for node in element.reduced_word()]
                for element in self.schubert_indices
            ],
            "matrix": _matrix_dict(self.coefficient_matrix),
            "image_lattice_basis": _matrix_dict(self.image_lattice_basis),
            "saturation_lattice_basis": _matrix_dict(self.saturation_lattice_basis),
        }


class IntegralMixin:
    """Integral-lattice analysis for a ``GeneratorAlgebra`` with fixed inputs."""

    def _integral_generators(self):
        if getattr(self, "_integral_generator_images", None) is None:
            images = []
            for name, polynomial in zip(self.generator_names, self.generator_polynomials):
                expansion = self.borel_map.expand(polynomial)
                nonintegral = next(((word, coefficient)
                                    for word, coefficient in expansion.items()
                                    if QQ(coefficient).denominator() != 1), None)
                if nonintegral is not None:
                    word, coefficient = nonintegral
                    raise ValueError(
                        "generator {} has a nonintegral Schubert coefficient {} at word {}; "
                        "integral images require integral generator classes, with the fixed "
                        "input normalization (denominators are not cleared)".format(
                            name, coefficient, tuple(int(i) for i in word.reduced_word())))
                images.append(expansion)
            self._integral_generator_images = tuple(images)
        return self._integral_generator_images

    def integral_degree(self, degree, include_ambient=False):
        r"""Compute the integral image lattice and Smith factors in one degree.

        The input polynomials retain their original normalization.  They may
        have rational Cartan coefficients, but every selected generator must
        have an integral Schubert image; this is checked even when requesting
        degree zero or a degree with no formal monomials.  An invalid generator
        raises ValueError, rather than silently clearing its denominators.

        The rational degree report supplies the same formal monomial columns
        and Schubert-support rows.  Their integer matrix M has Smith form
        U*M*V=D.  The nonzero diagonal entries are the invariant factors of
        ``saturation(image)/image``; they do not imply that the rational span
        is all of H^(2d)(G/P).  Set ``include_ambient=True`` to obtain that
        ambient dimension and the free cokernel rank.  The default does not
        enumerate the ambient Schubert basis.
        """
        # Lazy import avoids a cycle: GeneratorAlgebra inherits this mixin.
        from .algebra import _degree

        d = _degree(degree)
        if not isinstance(include_ambient, bool):
            raise ValueError("include_ambient must be a boolean")
        self._integral_generators()
        if not hasattr(self, "_integral_degree_cache"):
            self._integral_degree_cache = {}
        key = (d, include_ambient)
        if key in self._integral_degree_cache:
            return self._integral_degree_cache[key]

        if include_ambient and (d, False) in self._integral_degree_cache:
            dimension = self.degree(d, include_ambient=True).ambient_dimension
            result = replace(self._integral_degree_cache[d, False],
                             ambient_dimension=dimension)
            self._integral_degree_cache[key] = result
            return result

        rational = self.degree(d, include_ambient=include_ambient)
        try:
            coefficients = matrix(ZZ, rational.coefficient_matrix)
        except (TypeError, ValueError) as error:
            raise ArithmeticError("products of integral generator classes unexpectedly have "
                                  "nonintegral Schubert coefficients") from error
        rows, columns = coefficients.nrows(), coefficients.ncols()
        if rows and columns:
            diagonal, left, right = coefficients.smith_form()
        else:
            # Preserve empty matrix shapes; the two changes of basis remain
            # square identities even when there are no monomials or no support.
            diagonal = matrix(ZZ, rows, columns)
            left = identity_matrix(ZZ, rows)
            right = identity_matrix(ZZ, columns)
        if left * coefficients * right != diagonal:
            raise ArithmeticError("Smith decomposition failed its exact transformation check")
        factors = tuple(abs(ZZ(diagonal[i, i]))
                        for i in range(min(rows, columns)) if diagonal[i, i])
        rank = len(factors)
        if rank != rational.image_dimension:
            raise ArithmeticError("integer and rational image ranks disagree")

        # Normalize diagonal signs if the backend returns a negative entry.
        # Multiplying row i of U by -1 preserves unimodularity and U*M*V=D.
        for i in range(rank):
            if diagonal[i, i] < 0:
                left.rescale_row(i, -1)
                diagonal.rescale_row(i, -1)
        inverse_left = matrix(ZZ, left.inverse())
        saturation = inverse_left.matrix_from_columns(range(rank))
        image = saturation * matrix(ZZ, rank, rank,
                                    {(i, i): factor for i, factor in enumerate(factors)})
        for value in (coefficients, left, right, image, saturation):
            value.set_immutable()
        result = IntegralDegreeAnalysis(
            degree=d,
            domain_dimension=rational.domain_dimension,
            image_dimension=rank,
            monomials=rational.monomials,
            schubert_indices=rational.schubert_indices,
            coefficient_matrix=coefficients,
            smith_factors=factors,
            torsion_factors=tuple(factor for factor in factors if factor > 1),
            image_lattice_basis=image,
            saturation_lattice_basis=saturation,
            _smith_left=left,
            _smith_right=right,
            _borel_map=self.borel_map,
            _domain=self.domain,
            ambient_dimension=rational.ambient_dimension,
        )
        self._integral_degree_cache[key] = result
        if include_ambient and (d, False) not in self._integral_degree_cache:
            self._integral_degree_cache[d, False] = replace(result, ambient_dimension=None)
        return result
