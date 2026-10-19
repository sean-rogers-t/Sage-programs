"""Exact operations on degreewise integral Schubert image lattices.

All bases have columns in a shared integral Schubert coordinate system.
Smith transformations compute *integer* kernels and bases: clearing the
denominators of an arbitrary rational kernel basis would not suffice.
"""

from dataclasses import dataclass, field

from sage.all import QQ, ZZ, identity_matrix, matrix, prod, vector

from .expansion import Expansion
from .integral import _matrix_dict


def _smith(value):
    if value.nrows() and value.ncols():
        diagonal, left, right = value.smith_form()
    else:
        diagonal = matrix(ZZ, value.nrows(), value.ncols())
        left = identity_matrix(ZZ, value.nrows())
        right = identity_matrix(ZZ, value.ncols())
    rank = sum(bool(diagonal[i, i]) for i in range(min(value.dimensions())))
    return diagonal, left, right, rank


def _immutable(value):
    value.set_immutable()
    return value


@dataclass(frozen=True)
class IntegralLatticeQuotient:
    """The abstract finitely generated group ``larger / smaller``.

    ``index`` is None when there is a nonzero free part; otherwise it is
    the finite index, including 1 for the zero group. Unit factors are omitted.
    """

    free_rank: int
    torsion_factors: tuple

    @property
    def index(self):
        return None if self.free_rank else ZZ(prod(self.torsion_factors))

    def as_dict(self):
        return {"free_rank": int(self.free_rank),
                "torsion_factors": [int(x) for x in self.torsion_factors],
                "index": None if self.index is None else int(self.index)}


@dataclass(frozen=True, eq=False)
class IntegralImageLattice:
    """An embedded integral lattice, with no preferred source polynomial lift.

    Obtain this from ``IntegralDegreeAnalysis.lattice()`` or a comparison's
    ``sum`` / ``intersection``. Matrices are immutable, with basis columns
    indexed by ``schubert_indices``. Degree d means cohomology H^(2d).
    """

    degree: int
    schubert_indices: tuple
    basis: object
    smith_factors: tuple
    saturation_basis: object
    _borel_map: object = field(repr=False)

    @classmethod
    def _from_columns(cls, degree, indices, columns, borel_map):
        columns = matrix(ZZ, columns)
        if columns.nrows() != len(indices):
            raise ValueError("Schubert indices and matrix rows do not agree")
        # Remove unused coordinates, also giving the zero lattice a unique
        # empty support. Retain deterministic ordering supplied by the caller.
        used = [i for i in range(columns.nrows()) if any(columns.row(i))]
        indices = tuple(indices[i] for i in used)
        columns = columns.matrix_from_rows(used)
        diagonal, left, _, rank = _smith(columns)
        saturation = matrix(ZZ, left.inverse()).matrix_from_columns(range(rank))
        factors = tuple(abs(ZZ(diagonal[i, i])) for i in range(rank))
        basis = saturation * matrix.diagonal(ZZ, factors)
        return cls(degree, indices, _immutable(basis), factors,
                   _immutable(saturation), borel_map)

    @property
    def rank(self):
        return self.basis.ncols()

    @property
    def saturation_index(self):
        return ZZ(prod(self.smith_factors))

    @property
    def is_saturated(self):
        return self.saturation_index == 1

    def contains(self, expansion):
        if not isinstance(expansion, Expansion):
            raise TypeError("expected a Schubert Expansion")
        self._borel_map._expansion()._check_context(expansion)
        if any(w.length() != self.degree for w in expansion):
            raise ValueError("expansion has the wrong degree")
        if any(w not in self.schubert_indices or c not in ZZ
               for w, c in expansion.items()):
            return False
        target = vector(QQ, [expansion[w] for w in self.schubert_indices])
        try:
            coordinates = self.basis.change_ring(QQ).solve_right(target)
        except ValueError:
            return False
        return all(c in ZZ for c in coordinates)

    def compare(self, other):
        return compare_integral_images(self, other)

    def quotient_by(self, smaller):
        """Return self / smaller; raise ValueError unless smaller is contained."""
        smaller = _as_lattice(smaller)
        _, larger_basis, smaller_basis = _align(self, smaller)
        coordinates = _inclusion_coordinates(smaller_basis, larger_basis)
        if coordinates is None:
            raise ValueError("the proposed denominator is not a sublattice")
        diagonal, _, _, rank = _smith(coordinates)
        torsion = tuple(abs(ZZ(diagonal[i, i])) for i in range(rank)
                        if abs(diagonal[i, i]) > 1)
        return IntegralLatticeQuotient(self.rank - smaller.rank, torsion)

    def as_dict(self):
        return {"degree": int(self.degree), "rank": int(self.rank),
                "cartan_type": str(self._borel_map.cartan_type),
                "omitted_nodes": [int(i) for i in self._borel_map.omitted_nodes],
                "schubert_indices": [[int(i) for i in w.reduced_word()]
                                     for w in self.schubert_indices],
                "basis": _matrix_dict(self.basis),
                "smith_factors": [int(x) for x in self.smith_factors],
                "saturation_index": int(self.saturation_index),
                "saturation_basis": _matrix_dict(self.saturation_basis)}


def _as_lattice(value):
    from .integral import IntegralDegreeAnalysis
    if isinstance(value, IntegralDegreeAnalysis):
        return value.lattice()
    if not isinstance(value, IntegralImageLattice):
        raise TypeError("expected IntegralDegreeAnalysis or IntegralImageLattice")
    return value


def _align(left, right):
    left._borel_map._expansion()._check_context(right._borel_map._expansion())
    if left.degree != right.degree:
        raise ValueError("integral image comparisons require the same degree")
    indices = tuple(sorted(set(left.schubert_indices) | set(right.schubert_indices),
                           key=lambda w: tuple(w.reduced_word())))
    positions = {w: i for i, w in enumerate(indices)}

    def embed(lattice):
        return matrix(ZZ, len(indices), lattice.rank,
                      {(positions[w], j): lattice.basis[i, j]
                       for i, w in enumerate(lattice.schubert_indices)
                       for j in range(lattice.rank)})

    return indices, embed(left), embed(right)


def _inclusion_coordinates(smaller, larger):
    # Both have independent columns. Solve over Q then require Z coordinates.
    try:
        coordinates = larger.change_ring(QQ).solve_right(smaller)
    except ValueError:
        return None
    if any(c not in ZZ for c in coordinates.list()):
        return None
    return matrix(ZZ, coordinates)


@dataclass(frozen=True)
class IntegralImageComparison:
    left: IntegralImageLattice
    right: IntegralImageLattice
    sum: IntegralImageLattice
    intersection: IntegralImageLattice
    left_in_right: bool
    right_in_left: bool

    @property
    def equal(self):
        return self.left_in_right and self.right_in_left

    @property
    def same_rational_span(self):
        return self.left.rank == self.right.rank == self.sum.rank

    def as_dict(self):
        return {"left": self.left.as_dict(), "right": self.right.as_dict(),
                "sum": self.sum.as_dict(), "intersection": self.intersection.as_dict(),
                "left_in_right": self.left_in_right, "right_in_left": self.right_in_left,
                "equal": self.equal, "same_rational_span": self.same_rational_span,
                "sum_mod_left": self.sum.quotient_by(self.left).as_dict(),
                "sum_mod_right": self.sum.quotient_by(self.right).as_dict(),
                "left_mod_intersection": self.left.quotient_by(self.intersection).as_dict(),
                "right_mod_intersection": self.right.quotient_by(self.intersection).as_dict()}


def compare_integral_images(left, right):
    """Compare two embedded image lattices on the same flag in the same degree.

    Accept integral_degree reports or lattices returned by earlier operations.
    This computes degreewise subgroup sums, not the algebra generated by both
    sources: mixed products must be supplied through a joint GeneratorAlgebra.
    """
    left, right = _as_lattice(left), _as_lattice(right)
    indices, a, b = _align(left, right)
    total = IntegralImageLattice._from_columns(
        left.degree, indices, a.augment(b), left._borel_map)
    # For U [A,-B] V = D, the final columns of V are a Z-basis of
    # the integer kernel. Project to A-coordinates, then multiply by A.
    joined = a.augment(-b)
    _, _, transform, rank = _smith(joined)
    kernel = transform.matrix_from_columns(range(rank, joined.ncols()))
    common = a * kernel.matrix_from_rows(range(left.rank))
    intersection = IntegralImageLattice._from_columns(
        left.degree, indices, common, left._borel_map)
    return IntegralImageComparison(
        left, right, total, intersection,
        _inclusion_coordinates(a, b) is not None,
        _inclusion_coordinates(b, a) is not None)
