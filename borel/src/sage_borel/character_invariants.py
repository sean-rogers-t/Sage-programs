"""Degreewise integral characters inside a specified rational invariant algebra."""
from dataclasses import dataclass, field
from math import comb

from sage.all import QQ, ZZ, diagonal_matrix, gcd, matrix, prod, vector

from .algebra import GeneratorAlgebra, _degree
from .characters import SpringerCharacterMap, _integer, springer_character_map
from .lattices import IntegralImageLattice


def _integral_preimage(M):
    """Exact preimage of Z^rows; M must have independent columns.

    Column-content preconditioning and small PARI HNF batches avoid the large
    intermediate integers of a single unconditioned Hermite calculation.
    The early exit tests every original constraint, not only the current batch.
    """
    if not M.ncols():
        return matrix(QQ, 0, 0)
    scales = diagonal_matrix(QQ, [1/gcd(col) for col in M.columns()])
    A = matrix(ZZ, M*scales)
    rows = sorted(A.rows(), key=lambda r: (max(abs(c) for c in r),
                                           sum(abs(c) for c in r)))
    H = matrix(ZZ, 0, A.ncols())
    for i in range(0, len(rows), 16):
        H = H.stack(matrix(ZZ, rows[i:i+16])).hermite_form(
            algorithm="pari4", include_zero_rows=False)
        if H.nrows() == A.ncols() and i % 128 == 0:
            B = scales*H.change_ring(QQ).inverse()
            if (M*B).denominator() == 1:
                return B
    if H.nrows() != M.ncols():
        raise ValueError("rational source monomials must be independent")
    B = scales*H.change_ring(QQ).inverse()
    assert (M*B).denominator() == 1
    return B


@dataclass(frozen=True)
class CharacterInvariantDegree:
    """A full Z-basis of the selected character source in one degree.

    Columns of basis_matrix give source_basis in the formal monomials.
    These are source bases, including elements with zero Borel image.
    They are not a claim that a fixed finite list generates the integral ring.
    """
    degree: int
    monomials: tuple
    basis_matrix: object
    source_basis: tuple
    _algebra: object = field(repr=False, compare=False)

    def image(self):
        """Return the embedded Schubert image, supporting lattice comparisons."""
        rational = self._algebra.degree(self.degree)
        assert rational.monomials == self.monomials
        columns = rational.matrix*self.basis_matrix
        assert columns.denominator() == 1
        return IntegralImageLattice._from_columns(
            self.degree, rational.schubert_indices, columns, self._algebra.borel_map)


class IntegralCharacterSource:
    """Intersect a rational polynomial invariant algebra with integral characters.

    The trace-source theorem certifies the torus source by default. An explicit
    basis_change requests the older leading-term certificate. The supplied
    GeneratorAlgebra must have algebraically independent Cartan polynomials.
    This computes the intersection in THAT rational algebra; fullness of the
    Levi invariant algebra is not inferred for arbitrary user inputs.
    """

    def __init__(self, springer_map, rational_algebra, basis_change=None):
        if not isinstance(springer_map, SpringerCharacterMap):
            raise TypeError("springer_map must be a SpringerCharacterMap")
        if not isinstance(rational_algebra, GeneratorAlgebra):
            raise TypeError("rational_algebra must be a GeneratorAlgebra")
        if springer_map.cartan_type != rational_algebra.borel_map.cartan_type:
            raise ValueError("Springer map and invariant algebra must have the same type")
        certificate = (springer_map.certify_trace_source() if basis_change is None
                       else springer_map.certify_polynomial_source(basis_change))
        if not certificate.certified:
            raise ValueError("an all-primes polynomial-source certificate is required")
        fs = rational_algebra.generator_polynomials
        variables = rational_algebra.borel_map.weights
        jacobian = matrix([[f.derivative(x) for x in variables] for f in fs])
        if jacobian.rank() != len(fs):
            raise ValueError("rational generators must be algebraically independent")
        self.springer_map = springer_map
        self.rational_algebra = rational_algebra
        self.torus_certificate = certificate
        # Work in Lambda/c to remove a harmless common homogeneous factor.
        self._scale = ZZ(gcd(springer_map.linear_basis_matrix.list()))
        B = springer_map.linear_basis_matrix/self._scale
        # Columns of B are linear forms: u = transpose(B)*w.
        substitutions = B.transpose().inverse()*vector(variables)
        self._inputs = tuple(f(*substitutions) for f in fs)
        self._coefficient_weights = (1,)*springer_map.rank
        self._cache = {}

    def degree(self, degree, max_monomials=500, max_matrix_entries=2000000):
        """Compute the exact homogeneous source lattice, retaining kernel inputs.

        Budgets bound the number of rational monomials and the full possible
        coefficient-matrix size before expansion. They are not time limits.
        """
        d = _degree(degree)
        bound = _integer(max_monomials, "max_monomials", 1)
        entries = _integer(max_matrix_entries, "max_matrix_entries", 1)
        algebra = self.rational_algebra
        monomials = algebra.monomials(d)
        if len(monomials) > bound:
            raise ValueError("monomial budget exceeded; raise max_monomials explicitly")
        weights = self._coefficient_weights
        if all(w == 1 for w in weights):
            possible_rows = comb(d+len(weights)-1, len(weights)-1)
        else:
            counts = [1]+[0]*d
            for w in weights:
                for k in range(w, d+1):
                    counts[k] += counts[k-w]
            possible_rows = counts[d]
        if possible_rows*len(monomials) > entries:
            raise ValueError("coefficient matrix budget exceeded; raise max_matrix_entries explicitly")
        if d not in self._cache:
            fs = [prod(g**int(e) for g, e in zip(self._inputs, m.exponents()[0]))
                  for m in monomials]
            dictionaries = [f.dict() for f in fs]
            support = sorted(set().union(*(set(D) for D in dictionaries)))
            rows = sorted(set(tuple(D.get(e, 0) for D in dictionaries) for e in support))
            M = matrix(QQ, rows) if rows else matrix(QQ, 0, len(monomials))
            B = self._scale**d*_integral_preimage(M)
            B.set_immutable()
            basis = tuple(algebra.domain(sum(B[i, j]*m for i, m in enumerate(monomials)))
                          for j in range(B.ncols()))
            self._cache[d] = CharacterInvariantDegree(d, monomials, B, basis, algebra)
        return self._cache[d]


def f4_character_source(parabolic, *, representation="minimal"):
    """Full F4/P_i Levi character source for the minimal or adjoint module.

    ``representation="minimal"`` selects V(omega4), dimension 26 (default).
    ``representation="adjoint"`` selects V(omega1), dimension 52.
    The adjoint degree-d source and image are 3^d times the minimal ones.
    Completeness follows from the all-primes torus certificate and the
    published rational invariant generators. No finite Z-generator list
    or integral kernel presentation is claimed; use degree(d).image().
    """
    from .exceptional import exceptional_borel
    if representation not in ("minimal", "adjoint"):
        raise ValueError("representation must be 'minimal' or 'adjoint'")
    if representation == "minimal":
        highest = (0, 0, 0, 1)
        change = [[0, 1, 0, -1], [0, 1, -1, 1], [1, -1, 0, 0], [1, -2, 1, 0]]
    else:
        highest = (1, 0, 0, 0)
        change = [[1, -1, 0, 1], [0, 1, -1, 1], [0, 0, 0, 1], [0, 0, 1, -2]]
    borel = exceptional_borel("F4", parabolic)
    entries = {}
    for info in borel.generators():
        f = borel.generator_polynomial(info.name)
        # The middle-node catalogues repeat the same linear generator.
        if f not in entries.values():
            entries[info.name] = f
    return IntegralCharacterSource(
        springer_character_map("F4", highest),
        GeneratorAlgebra(borel, entries), change)
