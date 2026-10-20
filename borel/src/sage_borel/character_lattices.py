"""Closed linear Springer lattices for irreducible highest-weight modules.

Only the root system and highest weight are used. No representation weights
or Laurent expansions are enumerated. See IRREDUCIBLE_CHARACTER_LATTICES.md.
"""
from collections import deque
from dataclasses import dataclass

from sage.all import (CartanMatrix, CartanType, PolynomialRing, QQ, RootSystem,
                      ZZ, diagonal_matrix, lcm, matrix, prod, vector)

from .characters import _integer


@dataclass(frozen=True)
class IrreducibleCharacterLattice:
    """Exact source lattice for one nontrivial irreducible representation.

    Matrix columns use fundamental weights. ``trace_index`` uses long roots
    of squared length two, so the adjoint index is twice the dual Coxeter
    number. ``lattice_type`` names the exceptional lattice before scaling by
    the index; ``dual(Q+Z*lambda)`` describes other types via the matrix formula.
    This object does not construct the Springer Laurent formulas.
    """
    cartan_type: object
    highest_weight: tuple
    dimension: int
    trace_index: object
    central_character_order: int
    lattice_type: str
    trace_gram: object
    moment_lattice_basis: object
    linear_basis_matrix: object
    multiplicity_one_paths: tuple
    ring: object

    @property
    def linear_generators(self):
        """A complete Z-basis as rational fundamental-weight polynomials."""
        B = self.linear_basis_matrix
        return tuple(sum(B[i, j]*self.ring.gen(i) for i in range(B.nrows()))
                     for j in range(B.ncols()))


def _multiplicity_one_paths(C, highest):
    """Shortest paths from the support of lambda to each target node."""
    n = len(highest)
    paths = []
    for target in range(n):
        pending = deque([(target,)])
        while pending:
            path = pending.popleft()
            j = path[-1]
            if highest[j]:
                paths.append(tuple(i+1 for i in reversed(path)))
                break
            pending.extend(path+(i,) for i in range(n)
                           if i != j and C[i, j] and i not in path)
    assert len(paths) == n
    return tuple(paths)


def irreducible_character_lattice(cartan_type, highest_weight):
    """Compute Lambda_lambda from the highest weight without listing weights.

    For any finite simple type, D_lambda=sum Z*(m_mu*mu)=Q+Z*lambda and
    Lambda_lambda=M_lambda*D_lambda^#, where # is ordinary coordinate duality.
    M_lambda is obtained from the Weyl dimension and Dynkin index formulas.

    For G2/F4 this is index*Q_long. For E6/E7 it is index*P when lambda is
    in Q, and index*Q otherwise. For E8 it is index*P=index*Q. Direct sums
    are deliberately not accepted by this irreducible-only interface.
    """
    kind = CartanType(cartan_type)
    if not (kind.is_finite() and kind.is_irreducible()
            and kind.is_crystallographic()):
        raise ValueError("a finite irreducible crystallographic type is required")
    n = int(kind.rank())
    if kind != CartanType([kind.type(), n]) or tuple(kind.index_set()) != tuple(range(1, n+1)):
        raise ValueError("standard Dynkin numbering is required")
    highest = tuple(_integer(a, "highest-weight coordinate") for a in highest_weight)
    if len(highest) != n or not any(highest):
        raise ValueError("a nonzero dominant highest weight of the Cartan rank is required")
    C = matrix(ZZ, CartanMatrix(kind).rows(), sparse=False)
    symmetrizer = CartanMatrix(kind).symmetrizer()
    scale = max(symmetrizer.values())
    lengths = [QQ(symmetrizer[i+1])/scale for i in range(n)]
    # C_ij=<alpha_j,alpha_i^vee>; H is the simple-coroot Gram matrix.
    H = C*diagonal_matrix(QQ, [1/d for d in lengths])
    assert H.is_symmetric()
    positive = [tuple(a.to_vector()) for a in RootSystem(kind).root_lattice().positive_roots()]
    dimension = prod(QQ(sum(d*a*(v+1) for d,a,v in zip(lengths,alpha,highest)))
                     / sum(d*a for d,a in zip(lengths,alpha)) for alpha in positive)
    assert dimension.denominator() == 1
    weight = vector(QQ, highest)
    casimir = weight*(H.inverse()*vector(QQ, [a+2 for a in highest]))
    index = QQ(dimension*casimir/(n+2*len(positive)))
    assert index > 0 and index.denominator() == 1
    central_order = int(lcm(a.denominator() for a in C.inverse()*weight))
    D = C.augment(matrix(ZZ, n, 1, highest)).transpose().hermite_form(
        include_zero_rows=False).transpose()
    M = index*H
    raw_basis = M*D.inverse().transpose()
    assert raw_basis.denominator() == 1
    B = matrix(ZZ, raw_basis).transpose().hermite_form().transpose()
    name = "dual(Q+Z*lambda)"
    if kind.type() in ("G", "F"):
        name = "Q_long"
    elif kind.type() == "E":
        name = "P" if central_order == 1 else "Q"
    for A in (M, D, B):
        A.set_immutable()
    return IrreducibleCharacterLattice(kind, highest, int(dimension), index,
        central_order, name, M, D, B, _multiplicity_one_paths(C, highest),
        PolynomialRing(QQ, n, names=["w{}".format(i+1) for i in range(n)]))
