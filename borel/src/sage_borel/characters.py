"""Integral polynomial characters for normalized trace-projection maps.

The group is simply connected, semisimple and of finite irreducible type.
Weights are fundamental-weight coordinates. A successful leading certificate
proves the full torus intersection in ALL degrees; an inconclusive certificate
or a bounded search without a relation makes no such assertion.
"""
from dataclasses import dataclass
from math import comb
from operator import index

from sage.all import (CartanType, GF, LaurentPolynomialRing, PolynomialRing,
                      QQ, RootSystem, WeylCharacterRing, ZZ, identity_matrix,
                      matrix, prod, vector)
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing_generic


def _integer(value, name, minimum=0):
    if isinstance(value, bool):
        raise ValueError("{} must be an integer".format(name))
    try:
        value = index(value)
    except (TypeError, ValueError) as error:
        raise ValueError("{} must be an integer".format(name)) from error
    if value < minimum:
        raise ValueError("{} must be at least {}".format(name, minimum))
    return int(value)


def _integral_laurents(pullbacks):
    fs = tuple(pullbacks)
    if (not fs or not hasattr(fs[0], "parent")
            or not isinstance(fs[0].parent(), LaurentPolynomialRing_generic)):
        raise TypeError("supply nonempty Laurent polynomials over QQ or ZZ")
    parent = fs[0].parent()
    if parent.base_ring() not in (QQ, ZZ):
        raise ValueError("Laurent coefficients must be exact rationals or integers")
    if any(not hasattr(f, "parent") or f.parent() != parent for f in fs):
        raise ValueError("all pullbacks must have the same Laurent parent")
    if any(QQ(c).denominator() != 1 for f in fs for c in f.dict().values()):
        raise ValueError("pullbacks must have integral coefficients")
    if any(not f for f in fs):
        raise ValueError("pullbacks must be nonzero")
    return fs


@dataclass(frozen=True)
class LeadingTermCertificate:
    """Finite sufficient certificate for independence over every prime field.

    ``cases`` contains (characteristic, leading exponents, coefficients,
    exponent rank). Characteristic zero also covers every prime not listed.
    Rank is computed over QQ, including for positive-characteristic cases.
    """
    certified: bool
    exceptional_primes: tuple
    cases: tuple
    basis_change: tuple


def polynomial_source_certificate(pullbacks, basis_change=None):
    """Certify an integral Laurent polynomial subring by leading monomials.

    ``basis_change`` is an optional unimodular integer row change on the inputs.
    False means this sufficient test is inconclusive, NOT a dependence claim.
    """
    fs = _integral_laurents(pullbacks)
    n = len(fs)
    U = identity_matrix(ZZ, n)
    if basis_change is not None:
        raw = matrix(QQ, basis_change)
        if raw.dimensions() != (n, n) or raw.denominator() != 1:
            raise ValueError("basis_change must be a square integer matrix")
        U = matrix(ZZ, raw)
        if abs(U.det()) != 1:
            raise ValueError("basis_change must be unimodular")
    transformed = tuple(sum(U[i, j]*fs[j] for j in range(n)) for i in range(n))
    if any(not f for f in transformed):
        raise ValueError("basis change produced a zero pullback")
    dictionaries = [f.dict() for f in transformed]
    initial = [D[max(D)] for D in dictionaries]
    primes = tuple(sorted(set(p for c in initial
                              for p in ZZ(c).prime_divisors())))
    cases = []
    for p in (0,) + primes:
        exponents, coefficients = [], []
        for D in dictionaries:
            support = [e for e, c in D.items() if not p or ZZ(c) % p]
            if not support:
                exponents.append((0,)*fs[0].parent().ngens())
                coefficients.append(0)
            else:
                e = max(support)
                exponents.append(tuple(map(int, e)))
                coefficients.append(int(D[e]))
        rank = int(matrix(QQ, exponents).rank())
        cases.append((int(p), tuple(exponents), tuple(coefficients), rank))
    certified = all(case[3] == n and all(case[2]) for case in cases)
    return LeadingTermCertificate(certified, tuple(map(int, primes)), tuple(cases),
                                  tuple(tuple(map(int, row)) for row in U))


def _exponents(n, degree):
    if n == 1:
        yield (degree,)
    else:
        for first in range(degree + 1):
            for tail in _exponents(n - 1, degree - first):
                yield (first,) + tail


def modular_character_relation(pullbacks, prime, max_degree=3,
                               max_monomials=1000, max_terms=100000,
                               max_matrix_entries=2000000):
    """Return an integer lift of a nonzero mod-p relation, or None.

    The relation uses variables u1,... in input order and includes all degrees
    through max_degree. None asserts only absence within that finite bound.
    For arbitrary inputs this certifies a relation, not characteristic-zero
    independence. SpringerCharacterMap supplies the latter separately.
    Resource bounds cap monomials, the union of expansion supports (checked
    after each expansion), and the number of coefficient-matrix entries.
    """
    fs = _integral_laurents(pullbacks)
    p = _integer(prime, "prime", 2)
    if not ZZ(p).is_prime():
        raise ValueError("prime must be prime")
    d = _integer(max_degree, "max_degree")
    max_monomials = _integer(max_monomials, "max_monomials", 1)
    max_terms = _integer(max_terms, "max_terms", 1)
    max_matrix_entries = _integer(max_matrix_entries, "max_matrix_entries", 1)
    n = len(fs)
    if comb(n+d, n) > max_monomials:
        raise ValueError("monomial budget exceeded; raise max_monomials explicitly")
    powers = [e for k in range(d+1) for e in _exponents(n, k)]
    polys, support = [], set()
    for e in powers:
        f = prod(g**a for g, a in zip(fs, e))
        D = f.dict()
        support.update(D)
        if len(support) > max_terms:
            raise ValueError("Laurent support budget exceeded; raise max_terms explicitly")
        if len(support)*len(powers) > max_matrix_entries:
            raise ValueError("coefficient matrix budget exceeded; raise max_matrix_entries explicitly")
        polys.append(D)
    support = sorted(support)
    A = matrix(GF(p), [[D.get(e, 0) for D in polys] for e in support])
    kernel = A.right_kernel()
    if not kernel.dimension():
        return None
    coefficients = list(map(ZZ, kernel.basis()[0]))
    R = PolynomialRing(QQ, n, names=["u{}".format(i+1) for i in range(n)])
    relation = R(sum(c*prod(x**a for x, a in zip(R.gens(), e))
                     for c, e in zip(coefficients, powers)))
    assert relation and any(ZZ(c) % p for c in relation.coefficients())
    assert all(ZZ(c) % p == 0 for c in relation(*fs).dict().values())
    return relation


@dataclass(frozen=True)
class CharacterObstruction:
    """A verified extra integral polynomial F(linear_generators)/prime.

    It lies outside the symmetric algebra of the degree-one lattice.
    It need not be homogeneous or invariant under a specified Levi.
    """
    prime: int
    relation: object
    polynomial: object
    pullback: object


class SpringerCharacterMap:
    """Normalized Springer map from a sum of genuine highest-weight modules.

    Construct using springer_character_map(type, highest_weight), or pass a
    nonempty sequence of highest weights here to specify a direct sum.
    The torus character lattice is the full weight lattice of the simply
    connected group. Other central quotients are not silently inferred.
    """

    def __init__(self, cartan_type, highest_weights):
        kind = CartanType(cartan_type)
        if not (kind.is_finite() and kind.is_irreducible()
                and kind.is_crystallographic()):
            raise ValueError("a finite irreducible crystallographic type is required")
        n = int(kind.rank())
        if kind != CartanType([kind.type(), n]) or tuple(kind.index_set()) != tuple(range(1, n+1)):
            raise ValueError("standard Dynkin numbering is required")
        highest_weights = tuple(tuple(_integer(v, "highest-weight coordinate")
                                      for v in weight) for weight in highest_weights)
        if not highest_weights or any(len(w) != n for w in highest_weights):
            raise ValueError("each highest weight must have the Cartan rank")
        if not any(any(w) for w in highest_weights):
            raise ValueError("the representation must be infinitesimally faithful")
        self.cartan_type = kind
        self.highest_weights = highest_weights
        self.rank = n
        self.ring = PolynomialRing(QQ, n, names=["w{}".format(i+1) for i in range(n)])
        self.torus_ring = LaurentPolynomialRing(QQ, n,
                                               names=["t{}".format(i+1) for i in range(n)])
        WCR = WeylCharacterRing(kind, style="coroots")
        # A tuple avoids rank-one WCR(1) being interpreted as the scalar 1.
        character = sum(WCR(w) for w in highest_weights)
        coroots = RootSystem(kind).ambient_space().simple_coroots()
        weights = {}
        for w, m in character.weight_multiplicities().items():
            coordinates = tuple(ZZ(w.inner_product(coroots[i+1])) for i in range(n))
            weights[coordinates] = weights.get(coordinates, ZZ.zero()) + ZZ(m)
        self.weight_multiplicities = tuple(sorted((tuple(map(int, w)), int(m))
                                                 for w, m in weights.items()))
        self.dimension = int(sum(weights.values()))
        gram = matrix(QQ, n, n, lambda i, j: sum(m*w[i]*w[j] for w, m in weights.items()))
        assert gram.det()
        gram.set_immutable()
        self.trace_gram = gram
        ts = self.torus_ring.gens()
        moments = vector(self.torus_ring, [
            sum(m*w[i]*prod(t**e for t, e in zip(ts, w)) for w, m in weights.items())
            for i in range(n)])
        self.theta = tuple(gram.inverse()*moments)
        # Hermite row lattice gives all rational linear combinations with
        # integral coefficients, not just a denominator-cleared sublattice.
        dictionaries = [f.dict() for f in self.theta]
        support = sorted(set().union(*(set(D) for D in dictionaries)))
        A = matrix(QQ, [[D.get(e, 0) for D in dictionaries] for e in support])
        denominator = A.denominator()
        H = matrix(ZZ, denominator*A).hermite_form(include_zero_rows=False)
        assert H.dimensions() == (n, n)
        B = denominator*H.change_ring(QQ).inverse()
        assert (A*B).denominator() == 1
        # d(theta)_1 = identity: differentiating an integral Laurent pullback
        # gives integral fundamental-weight coefficients for its linear input.
        assert B.denominator() == 1
        B.set_immutable()
        self.linear_basis_matrix = B
        self.linear_generators = tuple(sum(B[i, j]*self.ring.gen(i) for i in range(n))
                                       for j in range(n))
        self.linear_pullbacks = tuple(self.pullback(f) for f in self.linear_generators)
        assert all(self.is_integral(f) for f in self.linear_generators)

    def pullback(self, polynomial):
        """Evaluate a rational fundamental-weight polynomial under theta*."""
        return self.torus_ring(self.ring(polynomial)(*self.theta))

    def is_integral(self, polynomial):
        """Exact finite Laurent test; no Levi-invariance assumption."""
        return all(QQ(c).denominator() == 1
                   for c in self.pullback(polynomial).dict().values())

    def certify_polynomial_source(self, basis_change=None):
        """Try the all-primes sufficient certificate for Sym_Z(Lambda_rho)."""
        return polynomial_source_certificate(self.linear_pullbacks, basis_change)

    def modular_obstruction(self, prime, max_degree=3, **budgets):
        """Find a certified extra integral polynomial, if visible in the bound."""
        relation = modular_character_relation(self.linear_pullbacks, prime,
                                              max_degree, **budgets)
        if relation is None:
            return None
        witness = self.ring(relation(*self.linear_generators)/prime)
        pullback = self.pullback(witness)
        assert all(QQ(c).denominator() == 1 for c in pullback.dict().values())
        return CharacterObstruction(int(prime), relation, witness, pullback)


def springer_character_map(cartan_type, highest_weight):
    """Construct the normalized map for one irreducible highest-weight module."""
    return SpringerCharacterMap(cartan_type, (highest_weight,))


def g2_character_algebra(parabolic):
    """Full seven-dimensional G2 character source for P1 or P2, then Borel.

    In x=6*w1, y=2*w2, the integral invariant rings are
    Z[x,y*(x-y)] for P1 and Z[y,x*(3*y-x)] for P2.
    These are literal character forms, not primitive-image rescalings.
    """
    from .algebra import GeneratorAlgebra
    from .exceptional import exceptional_borel
    borel = exceptional_borel("G2", parabolic)
    x, y = 6*borel.weights[0], 2*borel.weights[1]
    if borel.omitted_nodes == (1,):
        return GeneratorAlgebra(borel, {"x": x, "q": y*(x-y)})
    return GeneratorAlgebra(borel, {"y": y, "p": x*(3*y-x)})
