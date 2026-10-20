"""Adjoint quotients and differentials of normalized trace Springer maps.

All coordinates use the simply connected weight lattice. No Borel quotient
is taken: positive-degree Weyl invariants are retained here.
"""
from dataclasses import dataclass

from sage.all import (PolynomialRing, QQ, ZZ, RootSystem, WeylCharacterRing,
                      matrix, prod, vector)

from .characters import SpringerCharacterMap, _integer, springer_character_map


def _bounded(f, limit):
    if len(f.dict()) > limit:
        raise ValueError("Laurent support budget exceeded; raise max_terms explicitly")
    return f


def _torus_point(point, rank):
    values = vector(point)
    if len(values) != rank or any(x == 0 for x in values):
        raise ValueError("supply one nonzero torus coordinate per Dynkin node")
    if values.base_ring() is ZZ:
        values = values.change_ring(QQ)
    return tuple(values)


def _evaluate_laurent(f, values):
    # Sage does not implicitly reduce QQ coefficients into finite fields.
    field = values[0].parent()
    if field.characteristic():
        return sum(field(c)*prod(t**a for t,a in zip(values,e)) for e,c in f.dict().items())
    return f(*values)


@dataclass(frozen=True)
class PowerSumPullback:
    """Exact sum of powers of Laurent forms, without expanding the powers."""
    linear_forms: tuple
    degree: int
    normalization: object

    def evaluate(self, point):
        values = _torus_point(point, self.linear_forms[0].parent().ngens())
        return sum(_evaluate_laurent(f, values)**self.degree for f in self.linear_forms)/values[0].parent()(self.normalization)

    def expand(self, *, max_terms=50000):
        limit = _integer(max_terms, "max_terms", 1)
        L = self.linear_forms[0].parent()
        result = L.zero()
        for f in self.linear_forms:
            power = L.one()
            for _ in range(self.degree):
                power = _bounded(power*f, limit)
            result = _bounded(result+power, limit)
        return result/self.normalization


@dataclass(frozen=True)
class FactoredCayleyDeterminant:
    """Full determinant retained as the Cartan block and every root factor.

    evaluate() uses exact matrix/scalar operations and avoids Laurent
    expansion. torus_character_polynomial, when present, is the compact
    fundamental-character expression for the Cartan-block determinant.
    """
    torus_jacobian_matrix: object
    root_factors: tuple
    torus_character_polynomial: object = None

    def evaluate(self, point):
        values = _torus_point(point, self.torus_jacobian_matrix.nrows())
        J = matrix([[_evaluate_laurent(f, values) for f in row] for row in self.torus_jacobian_matrix])
        result = J.det()
        for alpha, factor in self.root_factors:
            result *= prod(t**a for t, a in zip(values, alpha))*_evaluate_laurent(factor, values)**2
        return result

    def expand(self, *, max_terms=50000):
        limit = _integer(max_terms, "max_terms", 1)
        result = _bounded(self.torus_jacobian_matrix.det(), limit)
        ts = result.parent().gens()
        for alpha, factor in self.root_factors:
            result = _bounded(result*prod(t**a for t, a in zip(ts, alpha))*factor*factor, limit)
        return result


def _exact_laurent_quotient(f, g):
    """Polynomial division after removing each Laurent monomial shift."""
    L = f.parent()
    n = L.ngens()
    def shifted(h):
        D = h.dict()
        shift = tuple(min(e[i] for e in D) for i in range(n))
        P = L.polynomial_ring()({tuple(e[i]-shift[i] for i in range(n)): c
                                 for e, c in D.items()})
        return P, shift
    if not f:
        return L.zero()
    F, a = shifted(f)
    G, b = shifted(g)
    q, r = F.quo_rem(G)
    if r:
        raise ArithmeticError("root-direction factor did not divide exactly")
    return L(q)*prod(t**(x-y) for t, x, y in zip(L.gens(), a, b))


@dataclass(frozen=True)
class AdjointQuotient:
    """Basic Lie invariants and their pullbacks in fundamental characters.

    ``character_polynomials`` use chi_i = character of V(omega_i).
    ``torus_pullbacks`` are Laurent polynomials, or PowerSumPullback objects
    in the E6/E7 helpers to avoid enormous expansions. These are
    characteristic-zero quotient coordinates, not integral-source
    generators. Basicness is checked by degrees and algebraic independence.
    """
    invariants: tuple
    degrees: tuple
    torus_pullbacks: tuple
    character_polynomials: tuple

    def character_jacobian(self):
        """Jacobian in quotient coordinates chi_i, distinct from d theta."""
        xs = self.character_polynomials[0].parent().gens()
        return matrix([[f.derivative(x) for x in xs] for f in self.character_polynomials])

    def generic_degree(self, *, max_rank=2):
        """Exact generic fiber length over QQ(a_1,...,a_rank).

        Uses Groebner bases over a rational function field, rather than a
        specialized fiber or a numerical root count. The rank guard bounds
        scope, not running time; high-rank elimination may be impractical.
        """
        n = len(self.invariants)
        if n > _integer(max_rank, "max_rank", 1):
            raise ValueError("rank budget exceeded; raise max_rank explicitly")
        K = PolynomialRing(QQ, n, names=["a%d" % (i+1) for i in range(n)]).fraction_field()
        R = PolynomialRing(K, n, names=["x%d" % (i+1) for i in range(n)])
        ideal = R.ideal([R(f(*R.gens()))-a for f, a in zip(self.character_polynomials, K.gens())])
        dimension = ideal.vector_space_dimension()
        if dimension < 1:
            raise ArithmeticError("generic fiber did not have positive finite length")
        return int(dimension)


class SpringerGeometry:
    """Exact geometry attached to an existing SpringerCharacterMap.

    Character conversion and the full determinant can be expensive. Support
    bounds are checked after each multiplication; they are not time/memory
    guarantees for the individual Sage operations.
    """
    def __init__(self, springer_map):
        if not isinstance(springer_map, SpringerCharacterMap):
            raise TypeError("springer_map must be a SpringerCharacterMap")
        self.springer_map = springer_map
        self._cartan = matrix(QQ, springer_map.cartan_type.cartan_matrix())
        self.character_ring = PolynomialRing(QQ, springer_map.rank,
                                            names=["chi%d" % (i+1) for i in range(springer_map.rank)])
        self._fundamentals = {}

    def torus_jacobian(self):
        """d theta in logarithmic torus coordinates; equals identity at 1."""
        s = self.springer_map
        ts = s.torus_ring.gens()
        J = matrix(s.torus_ring, s.rank, s.rank,
                   lambda i, j: ts[j]*s.theta[i].derivative(ts[j]))
        J.set_immutable()
        return J

    def root_factors(self):
        """Return (alpha, r_alpha) with r_alpha=alpha(theta)/(e^alpha-1).

        alpha is in fundamental-weight coordinates. The left-translated
        derivative on g_alpha is e^alpha*r_alpha; the opposite root pair
        contributes e^alpha*r_alpha^2. Every division is checked exactly.
        """
        s = self.springer_map
        ts = s.torus_ring.gens()
        result = []
        for root in RootSystem(s.cartan_type).root_space().positive_roots():
            a = self._cartan*vector(QQ, root.to_vector())
            a = tuple(map(int, a))
            e = prod(t**x for t, x in zip(ts, a))
            h = sum(x*f for x, f in zip(a, s.theta))
            result.append((a, _exact_laurent_quotient(h, e-1)))
        return tuple(result)

    def cayley_determinant(self, *, max_rank=4, max_terms=50000):
        """Full dim(g)-dimensional determinant, not just the torus determinant.

        Left translation identifies T_g G with g; d theta_1=id fixes the
        normalization. Its zero set on T determines the invariant critical
        divisor on G. The rank guard must be raised explicitly for E types.
        """
        rank_bound = _integer(max_rank, "max_rank", 1)
        limit = _integer(max_terms, "max_terms", 1)
        s = self.springer_map
        if s.rank > rank_bound:
            raise ValueError("rank budget exceeded; raise max_rank explicitly")
        f = _bounded(self.torus_jacobian().det(), limit)
        for a, r in self.root_factors():
            f = _bounded(f*prod(t**x for t, x in zip(s.torus_ring.gens(), a))*r*r, limit)
        return f

    def factored_cayley_determinant(self):
        """Return the full determinant as an exact, unexpanded product."""
        return FactoredCayleyDeterminant(self.torus_jacobian(), self.root_factors())

    def _is_invariant(self, D):
        n = self.springer_map.rank
        D = {tuple(e): c for e, c in D.items()}
        for i in range(n):
            reflected = {tuple(e[j]-e[i]*int(self._cartan[j, i]) for j in range(n)): c
                         for e, c in D.items()}
            if reflected != D:
                return False
        return True

    def character_coordinates(self, laurent, *, max_terms=50000, max_steps=10000,
                              max_character_dimension=100000):
        """Express a W-invariant Laurent polynomial in fundamental characters.

        Highest-weight triangular subtraction, with exact reconstruction.
        Simply connectedness is essential to these polynomial coordinates.
        """
        s = self.springer_map
        limit = _integer(max_terms, "max_terms", 1)
        steps = _integer(max_steps, "max_steps", 1)
        dimension_bound = _integer(max_character_dimension, "max_character_dimension", 1)
        f = _bounded(s.torus_ring(laurent), limit)
        if not self._is_invariant(f.dict()):
            raise ValueError("Laurent input must be Weyl invariant")
        WCR = WeylCharacterRing(s.cartan_type, style="coroots")
        coroots = RootSystem(s.cartan_type).ambient_space().simple_coroots()
        height = vector(QQ, [1]*s.rank)*self._cartan.inverse()
        remainder, answer = f, self.character_ring.zero()
        powers = [{0: s.torus_ring.one()} for _ in range(s.rank)]
        for step in range(steps+1):
            if not remainder:
                assert s.torus_ring(answer(*[self._fundamentals.get(i, s.torus_ring.zero())
                                              for i in range(s.rank)])) == f
                return answer
            if step == steps:
                raise ValueError("character conversion budget exceeded; raise max_steps explicitly")
            D = remainder.dict()
            dominant = [e for e in D if all(x >= 0 for x in e)]
            if not dominant:
                raise ArithmeticError("Weyl-invariant remainder has no dominant weight")
            e = max(dominant, key=lambda a: (height.dot_product(vector(a)), tuple(a)))
            term = s.torus_ring.one()
            for i, exponent in enumerate(e):
                if exponent and 1 not in powers[i]:
                    character = WCR(tuple(int(i == j) for j in range(s.rank)))
                    if character.degree() > dimension_bound:
                        raise ValueError("fundamental character dimension budget exceeded; "
                                         "raise max_character_dimension explicitly")
                    if i not in self._fundamentals:
                        self._fundamentals[i] = _bounded(s.torus_ring({
                            tuple(int(w.inner_product(coroots[j+1])) for j in range(s.rank)): QQ(m)
                            for w, m in character.weight_multiplicities().items()}), limit)
                    powers[i][1] = self._fundamentals[i]
                for k in range(2, int(exponent)+1):
                    if k not in powers[i]:
                        powers[i][k] = _bounded(powers[i][k-1]*self._fundamentals[i], limit)
                term = _bounded(term*powers[i][int(exponent)], limit)
            coefficient = D[e]
            answer += coefficient*prod(x**a for x, a in zip(self.character_ring.gens(), e))
            remainder = _bounded(remainder-coefficient*term, limit)

    def adjoint_quotient(self, invariants, *, max_terms=50000):
        """Certify homogeneous basic invariants and compute their pullbacks.

        Checks W invariance, fundamental degrees and a nonzero Jacobian.
        Matching Hilbert series then certifies generation over QQ.
        """
        s = self.springer_map
        fs = tuple(s.ring(f) for f in invariants)
        if len(fs) != s.rank or any(not f or not f.is_homogeneous() for f in fs):
            raise ValueError("supply rank-many nonzero homogeneous basic invariants")
        expected = RootSystem(s.cartan_type).root_space().weyl_group().degrees()
        degrees = tuple(int(f.total_degree()) for f in fs)
        if sorted(degrees) != sorted(expected):
            raise ValueError("invariant degrees must be the fundamental Weyl degrees")
        ws = s.ring.gens()
        for i in range(s.rank):
            images = list(ws)
            images[i] -= sum(self._cartan[j, i]*ws[j] for j in range(s.rank))
            if any(f(*images) != f for f in fs):
                raise ValueError("Lie inputs must be Weyl invariant")
        if not matrix([[f.derivative(w) for w in ws] for f in fs]).det():
            raise ValueError("Lie invariants must be algebraically independent")
        pulled = tuple(s.pullback(f) for f in fs)
        coordinates = tuple(self.character_coordinates(f, max_terms=max_terms) for f in pulled)
        return AdjointQuotient(fs, degrees, pulled, coordinates)


def g2_adjoint_quotient(*, representation="minimal"):
    """G2 quotient coordinates for the 7-dimensional or adjoint trace map.

    The fixed target generators are q=3w1^2-3w1w2+w2^2 and
    z=w1^2(w1-w2)^2(2w1-w2)^2, of degrees 2 and 6.
    """
    if representation not in ("minimal", "adjoint"):
        raise ValueError("representation must be 'minimal' or 'adjoint'")
    s = springer_character_map("G2", (1, 0) if representation == "minimal" else (0, 1))
    w, v = s.ring.gens()
    return SpringerGeometry(s).adjoint_quotient((3*w*w-3*w*v+v*v,
                                                w*w*(w-v)**2*(2*w-v)**2))


def quadratic_character_pullback(cartan_type, highest_weight, *, max_dimension=1000):
    """Pull back q(h)=(h,h)/2 by an irreducible trace map, via Casimirs.

    Long roots have squared length two. The return value is a rational Weyl
    character-ring element, expressed in irreducible characters, not the
    fundamental-character polynomial coordinates of AdjointQuotient.
    This computes only the degree-two invariant. Tensor-square decomposition
    is exact but can be costly; the input dimension is guarded before it runs.
    """
    from .character_lattices import irreducible_character_lattice
    bound = _integer(max_dimension, "max_dimension", 1)
    info = irreducible_character_lattice(cartan_type, highest_weight)
    if info.dimension > bound:
        raise ValueError("representation dimension budget exceeded; raise max_dimension explicitly")
    kind, weight, n = info.cartan_type, info.highest_weight, len(info.highest_weight)
    WCR = WeylCharacterRing(kind, base_ring=QQ, style="coroots")
    character = WCR(weight)
    weight_metric = (info.trace_gram/info.trace_index).inverse()
    def casimir(a):
        a = vector(QQ, a)
        return a*weight_metric*(a+vector(QQ, [2]*n))
    c = casimir(weight)
    coroots = WCR.space().simple_coroots()
    answer = WCR.zero()
    for nu, multiplicity in (character*character).monomial_coefficients().items():
        a = tuple(int(nu.inner_product(coroots[j+1])) for j in range(n))
        answer += multiplicity*(casimir(a)-2*c)/(4*info.trace_index**2)*WCR(a)
    return answer
