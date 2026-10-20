"""Certified literal E6/P1 and E7/P7 character sources.

The symmetric-coordinate model for E7 is integral, not merely rational.
Full-source proofs and independent checks are in MINUSCULE_CHARACTER_API.md
and E6_REPRESENTATION_COMPARISON.md. No saved image tables are imported.
"""
from sage.all import (CartanMatrix, PolynomialRing, QQ, ZZ, binomial,
                      matrix, prod, vector)

from .algebra import GeneratorAlgebra
from .character_invariants import IntegralCharacterSource
from .characters import _integer, springer_character_map
from .exceptional import exceptional_borel
from ._fixed_point import FixedPointBorel


class _MinusculeAlgebra(GeneratorAlgebra):
    """Evaluate products at fixed points before ordinary interpolation."""

    def _evaluate_monomial(self, exponents):
        if exponents not in self._image_cache:
            if not hasattr(self, "_fixed_point"):
                b = self.borel_map
                self._fixed_point = FixedPointBorel(b.exceptional_type, b.parabolic)
                self._fixed_values = tuple(self._fixed_point.polynomial_values(f)
                                           for f in self.generator_polynomials)
            loc = self._fixed_point
            d = sum(e*w for e, w in zip(exponents, self.generator_degrees))
            if d > loc.top:
                answer = {}
            else:
                values = [prod(v[i]**e for v, e in zip(self._fixed_values, exponents))
                          for i in range(len(loc.basis))]
                answer = {self.borel_map.W.from_reduced_word(w.reduced_word()): c
                          for w, c in loc.coefficients(values, d).items()}
            self._image_cache[exponents] = self.borel_map._expansion(answer)
        return self._image_cache[exponents]


def _certified_source(s, algebra, change, inputs, scale, weights):
    """Internal constructor for the explicitly proved coefficient models."""
    certificate = s.certify_polynomial_source(change)
    if not certificate.certified:
        raise ValueError("the named source failed its all-primes certificate")
    source = IntegralCharacterSource.__new__(IntegralCharacterSource)
    source.springer_map = s
    source.rational_algebra = algebra
    source.torus_certificate = certificate
    source._scale = ZZ(scale)
    source._inputs = tuple(inputs)
    source._coefficient_weights = tuple(weights)
    source._cache = {}
    return source


def e6_character_source(parabolic=1, *, representation="minimal"):
    """Full character source on E6/P1, for V(omega1) or the adjoint V(omega2).

    The torus sources are Sym_Z(6Q) and Sym_Z(24P), respectively, for root
    lattice Q and weight lattice P of simply connected E6. They are not
    related by uniform scalar rescaling. Only P1 is currently packaged.
    """
    if _integer(parabolic, "parabolic", 1) != 1:
        raise ValueError("e6_character_source currently supports only P1")
    if representation not in ("minimal", "adjoint"):
        raise ValueError("representation must be 'minimal' or 'adjoint'")
    C = CartanMatrix("E6")
    b = exceptional_borel("E6", 1)
    names = ("y1", "e1", "e2", "euler", "e3", "e4")
    algebra = _MinusculeAlgebra(b, {n: b.generator_polynomial(n) for n in names})
    if representation == "minimal":
        s = springer_character_map("E6", (1, 0, 0, 0, 0, 0))
        U = matrix(ZZ, [[0,-2,-1,-3,-2,-1], [1,-1,1,0,0,0],
            [0,-2,-1,-2,-2,-1], [0,-2,-1,-2,-1,-1],
            [0,-2,-1,-2,-1,0], [0,-3,-1,-3,-2,-1]])
        change = U*(6*C)*s.linear_basis_matrix.transpose().inverse()
        substitutions = C.inverse()*vector(b.weights)
        inputs = tuple(f(*substitutions) for f in algebra.generator_polynomials)
        scale = 6
    else:
        s = springer_character_map("E6", (0, 1, 0, 0, 0, 0))
        change = [[1,-1,-1,1,0,0], [0,0,0,1,-1,0], [0,0,0,0,1,-1],
                  [0,1,0,-1,0,1], [0,1,0,-1,0,0], [0,0,1,-1,0,0]]
        assert s.linear_basis_matrix == 24*matrix.identity(ZZ, 6)
        inputs = algebra.generator_polynomials
        scale = 24
    return _certified_source(s, algebra, change, inputs, scale, (1,)*6)


def _e7_symmetric_inputs():
    """Newton formulas in Z[e1(v),...,e6(v),u], the integral A5 model."""
    R = PolynomialRing(QQ, ["E"+str(i) for i in range(1, 7)]+["u"])
    es = R.gens()[:6]
    u = R.gen(6)
    E = es[0]
    z = 3*E-6*u
    power = [R(6)]
    for k in range(1, 13):
        value = sum((-1)**(i-1)*es[i-1]*power[k-i] for i in range(1, min(k, 7)))
        if k <= 6:
            value += (-1)**(k-1)*k*es[k-1]
        power.append(R(value))
    X = [R(sum(binomial(k,j)*6**j*(-E)**(k-j)*power[j] for j in range(k+1)))
         for k in range(13)]
    inputs = [(3*u-E)/2]
    for d in (2, 5, 6, 8, 9, 12):
        first = 2*sum(binomial(d,j)*X[d-j]*z**j for j in range(0, d+1, 2))
        pairs = QQ((-1)**d)/2*(sum(binomial(d,j)*X[j]*X[d-j] for j in range(d+1))-2**d*X[d])
        inputs.append(R(first+pairs))
    return tuple(inputs)


def e7_character_source(parabolic=7, *, representation="minimal"):
    """Full minimal or adjoint character source on E7/P7.

    V(omega7), dimension 56, has torus form Sym_Z(12Q); the adjoint
    V(omega1), dimension 133, has form Sym_Z(36P). Uses integral A5
    symmetric coordinates for Q or P. Only P7 is currently packaged.
    """
    if _integer(parabolic, "parabolic", 1) != 7:
        raise ValueError("e7_character_source currently supports only P7")
    if representation not in ("minimal", "adjoint"):
        raise ValueError("representation must be 'minimal' or 'adjoint'")
    inputs = _e7_symmetric_inputs()
    if representation == "minimal":
        s = springer_character_map("E7", (0, 0, 0, 0, 0, 0, 1))
        C = CartanMatrix("E7")
        U = matrix(ZZ, [[-1,-2,-3,-4,-3,-2,-1], [0,0,-1,0,0,0,0],
            [0,-1,-1,-1,-1,-1,-1], [0,-1,-1,-1,-1,-1,0],
            [0,-1,-1,-1,-1,0,0], [0,-1,-1,-1,0,0,0], [-1,-3,-4,-4,-3,-2,-1]])
        change = U*(12*C)*s.linear_basis_matrix.transpose().inverse()
        scale = 12
    else:
        s = springer_character_map("E7", (1, 0, 0, 0, 0, 0, 0))
        assert s.linear_basis_matrix == 36*matrix.identity(ZZ, 7)
        change = [[1,-1,-1,1,0,0,0], [0,0,0,1,-1,0,0],
                  [0,0,0,0,1,-1,0], [0,0,0,0,0,1,-1],
                  [0,1,0,-1,0,0,1], [0,1,0,-1,0,0,0], [0,0,1,-1,0,0,0]]
        # P has integral basis v1,...,v6,z, with z=(u-E1)/2.
        # A5 permutes the v's and fixes z; u=2z+E1.
        R = inputs[0].parent()
        es, z = R.gens()[:6], R.gen(6)
        inputs = tuple(f(*(es+(2*z+es[0],))) for f in inputs)
        scale = 36
    b = exceptional_borel("E7", 7)
    names = ("omega7", "psi2", "psi5", "psi6", "psi8", "psi9", "psi12")
    algebra = _MinusculeAlgebra(b, {n: b.generator_polynomial(n) for n in names})
    return _certified_source(s, algebra, change, inputs, scale,
                             (1, 2, 3, 4, 5, 6, 1))
