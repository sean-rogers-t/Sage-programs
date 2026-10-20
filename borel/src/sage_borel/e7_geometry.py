"""Certified higher E7 quotient coordinates and full factored determinants."""
from functools import lru_cache
from importlib.resources import files
import json

from sage.all import QQ, prod

from .characters import springer_character_map
from .geometry import (AdjointQuotient, FactoredCayleyDeterminant,
                       PowerSumPullback, SpringerGeometry, _torus_point)


@lru_cache(maxsize=1)
def _data():
    return json.loads(files('sage_borel').joinpath('data/e7_geometry.json').read_text())


def _representation(name):
    if name not in ('minimal', 'adjoint'):
        raise ValueError("representation must be 'minimal' or 'adjoint'")
    return (0,0,0,0,0,0,1) if name == 'minimal' else (1,0,0,0,0,0,0)


def _polynomial(R, record):
    return R({tuple(e): QQ(c)/QQ(record['scale']) for e,c in record['terms']})


def e7_fundamental_characters(point):
    """Evaluate chi_1,...,chi_7 in Bourbaki/Sage numbering.

    Uses 56 minuscule weights, 126 roots, the 912-dimensional character,
    and exterior characters. Supports QQ and finite fields of characteristic
    greater than 3, provided all torus coordinates are nonzero.
    """
    values = _torus_point(point, 7)
    data = _data()
    monomial = lambda w: prod(t**a for t,a in zip(values,w))
    weights = [monomial(w) for w in data['weights_56']]
    p1,p2,p3,p4 = (sum(v**k for v in weights) for k in (1,2,3,4))
    a1 = 7+sum(monomial(w) for w in data['roots'])
    a2 = 7+sum(monomial(tuple(2*a for a in w)) for w in data['roots'])
    e2 = (p1*p1-p2)/2
    e3 = (p1**3-3*p1*p2+2*p3)/6
    e4 = (p1**4-6*p1*p1*p2+3*p2*p2+8*p1*p3-6*p4)/24
    return (a1, sum(m*monomial(w) for w,m in data['weights_912']),
            (a1*a1-a2)/2-a1, e4-e2, e3-p1, e2-1, p1)


@lru_cache(maxsize=1)
def _target_invariants():
    # The same target coordinates are used for both source representations.
    s = springer_character_map('E7', _representation('minimal'))
    forms = tuple(sum(a*w for a,w in zip(mu,s.ring.gens())) for mu in _data()['weights_56'])
    invariants = tuple(sum(f**d for f in forms)/QQ(24 if d == 2 else 1)
                       for d in _data()['degrees'])
    return forms, invariants


def e7_adjoint_quotient(*, representation='minimal'):
    """All seven E7 quotient coordinates for the 56- or 133-dimensional map.

    Target: p_2/24,p_6,p_8,p_10,p_12,p_14,p_18, with p_d(h)=tr_56(h^d).
    Returns certified exact fundamental-character polynomials. Torus pullbacks
    stay as PowerSumPullback objects with evaluate() and guarded expand().
    This constructor does not compute the generic degree of the E7 map.
    """
    weight = _representation(representation)
    data = _data()
    s = springer_character_map('E7', weight)
    g = SpringerGeometry(s)
    forms, invariants = _target_invariants()
    pulled_forms = tuple(s.pullback(f) for f in forms)
    pulled, coordinates = [], []
    for record in data['representations'][representation]['coordinates']:
        degree = int(record['degree'])
        normalization = QQ(24 if degree == 2 else 1)
        pulled.append(PowerSumPullback(pulled_forms, degree, normalization))
        coordinates.append(_polynomial(g.character_ring,record)/normalization)
    return AdjointQuotient(invariants, tuple(data['degrees']), tuple(pulled), tuple(coordinates))


def e7_cayley_determinant(*, representation='minimal'):
    """Full 133-dimensional Cayley determinant, retained as an exact product.

    Includes the 7-dimensional Cartan block and all 63 opposite-root pairs.
    Also supplies the Cartan determinant as a fundamental-character polynomial.
    """
    s = springer_character_map('E7', _representation(representation))
    geometry = SpringerGeometry(s)
    base = geometry.factored_cayley_determinant()
    polynomial = _polynomial(geometry.character_ring,
        _data()['representations'][representation]['torus_determinant'])
    return FactoredCayleyDeterminant(base.torus_jacobian_matrix, base.root_factors, polynomial)
