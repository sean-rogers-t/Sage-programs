"""Certified higher E6 quotient coordinates and full factored determinants."""
from functools import lru_cache
from importlib.resources import files
import json

from sage.all import QQ, prod

from .characters import springer_character_map
from .degrees import exceptional_springer_degree
from .geometry import (AdjointQuotient, FactoredCayleyDeterminant,
                       PowerSumPullback, SpringerGeometry, _torus_point)


@lru_cache(maxsize=1)
def _data():
    return json.loads(files('sage_borel').joinpath('data/e6_geometry.json').read_text())


def _representation(name):
    if name not in ('minimal', 'adjoint'):
        raise ValueError("representation must be 'minimal' or 'adjoint'")
    return (1,0,0,0,0,0) if name == 'minimal' else (0,1,0,0,0,0)


def _polynomial(R, record):
    return R({tuple(e): QQ(c)/QQ(record['scale']) for e,c in record['terms']})


def e6_fundamental_characters(point):
    """Evaluate chi_1,...,chi_6 using 27 weights and exterior/dual characters.

    The numbering is Bourbaki/Sage: chi_3=Lambda^2(27), chi_4=Lambda^3(27),
    chi_5=chi_3^dual, chi_6=chi_1^dual, chi_2=adjoint (dimension 78).
    Rational and finite-field inputs are supported when 2 and 3 are invertible.
    """
    values = _torus_point(point, 6)
    data = _data()
    monomial = lambda w: prod(t**a for t,a in zip(values,w))
    weights = [monomial(w) for w in data['weights_27']]
    p1, p2, p3 = (sum(v**k for v in weights) for k in (1,2,3))
    dual = [1/v for v in weights]
    d1, d2 = (sum(v**k for v in dual) for k in (1,2))
    return (p1, 6+sum(monomial(w) for w in data['roots']), (p1*p1-p2)/2,
            (p1**3-3*p1*p2+2*p3)/6, (d1*d1-d2)/2, d1)


def e6_adjoint_quotient(*, representation='minimal'):
    """All six E6 quotient coordinates for the 27- or 78-dimensional trace map.

    Target: p_2/12,p_5,p_6,p_8,p_9,p_12, where p_d(h)=tr_27(h^d).
    Returns exact fundamental-character polynomials. The torus pullbacks are
    PowerSumPullback objects with evaluate() and explicitly bounded expand().
    generic_degree() uses the proved toric-intersection degree certificate.
    """
    weight = _representation(representation)
    data = _data()
    s = springer_character_map('E6', weight)
    g = SpringerGeometry(s)
    forms = tuple(sum(a*w for a,w in zip(mu,s.ring.gens())) for mu in data['weights_27'])
    pulled_forms = tuple(s.pullback(f) for f in forms)
    invariants, pulled, coordinates = [], [], []
    for record in data['representations'][representation]['coordinates']:
        degree = int(record['degree'])
        normalization = QQ(12 if degree == 2 else 1)
        invariants.append(sum(f**degree for f in forms)/normalization)
        pulled.append(PowerSumPullback(pulled_forms, degree, normalization))
        coordinates.append(_polynomial(g.character_ring,record)/normalization)
    invariants=tuple(invariants)
    degrees=tuple(int(r['degree']) for r in data['representations'][representation]['coordinates'])
    coordinates=tuple(coordinates)
    return AdjointQuotient(invariants,degrees,tuple(pulled),coordinates,
        exceptional_springer_degree('E6',representation=representation),
        (invariants,degrees,coordinates))


def e6_cayley_determinant(*, representation='minimal'):
    """Full 78-dimensional Cayley determinant, retained as an exact product.

    Contains the 6-dimensional Cartan block and all 36 opposite-root pairs.
    The Cartan determinant is also provided in fundamental characters.
    No expansion of the full root product is required for evaluation.
    """
    s = springer_character_map('E6', _representation(representation))
    geometry = SpringerGeometry(s)
    base = geometry.factored_cayley_determinant()
    polynomial = _polynomial(geometry.character_ring,
        _data()['representations'][representation]['torus_determinant'])
    return FactoredCayleyDeterminant(base.torus_jacobian_matrix, base.root_factors, polynomial)
