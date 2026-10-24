"""Research helper: certified unit-matrix elimination, not a public API.

The input is a polynomial in affine coordinates over the fraction field
of a separately specified Laurent ring. Every coefficient is checked to
belong to that Laurent ring. Unsupported shapes and nonunit determinants
raise ValueError; this
helper never removes a boundary divisor by passing to a fraction field.
"""
from sage.all import matrix, vector


def as_laurent(R,value):
    """Check Laurent membership explicitly, also for multivariate Sage bases."""
    if value.parent() is R:return value
    numerator=R(value.numerator());denominator=R(value.denominator())
    if not denominator.is_unit():
        raise ValueError('Coefficient has a nonunit Laurent denominator')
    return numerator*denominator**(-1)


def eliminate(phi, multipliers, paired, laurent_ring):
    P=phi.parent(); R=laurent_ring
    try:
        for coefficient in phi.coefficients(): as_laurent(R,coefficient)
    except (TypeError,ValueError) as error:
        raise ValueError('Character coefficients must be Laurent polynomials') from error
    multipliers=tuple(multipliers); paired=tuple(paired)
    if not multipliers or len(multipliers)!=len(paired):
        raise ValueError('Require equally many nonempty multiplier and paired coordinates')
    selected=multipliers+paired
    if len(set(selected))!=len(selected) or any(v not in P.gens() for v in selected):
        raise ValueError('Selected variables must be distinct affine generators')
    constraints=[phi.derivative(t) for t in multipliers]
    if any(f.derivative(t) for f in constraints for t in multipliers):
        raise ValueError('Multiplier dependence is not affine linear')
    try:
        C=matrix(R,[[as_laurent(R,P.base_ring()(f.derivative(a))) for a in paired] for f in constraints])
    except (TypeError,ValueError) as error:
        raise ValueError('Paired coefficients must lie in the Laurent base ring') from error
    det=C.det()
    if not det.is_unit():
        raise ValueError('Constraint determinant is not a Laurent unit; boundary unresolved')
    f=vector(P,[sum(P(C[i,j])*paired[j] for j in range(len(paired)))-g
                for i,g in enumerate(constraints)])
    if any(g.derivative(a) for g in f for a in selected):
        raise ValueError('Constraint remainder depends on selected coordinates')
    G=phi.subs({t:0 for t in multipliers})
    assert phi==G+sum(t*g for t,g in zip(multipliers,constraints))
    Ci=matrix(R,[[as_laurent(R,v) for v in row] for row in C.inverse().rows()])
    alpha=Ci.change_ring(P)*f
    graph=dict(zip(paired,alpha))
    tau=-Ci.transpose().change_ring(P)*vector(P,[G.derivative(a).subs(graph) for a in paired])
    graph.update(zip(multipliers,tau))
    reduced=G.subs(graph)
    assert all(phi.derivative(v).subs(graph)==0 for v in selected)
    return dict(reduced=reduced,graph=graph,determinant=det,
                exact_graph_verified=True,removed_source_divisors=[])
