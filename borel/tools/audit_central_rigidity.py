"""Exact corroboration of central rigidity; no group enumeration.

The geometric proof, including constructibility and the residue argument,
is in CENTRAL_RIGIDITY_FIBERS.md. These controls do not prove smoothness.
"""
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def main():
    controls = {}
    t, c1, c2 = sp.symbols('t c1 c2')
    z1, z2 = t**2*(t-1), (t-1)**3/t
    log_form = sp.cancel(c1*sp.diff(z1,t)/z1+c2*sp.diff(z2,t)/z2)
    residues = [sp.residue(log_form,t,0), sp.residue(log_form,t,1)]
    assert residues == [2*c1-c2,c1+3*c2]
    assert sp.Matrix([[2,-1],[1,3]]).det() == 7
    at_infinity = sp.cancel(-log_form.subs(t,1/t)/t**2)
    assert sp.residue(at_infinity,t,0) == -3*c1-2*c2
    # Exact differentials have no residue; any pole of f gains an order.
    for m in range(1,8):
        derivative = sp.diff(t**(-m),t)
        assert derivative == -m*t**(-m-1)
        assert sp.residue(derivative,t,0) == 0
    controls['divisorial_residues_and_pole_order_controls'] = True

    # A resonant curve lies in the kernel of z1^3*z2^2.
    r1,r2 = t**2,t**-3
    assert sp.cancel(3*sp.diff(r1,t)/r1+2*sp.diff(r2,t)/r2) == 0
    assert sp.cancel(r1**3*r2**2) == 1
    assert sp.Matrix([[3,2]])*sp.Matrix([2,-3]) == sp.zeros(1,1)
    u,v,q,d = sp.symbols('u v q d', nonzero=True)
    monomial = u**3*v**2
    phi = monomial+1/monomial
    mu = [sp.expand(u*sp.diff(phi,u)),sp.expand(v*sp.diff(phi,v))]
    assert sp.cancel(2*mu[0]-3*mu[1]) == 0
    assert [sp.cancel(f.subs({u:2*t**2,v:t**-3})) for f in mu] == [sp.Rational(189,8),sp.Rational(63,4)]
    assert sp.cancel(phi.subs({u:2*t**2,v:t**-3})) == sp.Rational(65,8)
    assert sp.expand(q*(q-1/q-d)) == q*q-d*q-1
    controls['resonant_subtorus_and_rational_target_hyperplane'] = True

    # Full symbolic block identity, valid for arbitrary symmetric H,C.
    h = sp.symbols('h0:6')
    H = sp.Matrix([[h[0],h[1],h[2]],[h[1],h[3],h[4]],[h[2],h[4],h[5]]])
    B = sp.Matrix(3,2,sp.symbols('b0:6'))
    k = sp.symbols('k0:3')
    C = sp.Matrix([[k[0],k[1]],[k[1],k[2]]])
    P = sp.eye(3).row_join(sp.zeros(3,2))
    J = H.row_join(B)
    Z = sp.zeros(2,3).row_join(sp.eye(2))
    U = B.T.row_join(C)
    F = Z.col_join(U)
    omega = sp.zeros(2).row_join(-sp.eye(2)).col_join(
        sp.eye(2).row_join(sp.zeros(2)))
    assert F.T*omega*F == P.T*J-J.T*P
    controls['universal_mixed_Hessian_isotropy_identity'] = True

    x,y,z,c = sp.symbols('x y z c')
    field = sp.QQ.frac_field(c)
    potential = x**3+(z-1)*y
    equations = [sp.diff(potential,x),sp.diff(potential,y),z*sp.diff(potential,z)-c]
    gb = sp.groebner(equations,y,z,x,domain=field,order='lex')
    expected = sp.groebner([y-c,z-1,x*x],y,z,x,domain=field,order='lex')
    assert gb == expected
    assert gb.reduce(potential)[1] == 0
    assert gb.reduce(z-1)[1] == 0
    assert gb.reduce(x)[1] == x and gb.reduce(x*x)[1] == 0
    controls['generic_central_rigidity_does_not_remove_nilpotents'] = True

    # Positive-dimensional excess is allowed with constant central position.
    potential = x*x/2+z  # y remains free
    gb = sp.groebner([sp.diff(potential,x),sp.diff(potential,y),
                     z*sp.diff(potential,z)-c],z,x,y,domain=field,order='lex')
    assert gb == sp.groebner([z-c,x],z,x,y,domain=field,order='lex')
    assert gb.reduce(y)[1] == y and gb.reduce(potential)[1] == c
    controls['smooth_excess_with_fixed_central_position'] = True

    # Actual standard SL2 Springer fiber: h=diag(i,-i), (z-i)^2=0.
    eq = z*z-2*sp.I*z-1
    ring = sp.groebner([eq],z,extension=sp.I)
    assert sp.expand(eq-(z-sp.I)**2) == 0
    inverse = z-2*sp.I
    assert ring.reduce(z*inverse-1)[1] == 0
    character = z+inverse
    assert ring.reduce(character)[1] == 2*(z-sp.I)
    assert ring.reduce(character**2)[1] == 0
    assert ring.reduce(z-sp.I)[1] != 0
    assert ring.reduce((z-inverse)/2-sp.I)[1] == 0
    controls['SL2_scheme_level_constancy_is_false'] = True

    # Finite central images split the actual algebra, retaining nilpotents.
    p1,p2 = (t-1)**2,(t-2)**3
    modulus = sp.expand(p1*p2)
    e1 = sp.rem(p2*sp.invert(p2,p1),modulus,t)
    e2 = 1-e1
    assert sp.rem(e1*e1-e1,modulus,t) == 0
    assert sp.rem(e2*e2-e2,modulus,t) == 0
    assert sp.rem(e1*e2,modulus,t) == 0
    assert sp.rem(e1,p1,t) == 1 and sp.rem(e1,p2,t) == 0
    n1,n2 = sp.rem((t-1)*e1,modulus,t),sp.rem((t-2)*e2,modulus,t)
    assert n1 != 0 and n2 != 0
    assert sp.rem(n1**2,modulus,t) == 0
    assert sp.rem(n2**2,modulus,t) != 0 and sp.rem(n2**3,modulus,t) == 0
    assert sp.rem(t*sp.invert(t,modulus)-1,modulus,t) == 0
    controls['finite_central_algebra_splits_with_nilpotents_retained'] = True

    # Generic targets do not force the finite central image to be reduced.
    # This is a polynomial potential, NOT an adjoint representation.
    f = x**5+y**5+x**3*y**3
    auxiliary = sp.symbols('auxiliary')
    saturated = sp.groebner([sp.diff(f,x),sp.diff(f,y),
                            1-auxiliary*(25-9*x*y)],auxiliary,x,y)
    local = sp.groebner([p.as_expr() for p in saturated.polys
                        if not p.as_expr().has(auxiliary)],x,y,domain=sp.QQ)
    assert local.reduce(f)[1] == y**5/3
    assert local.reduce(f*f)[1] == 0
    assert local.reduce(x**6)[1] == local.reduce(y**6)[1] == 0
    leading = [p.LM(order=local.order).exponents for p in local.polys]
    standard = [(i,j) for i in range(6) for j in range(6)
                if not any(i>=a and j>=b for a,b in leading)]
    assert len(standard) == 16
    assert local.reduce(f+x**3*y**3/5)[1] == 0

    # Independent 16-dimensional cyclic module, without Groebner reduction.
    basis = [(i,j) for i in range(4) for j in range(4)]
    mx,my = sp.zeros(16),sp.zeros(16)
    for col,(i,j) in enumerate(basis):
        if i<3:
            mx[basis.index((i+1,j)),col] = 1
        elif j==0:
            mx[basis.index((2,3)),col] = -sp.Rational(3,5)
        if j<3:
            my[basis.index((i,j+1)),col] = 1
        elif i==0:
            my[basis.index((3,2)),col] = -sp.Rational(3,5)
    zero = sp.zeros(16)
    assert mx*my == my*mx
    assert 5*mx**4+3*mx**2*my**3 == zero
    assert 5*my**4+3*mx**3*my**2 == zero
    assert (25*sp.eye(16)-9*mx*my).det() == 25**16
    fm = mx**5+my**5+mx**3*my**3
    assert fm.rank() == 1 and fm**2 == zero
    assert fm == -mx**3*my**3/5
    cyclic = sp.eye(16)[:,0]
    assert sp.Matrix.hstack(*[mx**i*my**j*cyclic for i,j in basis]) == sp.eye(16)
    # For phi=z(1+f), mu=phi. At any c!=0: z=c(1-f).
    zm = c*(sp.eye(16)-fm)
    assert zm*(sp.eye(16)+fm) == c*sp.eye(16)
    assert (zm-c*sp.eye(16))**2 == zero
    assert (zm-c*sp.eye(16)).rank() == 1
    # The same obstruction occurs in affine degree four with three variables,
    # so a quartic degree bound alone cannot close the adjoint proof.
    w = sp.symbols('w')
    quartic = x**4+y**4+w**4+x*y*w
    saturated4 = sp.groebner([sp.diff(quartic,v) for v in (x,y,w)]+
        [1-auxiliary*(1+64*x*y*w)],auxiliary,x,y,w)
    local4 = sp.groebner([p.as_expr() for p in saturated4.polys
        if not p.as_expr().has(auxiliary)],x,y,w,domain=sp.QQ)
    assert local4.reduce(quartic)[1] == -w**4
    assert local4.reduce(quartic**2)[1] == 0
    leading4 = [p.LM(order=local4.order).exponents for p in local4.polys]
    standard4 = [(i,j,k) for i in range(5) for j in range(5) for k in range(5)
        if not any(all(a>=b for a,b in zip((i,j,k),lm)) for lm in leading4)]
    assert len(standard4) == 11
    # Independent basis: 1,x,x^2,x^3,y,y^2,y^3,w,w^2,w^3,q.
    matrices = []
    for which in range(3):
        mat = sp.zeros(11)
        start = 1+3*which
        mat[start,0] = 1
        mat[start+1,start] = mat[start+2,start+1] = 1
        mat[10,start+2] = 1
        for other in range(3):
            if other!=which:
                third = 3-which-other
                mat[3+3*third,1+3*other] = -4
        matrices.append(mat)
    X,Y,W = matrices
    assert all(a*b==b*a for a in matrices for b in matrices)
    zero11 = sp.zeros(11)
    assert 4*X**3+Y*W == 4*Y**3+X*W == 4*W**3+X*Y == zero11
    assert (sp.eye(11)+64*X*Y*W).det() == 1
    f4m = X**4+Y**4+W**4+X*Y*W
    assert f4m.rank() == 1 and f4m**2 == zero11
    assert f4m == -W**4
    vector = sp.eye(11)[:,0]
    generated = [vector]+[M**i*vector for M in matrices for i in range(1,4)]+[W**4*vector]
    assert sp.Matrix.hstack(*generated) == sp.eye(11)
    controls['generic_nonreduced_central_image_Milnor_algebra_control'] = True

    report = dict(result='passed',uses_Sage=False,control_count=len(controls),
        controls=controls,uniform_written_theorem=True,
        covers_all_simple_types=True,covers_nonregular_source_points=True,
        semisimple_targets_only=True,constancy_on_reduced_components=True,
        generic_abelianization_image_finite=True,
        finite_scheme_image_and_actual_fiber_decomposition_proved=True,
        finite_rational_hyperplane_existence_proved=True,
        hyperplanes_effectively_computed_in_general=False,
        independent_geometric_routes=['divisorial residues and constructibility',
                                      'isotropic critical image'],
        generic_smoothness_proved=False,scheme_level_constancy_claim=False,
        generic_reduced_central_image_follows_from_trace_identity=False,
        new_control_is_adjoint_counterexample=False,
        quartic_degree_bound_alone_suffices=False,
        new_Levi_coverage=False,literature_priority_claim=False,
        proof_note='CENTRAL_RIGIDITY_FIBERS.md',
        scope='Exact corroboration and controls; the geometric arguments are written proofs, not certified by these tests.')
    paths = [Path(__file__),ROOT/report['proof_note']]
    report['source_sha256'] = {p.relative_to(ROOT).as_posix():
        sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths}
    (ROOT/'central-rigidity-audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__ == '__main__':
    main()
