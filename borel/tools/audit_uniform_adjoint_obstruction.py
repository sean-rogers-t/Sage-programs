"""Exact counterexample to a proposed shortcut, not to generic smoothness.

Sp4 adjoint: a singular incidence point can have nonzero target and exactly
the intended Levi centralizer. Check the local scheme, defining matrices,
and independently the ten-dimensional adjoint trace construction.
"""
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
simp = sp.simplify


def coords(M):
    return sp.Matrix([M[0, 0], M[0, 1], M[1, 0], M[1, 1],
                      M[0, 2], M[0, 3], M[1, 3],
                      M[2, 0], M[2, 1], M[3, 1]])


def sp4_basis():
    basis = []
    for i, j in ((0, 0), (0, 1), (1, 0), (1, 1)):
        M = sp.zeros(4)
        M[i, j] = 1
        M[j+2, i+2] = -1
        basis.append(M)
    for block in ('B', 'C'):
        for i, j in ((0, 0), (0, 1), (1, 1)):
            M = sp.zeros(4)
            for a, b in {(i, j), (j, i)}:
                if block == 'B':
                    M[a, b+2] = 1
                else:
                    M[a+2, b] = 1
            basis.append(M)
    assert sp.Matrix.hstack(*(coords(M) for M in basis)) == sp.eye(10)
    J = sp.zeros(2).row_join(sp.eye(2)).col_join(
        (-sp.eye(2)).row_join(sp.zeros(2)))
    assert all(M.T*J+J*M == sp.zeros(4) for M in basis)
    return basis, J


def higher_order_audit():
    """Check formal elimination, not a new generic-smoothness assertion."""
    v, w, y1, y2, eps = sp.symbols('v w y1 y2 eps')
    H = sp.Matrix([[2, 1], [1, 3]])
    q = sp.symbols('q0:6')
    b = sp.Matrix([q[0]*v*v+q[1]*v*w+q[2]*w*w,
                   q[3]*v*v+q[4]*v*w+q[5]*w*w])
    t, s, mixed = sp.symbols('t s mixed')
    y = sp.Matrix([y1, y2])
    cubic = t*(v**3+v*w*w)
    quartic = s*(v**4+w**4)
    W = ((y.T*H*y)[0]/2 + cubic + (y.T*b)[0]/2
         + mixed*v*y1*y2 + y1**3 + quartic + w**3*y2)
    # Solve the degree-two complementary equations directly, independently
    # of the displayed tensor formula. Higher mixed terms cannot affect F4.
    c1, c2 = sp.symbols('c1 c2')
    scaled = {v: eps*v, w: eps*w, y1: eps**2*c1, y2: eps**2*c2}
    equations = [sp.expand(sp.diff(W, yi).subs(scaled)).coeff(eps, 2)
                 for yi in (y1, y2)]
    solution = sp.solve(equations, (c1, c2))
    expected_y2 = -H.inv()*b/2
    assert all(sp.expand(solution[ci]-expected_y2[i]) == 0
               for i, ci in enumerate((c1, c2)))
    substituted = sp.expand(W.subs(scaled).subs(solution))
    corrected = quartic-(b.T*H.inv()*b)[0]/8
    assert sp.expand(substituted.coeff(eps, 3)-cubic) == 0
    assert sp.expand(substituted.coeff(eps, 4)-corrected) == 0

    # Actual adjoint Sp4: expand the ORIGINAL logarithmic potential and
    # solve its central equation before computing the effective quartic.
    x, yy = sp.symbols('x yy')
    delta = 2*sp.I*sp.sqrt(3)
    z0 = (1+sp.I*sp.sqrt(3))/2
    # D^j A(z0)/j!, obtained from the original two Laurent monomials.
    A_series = sum(simp(2**j*z0**2+(-2)**j*z0**-2)*yy**j/sp.factorial(j)
                   for j in range(1, 5))
    Wsp = sp.expand(A_series*(x*x-1)+delta*yy)
    central_eq = sp.expand(sp.diff(Wsp, yy).subs(yy, c1*x*x))
    a2 = sp.solve(central_eq.coeff(x, 2), c1)[0]
    assert simp(a2+delta/4) == 0
    effective = sp.expand(Wsp.subs(yy, a2*x*x))
    assert simp(effective.coeff(x, 4)) == sp.Rational(3, 2)
    assert sp.expand(Wsp.subs(yy, 0)) == 0  # raw quartic is zero
    assert simp(sp.diff(effective, x).coeff(x, 3)) == 6

    # Known E6 A5 excess: all four complementary equations are independent
    # of x2,x4, and have an invertible Hessian at every stated target.
    x1, x2, x3, x4, x5, z, c = sp.symbols('x1 x2 x3 x4 x5 z c')
    phi = x1*x5+z*z+z**-2+(z+z**-1)*x3
    D = lambda f: z*sp.diff(f, z)
    gradients = [sp.diff(phi, xi) for xi in (x1, x5, x3)]+[D(phi)-c]
    point = {x1: 0, x5: 0, x3: c/(2*sp.I), z: sp.I}
    assert all(simp(f.subs(point)) == 0 for f in gradients)
    assert not phi.has(x2, x4)
    Hess = sp.Matrix([[sp.diff(f, xi) for xi in (x1, x5, x3)]+[D(f)]
                      for f in gradients]).subs(point).applyfunc(simp)
    assert Hess == sp.Matrix([[0,1,0,0],[1,0,0,0],
                              [0,0,0,2*sp.I],[0,0,2*sp.I,-8]])
    assert Hess.det() == -4

    # Previously audited altered-character control: exact X elimination.
    u, X, S = sp.symbols('u X S')
    control = -S*S-2+X*X*(1-u*v)+S*X*u*v
    Xsol = sp.solve(sp.diff(control, X), X)[0]
    eliminated = -S*S-2-S*S*u*u*v*v/(4*(1-u*v))
    assert simp(control.subs(X, Xsol)-eliminated) == 0
    residual = eliminated+S*S+2
    assert sp.series(residual.subs({u: eps*u, v: eps*v}), eps, 0, 5).removeO() == -S*S*u*u*v*v*eps**4/4
    S_z = z**3+z**-3
    assert simp(D(D(-S_z**2-2))+36*(z**6+z**-6)) == 0
    # The perturbation and its central derivative start at degree four,
    # so the implicit central shift starts at degree four and alters F
    # only from degree eight. No scalar specialization of nilpotents.
    central_perturbation = D(residual.subs(S, S_z))
    low = sp.series(central_perturbation.subs({u: eps*u, v: eps*v}), eps, 0, 4).removeO()
    assert low == 0

    # A passing four-jet is not an all-orders certificate.
    delayed = yy*yy/2+v**6
    residual_delayed = delayed.subs(yy, sp.solve(sp.diff(delayed, yy), yy)[0])
    assert sp.Poly(residual_delayed, v).degree() == 6
    assert sp.diff(residual_delayed, v) == 6*v**5
    return dict(corrected_quartic_symbolic_check=True,
                Sp4_raw_quartic='0', Sp4_effective_quartic='3*x^4/2',
                Sp4_implicit_expansion_matches_exact_local_ideal=True,
                E6_excess_residual_constant_to_all_orders=True,
                altered_Levi_control_effective_quartic='-S(z0)^2*u^2*v^2/4',
                four_jet_vanishing_not_sufficient=True,
                generic_bad_image_bound_proved=False)


def generic_target_audit():
    """Keep target parameters generic, without making source parameters generic."""
    x, z, c, u, v, inv = sp.symbols('x z c u v inv')
    field = sp.QQ.frac_field(c)

    # The SAME Sp4 model: the actual quotient fiber over Q(c).
    # Multiplication by z-powers is harmless: the ideal makes z a unit.
    f1 = (z**4+z*z+1)*x
    f2 = 2*(z**4-1)*(x*x-1)-c*z*z
    original = sp.groebner([f1, f2], x, z, domain=field)
    assert original.is_zero_dimensional
    assert sp.groebner([f1, f2, z], x, z, domain=field).reduce(sp.S.One)[1] == 0
    jac = sp.Matrix([f1, f2]).jacobian([x, z]).det()
    smoothness = sp.groebner([f1, f2, jac], x, z, domain=field)
    assert smoothness.reduce(sp.S.One)[1] == 0
    # Independent branch separation: x=0, or A=0. These branches are
    # disjoint over Q(c), since their intersection forces c^2+12=0.
    A = z**4+z*z+1
    P = 2*z**4+c*z*z-2
    assert sp.factor(sp.resultant(A, P, z)) == (c*c+12)**2
    assert sp.gcd(sp.Poly(A*P, z, domain=field),
                  sp.Poly(sp.diff(A*P, z), z, domain=field)).degree() == 0

    # Known E6 excess. Original derivatives, cleared only by Laurent units.
    x1, x2, x3, x4, x5 = sp.symbols('x1 x2 x3 x4 x5')
    variables = (x1, x5, x3, z, x2, x4)
    eq = [x5, x1, z*z+1,
          2*(z**4-1)+(z**3-z)*x3-c*z*z]
    simple = [x1, x5, z*z+1, 2*z*x3-c]
    actual = sp.groebner(eq, *variables, domain=field)
    expected = sp.groebner(simple, *variables, domain=field)
    assert actual == expected
    J = sp.Matrix(simple).jacobian([x1, x5, z, x3])
    assert sp.groebner(simple+[J.det()], *variables, domain=field).reduce(sp.S.One)[1] == 0
    assert actual.reduce(sp.S.One)[1] == 1
    # Inverting a coefficient which vanishes on this whole dominant branch
    # incorrectly deletes it. This is NOT the legitimate target localization.
    erased = sp.groebner(eq+[inv*(z*z+1)-1], inv, *variables, domain=field)
    assert erased.reduce(sp.S.One)[1] == 0

    # Previously proved altered-Levi local model over a generic target:
    # finite separable central image, but an embedded nilpotent at the
    # intersection of two generically smooth components.
    central = 6*z**12+c*z**6-6
    assert sp.gcd(sp.Poly(central, z, domain=field),
                  sp.Poly(sp.diff(central, z), z, domain=field)).degree() == 0
    assert sp.gcd(sp.Poly(central, z, domain=field),
                  sp.Poly(z**6+1, z, domain=field)).degree() == 0
    bad = sp.groebner([u*u*v, u*v*v], u, v, domain=field)
    assert bad.reduce(u*v)[1] == u*v
    assert bad.reduce((u*v)**2)[1] == 0
    assert bad.reduce(u*(u*v))[1] == bad.reduce(v*(u*v))[1] == 0
    # Primary decomposition (u) intersect (v) intersect (u^2,v^2).
    # Compute the latter intersection with (uv), using an elimination variable.
    t = sp.symbols('t')
    intersection = sp.groebner([t*u*v, (1-t)*u*u, (1-t)*v*v],
                               t, u, v, domain=field)
    eliminated = [p.as_expr() for p in intersection.polys
                  if not p.as_expr().has(t)]
    assert sp.groebner(eliminated, u, v, domain=field) == bad
    for unit, other in ((u, v), (v, u)):
        localized = sp.groebner([u*u*v, u*v*v, inv*unit-1],
                                inv, u, v, domain=field)
        assert localized.reduce(other)[1] == 0
    # At the intersection the Jacobian vanishes, with tangent dimension 2
    # and Krull dimension 1 (the two coordinate axes).
    assert sp.Matrix([u*u*v, u*v*v]).jacobian([u, v]).subs({u:0, v:0}) == sp.zeros(2)
    return dict(target_field='Q(c); certificates remain valid over C(c)',
                Sp4_generic_fiber_jacobian_unit=True,
                Sp4_independent_branch_separation=True,
                E6_original_generic_ideal_verified=True,
                E6_smooth_excess_dimension=2,
                source_coefficient_localization_can_erase_dominant_fiber=True,
                altered_Levi_central_polynomial_squarefree=True,
                altered_Levi_components_generically_smooth=True,
                altered_Levi_embedded_nilpotent='u*v; annihilator (u,v)',
                altered_Levi_primary_decomposition_verified=True,
                generic_adjoint_smoothness_proved=False,
                new_Levi_coverage=False)


def reducedness_audit():
    """Separate generic multiplicity, embedded structure and target torsion."""
    x, y, z, c = sp.symbols('x y z c')
    field = sp.QQ.frac_field(c)
    # Existing abstract control: a complete-intersection generic fiber can
    # still be nonreduced at its minimal prime. It has no embedded primes.
    phi = x**3+(z-1)*y
    fiber = [sp.diff(phi, x), sp.diff(phi, y), z*sp.diff(phi, z)-c]
    actual = sp.groebner(fiber, z, y, x, domain=field)
    expected = sp.groebner([z-1, y-c, x*x], z, y, x, domain=field)
    assert actual == expected
    assert actual.reduce(x)[1] == x
    assert actual.reduce(x*x)[1] == 0

    # Target-torsion diagnostic: the relative incidence need not itself
    # be reduced if its nilpotents are confined to special targets.
    vertical_phi = (z-1)**2*x+z
    incidence = sp.groebner([sp.diff(vertical_phi, x)], x, z)
    mu = z*sp.diff(vertical_phi, z)
    epsilon = z-1
    assert incidence.reduce(epsilon)[1] == epsilon
    assert incidence.reduce(epsilon**2)[1] == 0
    assert incidence.reduce(sp.expand((mu-1)*epsilon))[1] == 0
    assert incidence.reduce(sp.expand((mu-1)**2))[1] == 0
    generic = sp.groebner([sp.diff(vertical_phi, x), mu-c], x, z, domain=field)
    assert generic.reduce(sp.S.One)[1] == 0

    # Actual Sp4 relative incidence is reduced: one squarefree equation.
    # The earlier audit nevertheless proves a nonreduced SPECIAL fiber.
    A = z**4+z*z+1
    assert sp.gcd(sp.Poly(A, z), sp.Poly(sp.diff(A, z), z)).degree() == 0
    assert not A.has(x)

    # Actual E6 relative incidence has only three nonzero independent
    # equations. The full generic fiber has four, not six: excess remains.
    x1, x2, x3, x4, x5 = sp.symbols('x1 x2 x3 x4 x5')
    phi6 = x1*x5+z*z+z**-2+(z+z**-1)*x3
    relative = [sp.diff(phi6, xi) for xi in (x1, x2, x3, x4, x5)]
    assert relative == [x5, 0, z+z**-1, 0, x1]
    assert sp.gcd(sp.Poly(z*z+1, z), sp.Poly(2*z, z)).degree() == 0
    return dict(complete_intersection_does_not_force_R0=True,
                existing_abstract_control_generic_ring='Q(c)[x]/(x^2)',
                existing_embedded_control_R0_without_S1=True,
                target_torsion_control_annihilator='c-1',
                target_torsion_control_generic_ring_zero=True,
                Sp4_relative_incidence_squarefree=True,
                E6_relative_incidence_radical=True,
                E6_displayed_full_gradient_not_regular_sequence=True,
                all_adjoint_R0_proved=False,
                all_adjoint_S1_proved=False,
                uniform_generic_reducedness_proved=False,
                new_Levi_coverage=False)


def main():
    x, z, u = sp.symbols('x z u')
    S = z*z+z**-2
    A = S+1
    delta_z = 2*(z*z-z**-2)
    phi = A*x*x-S
    D = lambda f: sp.expand(z*sp.diff(f, z))
    assert simp(D(phi)-delta_z*(x*x-1)) == 0

    # Reconstruct the character from Sym^2 of the defining module.
    diagonal = [z*u, z/u, 1/(z*u), u/z]
    sym_square = (sum(diagonal)**2+sum(t*t for t in diagonal))/2
    assert simp(sym_square-phi.subs(x, u+1/u)) == 0

    z0 = (1+sp.I*sp.sqrt(3))/2
    delta = 2*sp.I*sp.sqrt(3)
    at = {x: 0, z: z0}
    assert simp(A.subs(at)) == 0
    assert simp(delta_z.subs(at)-delta) == 0
    assert simp(D(delta_z).subs(at)) == -4
    assert delta**2 == -12

    # Etale coordinates a=A(z), d=D A(z): d^2=4(a^2-2a-3).
    # Saturation isolates the local branch at (x,a,d)=(0,0,delta).
    a, d, inv = sp.symbols('a d inv')
    field = sp.QQ.algebraic_field(sp.I*sp.sqrt(3))
    equations = [a*x, d*d-4*(a*a-2*a-3),
                 d*(x*x-1)+delta, inv*(a-2)*(x*x-2)-1]
    gb = sp.groebner(equations, inv, d, a, x, domain=field, order='lex')
    eliminated = [p.as_expr() for p in gb.polys if not p.as_expr().has(inv)]
    actual = sp.groebner(eliminated, d, a, x, domain=field, order='lex')
    expected = sp.groebner([x**3, a-3*x*x, d-delta*(1+x*x)],
                           d, a, x, domain=field, order='lex')
    assert actual == expected
    assert actual.reduce(x*x)[1] == x*x
    assert actual.reduce(x**3)[1] == 0

    basis, J = sp4_basis()
    g = sp.diag(sp.I*z0, -sp.I*z0, 1/(sp.I*z0), -1/(sp.I*z0))
    g = g.applyfunc(simp)
    gi = g.inv().applyfunc(simp)
    H = sp.diag(1, 1, -1, -1)
    assert (g.T*J*g-J).applyfunc(simp) == sp.zeros(4)
    assert simp(sp.trace(g)) == 0
    lam = -delta/24
    theta = ((g*g-gi*gi)/12).applyfunc(simp)
    assert (theta-lam*H).applyfunc(simp) == sp.zeros(4)
    t1, t2 = g[0, 0], g[1, 1]
    root_values = [t1/t2, t2/t1, t1*t2, 1/(t1*t2),
                   t1*t1, t2*t2, t1**-2, t2**-2]
    assert all(simp(t-1) != 0 for t in root_values)

    # Construct Ad(g) and every ad(X_i) directly on sp4, without characters.
    Ad = sp.Matrix.hstack(*(coords(g*M*gi) for M in basis)).applyfunc(simp)
    ads = [sp.Matrix.hstack(*(coords(X*M-M*X) for M in basis)) for X in basis]
    gram = sp.Matrix(10, 10, lambda i, j: sp.trace(ads[i]*ads[j]))
    assert gram == sp.Matrix(10, 10, lambda i, j: 6*sp.trace(basis[i]*basis[j]))
    assert gram.det() != 0
    linear = sp.Matrix([simp(sp.trace(Ad*ad)) for ad in ads])
    assert (gram.inv()*linear-coords(theta)).applyfunc(simp) == sp.zeros(10, 1)

    trace_diff = sp.Matrix(10, 10,
        lambda j, i: simp(sp.trace(Ad*ads[i]*ads[j])))
    trace_diff = (gram.inv()*trace_diff).applyfunc(simp)
    matrix_columns = []
    for X in basis:
        derivative = (sp.trace(g*X)*(g-gi)
                      +g*X*g+g*g*X+X*gi*gi+gi*X*gi)/12
        matrix_columns.append(coords(derivative).applyfunc(simp))
    matrix_diff = sp.Matrix.hstack(*matrix_columns)
    assert (trace_diff-matrix_diff).applyfunc(simp) == sp.zeros(10)
    assert trace_diff.rank() == 7
    assert sp.Matrix.hstack(*(coords(g*X-X*g) for X in basis)).rank() == 8
    assert sp.Matrix.hstack(*(coords(theta*X-X*theta) for X in basis)).rank() == 6

    # Universal adjoint differential formula, computed independently from
    # the defining-matrix and trace differentials above. The Killing-form
    # Casimir is the identity on the adjoint module.
    inverse_gram = gram.inv()
    pairs = [(i,j,inverse_gram[i,j]) for i in range(10) for j in range(10)
             if inverse_gram[i,j] != 0]
    casimir = sum((v*ads[i]*ads[j] for i,j,v in pairs),sp.zeros(10))
    assert casimir == sp.eye(10)
    sandwich = sum((v*ads[i]*Ad*ads[j] for i,j,v in pairs),sp.zeros(10))
    assert (sandwich-trace_diff).applyfunc(simp) == sp.zeros(10)
    def bracket(P,Q): return P*Q-Q*P
    conjugation_casimir = sum((v*bracket(ads[i],bracket(ads[j],Ad))
                              for i,j,v in pairs),sp.zeros(10))
    assert (Ad-conjugation_casimir/2-trace_diff).applyfunc(simp) == sp.zeros(10)
    ad_theta = sp.Matrix.hstack(*(coords(theta*X-X*theta) for X in basis))
    assert (trace_diff*(sp.eye(10)-Ad.inv())-ad_theta).applyfunc(simp) == sp.zeros(10)
    assert (trace_diff-gram.inv()*trace_diff.T*gram-ad_theta).applyfunc(simp) == sp.zeros(10)
    assert (trace_diff*Ad-Ad*trace_diff).applyfunc(simp) == sp.zeros(10)

    # A kernel is not automatically closed under the Lie bracket. Use a
    # different point of the SAME Sp4/GL2 family, without a new Levi case.
    g2 = sp.diag(2*sp.I,-2*sp.I,-sp.I/2,sp.I/2)
    g2inv = g2.inv()
    theta2 = (g2*g2-g2inv*g2inv)/12
    assert theta2 == -sp.Rational(5,16)*H
    def differential2(X):
        return (sp.trace(g2*X)*(g2-g2inv)+g2*X*g2+g2*g2*X
                +X*g2inv*g2inv+g2inv*X*g2inv)/12
    E,F = basis[1],basis[2]
    coroot = E*F-F*E
    assert differential2(E) == sp.zeros(4)
    assert differential2(F) == sp.zeros(4)
    assert differential2(coroot) == -sp.Rational(7,4)*coroot
    Ad2 = sp.Matrix.hstack(*(coords(g2*M*g2inv) for M in basis))
    sandwich2 = sum((v*ads[i]*Ad2*ads[j] for i,j,v in pairs),sp.zeros(10))
    assert sandwich2*coords(E) == sp.zeros(10,1)
    assert sandwich2*coords(F) == sp.zeros(10,1)
    assert sandwich2*coords(coroot) == -sp.Rational(7,4)*coords(coroot)
    assert sp.Matrix.hstack(*(coords(g2*X-X*g2) for X in basis)).rank() == 8
    assert sp.Matrix.hstack(*(coords(theta2*X-X*theta2) for X in basis)).rank() == 6

    hessian = sp.Matrix([[sp.diff(phi, x, 2), D(sp.diff(phi, x))],
                         [D(sp.diff(phi, x)), D(D(phi))]])
    assert hessian.subs(at).applyfunc(simp) == sp.diag(0, 4)
    assert simp(6*sp.trace(theta*theta)) == -sp.Rational(1, 2)

    report = dict(
        result='passed', group='Sp4(C)', representation='adjoint',
        proposed_stronger_lemma_false=True, generic_conjecture_disproved=False,
        all_simple_adjoint_conjecture_proved=False, new_Levi_coverage=False,
        character_from_defining_Sym2=True, local_scheme_exact=True,
        transverse_local_algebra='C[[x]]/(x^3)', transverse_length=3,
        square_class_nonzero=True, incidence_local_equation='a*x=0',
        central_target='h=lambda*diag(1,1,-1,-1), lambda^2=-1/48',
        target_killing_norm='-1/2', target_centralizer='GL2',
        source_regular_semisimple=True, source_centralizer_dimension=2,
        target_centralizer_dimension=4, adjoint_differential_rank=7,
        quotient_hessian=[[0, 0], [0, 4]],
        trace_form_pairings=100, trace_definition_pairings=10,
        differential_pairings=100,
        adjoint_Casimir_sandwich_formula_verified=True,
        conjugation_Casimir_formula_verified=True,
        full_equivariance_operator_identity_verified=True,
        kernel_Lie_subalgebra_shortcut_false=True,
        kernel_control_target='-5/16*diag(1,1,-1,-1)',
        kernel_control_coroot_eigenvalue='-7/4',
        full_completed_local_ring='C[[u1,u2,x]]/(x^3)',
        higher_order=higher_order_audit(),
        generic_target=generic_target_audit(),
        reducedness=reducedness_audit(),
        proof_note='UNIFORM_ADJOINT_PROOF_TARGET.md')
    paths = [Path(__file__), ROOT/report['proof_note']]
    report['source_sha256'] = {
        str(p.relative_to(ROOT)).replace('\\', '/'):
        sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
        for p in paths}
    (ROOT/'uniform-adjoint-obstruction-audit.json').write_text(
        json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
