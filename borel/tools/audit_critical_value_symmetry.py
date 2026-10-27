"""Exact controls for the critical-value symmetry criterion (native SymPy).

Computations corroborate the written geometric proof; they do not formally
verify it. No Sage or existing mathematical audit helper is imported.
"""
from pathlib import Path
from itertools import combinations, product
from collections import Counter
from hashlib import sha256
import json
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def log_hessian(f, z):
    return sp.Matrix([[sp.cancel(v*sp.diff(u*sp.diff(f, u), v)) for v in z] for u in z])


def support_rank_rational(f, z):
    numerator, denominator = sp.cancel(f).as_numer_denom()
    exponents = set(sp.Poly(numerator, *z).monoms()) | set(sp.Poly(denominator, *z).monoms())
    first = min(exponents)
    return sp.Matrix([[a-b for a, b in zip(e, first)] for e in exponents]).rank()


def roots_from_cartan(kind, n):
    A = sp.eye(n)*2
    if kind == 'A':
        edges = [(i, i+1) for i in range(n-1)]
    elif kind == 'D':
        edges = [(i, i+1) for i in range(n-2)] + [(n-3, n-1)]
    else:
        edges = [(0, 2), (2, 3), (1, 3)] + [(i, i+1) for i in range(3, n-1)]
    for i, j in edges:
        A[i, j] = A[j, i] = -1
    roots = {tuple(int(i == j) for i in range(n)) for j in range(n)}
    todo = list(roots)
    while todo:
        root = todo.pop()
        for j in range(n):
            reflected = list(root)
            reflected[j] -= sum(root[i]*int(A[i, j]) for i in range(n))
            reflected = tuple(reflected)
            if reflected not in roots:
                roots.add(reflected)
                todo.append(reflected)
    return A, sorted(roots)


def direct_adjoint_matrix_controls():
    """Check trace projection itself, without exterior/symmetric trace formulas."""
    records = []
    u = sp.Symbol('u')
    for kind, size, tau, expected_character, scalar in (
        ('B4',9,-1,-1,sp.Rational(3,14)),
        ('C4',8,0,0,sp.Rational(-3,20)),
        ('D4',8,0,0,sp.Rational(1,4))):
        symplectic = kind.startswith('C')
        middle = int(size == 9)
        block = sp.Matrix([[0,0,0,1],[1,0,0,-tau],[0,1,0,-3],[0,0,1,tau]])
        inverse = block.inv()
        g = sp.diag(block, *([sp.ones(1)] if middle else []), inverse.T)
        form = sp.zeros(size)
        for i in range(4):
            form[i,size-4+i] = 1
            form[size-4+i,i] = -1 if symplectic else 1
        if middle:
            form[4,4] = 1
        assert g.T*form*g == form and g.det() == 1
        p = sp.Poly(u**4-tau*u**3+3*u**2+tau*u-1,u)
        reciprocal = sp.Poly(u**4-tau*u**3-3*u**2+tau*u-1,u)
        assert sp.gcd(p,p.diff()).degree() == sp.gcd(p,reciprocal).degree() == 0
        assert p.eval(1) != 0 and p.eval(-1) != 0
        pairs = [(i,j) for i in range(size) for j in range(i if symplectic else i+1,size)]
        form_inverse = form.inv()
        basis = []
        for i,j in pairs:
            S = sp.zeros(size)
            S[i,j] = 1
            S[j,i] = 1 if symplectic else -1
            basis.append(form_inverse*S)
        def coordinates(X):
            S = form*X
            return sp.Matrix([S[i,j] for i,j in pairs])
        g_inverse = g.inv()
        Ad = sp.Matrix.hstack(*(coordinates(g*X*g_inverse) for X in basis))
        assert sp.trace(Ad) == expected_character
        target = scalar*sp.diag(*([1]*4+[0]*middle+[-1]*4))
        ad_target = sp.Matrix.hstack(*(coordinates(target*X-X*target) for X in basis))
        difference = Ad-ad_target
        # tr(Ad(g) ad X) = tr(ad(target) ad X) for every Lie basis X.
        for X in basis:
            pairing = sum((difference.row(j)*coordinates(X*Y-Y*X))[0]
                          for j,Y in enumerate(basis))
            assert pairing == 0
        records.append(dict(type=kind, adjoint_dimension=len(basis),
            adjoint_character=expected_character, target_scalar=str(scalar),
            all_trace_projection_pairings_checked=True,
            regular_semisimple_matrix_source=True))
    return records


def main():
    z, w, v = sp.symbols('z w v', nonzero=True)
    rational_cases = [
        (z+w, (z, w)), (z/w+w/z, (z, w)), ((z+w)/(z-w), (z, w)),
        (1/(z+w), (z, w)), ((z+w+1)/(z*w), (z, w)),
        ((z*z-w*w)/(z-w), (z, w)), ((z+w)/(z+w), (z, w)),
        (z/w+v, (z, w, v)), ((z+w)/(z-w)+v, (z, w, v)),
        ((z+w+v)/(z*w*v+1), (z, w, v)),
        (z+w+1/z-1/w, (z, w)),
    ]
    rank_records = []
    for f, variables in rational_cases:
        expected = support_rank_rational(f, variables)
        rank = log_hessian(f, variables).rank(iszerofunc=lambda x: sp.cancel(x) == 0)
        assert rank == expected, f
        rank_records.append(dict(function=str(f), variables=list(map(str, variables)), rank=rank))

    # A nonlinear coefficient locus with constant critical value still has
    # a dominant conormal target map. The inverse is explicit on c1+c2 !=0.
    c1, c2, u = sp.symbols('c1 c2 u')
    inverse = {u: c1+c2, z: c1/(c1+c2), w: c2/(c1+c2)}
    assert sp.cancel((z+w-1).subs(inverse, simultaneous=True)) == 0
    assert sp.cancel((u*z-c1).subs(inverse, simultaneous=True)) == 0
    assert sp.cancel((u*w-c2).subs(inverse, simultaneous=True)) == 0
    # Positive-dimensional subtorus stabilizer: z*w=1 and mu=(u,u).
    assert sp.expand((u*z*w-u).subs(w, 1/z)) == 0

    # The clean rank-drop component supplies the generic target fiber.
    x = sp.symbols('x')
    phi = (z-1)*x*x/2
    assert sp.diff(phi, x) == x*(z-1)
    assert sp.expand((z*sp.diff(phi, z)).subs(z, 1)) == x*x/2

    # E6 contact: do not remove z=+/-i, the entire generic target source.
    a, b, f1, f2, c = sp.symbols('a b f1 f2 c')
    phi_contact = a*b+z*z+z**-2+(z+z**-1)*u
    for sign in (-1, 1):
        point = {a: 0, b: 0, z: sign*sp.I, u: c/(2*sign*sp.I)}
        assert all(sp.simplify(sp.diff(phi_contact, t).subs(point)) == 0 for t in (a, b, u, f1, f2))
        assert sp.simplify((z*sp.diff(phi_contact, z)-c).subs(point)) == 0
        assert sp.simplify(phi_contact.subs(point)) == -2

    # Actual SL4 adjoint control: a reduced-support-only criterion is false.
    y = sp.symbols('y')
    phi_sl4 = x*x+y*y+(z*z+z**-2)*x*y-1
    U, V = sp.symbols('U V')
    changed = sp.expand(phi_sl4.subs({x:(U+V)/2, y:(U-V)/2}, simultaneous=True))
    assert sp.cancel(changed-((z+1/z)**2*U*U/4-(z-1/z)**2*V*V/4-1)) == 0
    relative_v = sp.factor(sp.diff(changed, V))
    assert sp.cancel(relative_v/((z-1)**2*V)+(z+1)**2/(2*z*z)) == 0
    assert sp.simplify((z*sp.diff(changed, z)).subs({U:0, z:1})) == 0
    assert sp.simplify(changed.subs({U:0, z:1})) == -1

    # Componentwise hypotheses are essential: the union of two critical
    # value components has finite stabilizer while both images are proper.
    w = sp.symbols('w', nonzero=True)
    reducible = u*(z-1)*(w-1)
    mu1, mu2 = z*sp.diff(reducible,z), w*sp.diff(reducible,w)
    assert sp.cancel(mu1*mu2/(u*u*z*w)-(z-1)*(w-1)) == 0
    # A genuinely algebraic critical value: elimination must retain both
    # conjugate branches, rather than choosing a square root globally.
    t = sp.symbols('t')
    algebraic_phi = x**3/3-(z+w-1)*x
    critical_value_equation = sp.factor(sp.resultant(sp.diff(algebraic_phi,x), t-algebraic_phi, x))
    assert sp.expand(critical_value_equation-(t*t-4*(z+w-1)**3/9)) == 0

    # Character composition: an exact original-ideal certificate for
    # erased E6 contact components, valid before passing to a radical.
    a, b, u, c = sp.symbols('a b u c')
    F = z+1/z
    q = a*b+F*u+F*F
    Dq = z*sp.diff(q,z)
    qq, dd = sp.symbols('qq dd')
    composition_powers = []
    for m in range(2,7):
        composed = q**m
        combination = a*sp.diff(composed,a)+(u+F)*sp.diff(composed,u)
        assert sp.expand(combination-m*q**m) == 0
        new_mu = m*qq**(m-1)*dd
        telescoping = (c-new_mu)*sum(c**(m-1-j)*new_mu**j for j in range(m))
        certificate = telescoping+m**(m-1)*qq**(m*(m-2))*dd**m*(m*qq**m)
        # Compare as a polynomial in q,Dq first to avoid unnecessary large
        # expansion in central Laurent coordinates.
        assert sp.factor(certificate-c**m) == 0
        composition_powers.append(m)
    # The nonconstant SL3 adjoint critical value survives every tested
    # polynomial composition; this controls the other side of the rule.
    S, R = z**3+z**-3, z**3-z**-3
    sl3_value = -S*S/4
    for m in range(1,6):
        assert sp.cancel(z*sp.diff(sl3_value**m,z)-m*sl3_value**(m-1)*(-3*S*R/2)) == 0

    # Direct root sums, independent of the Killing-form splitting proof,
    # check the closed type-A involution Hessian determinant.
    involution_records = []
    for N in range(2,13):
        for p in range(N+1):
            q = N-p
            signs = [1]*p+[-1]*q
            H = sp.zeros(N-1)
            for i in range(N):
                for j in range(N):
                    if i == j:
                        continue
                    root = sp.Matrix([int(i==a)-int(j==a)-int(i==N-1)+int(j==N-1)
                                      for a in range(N-1)])
                    H += signs[i]*signs[j]*root*root.T
            expected = (-1)**q*2**(N-1)*N*N*(p-q)**(N-2)
            assert H.det() == expected
            assert H.rank() == (1 if p==q and N>2 else N-1)
        involution_records.append(dict(N=N, sign_multiplicities_checked=N+1))

    # Independent finite controls of the simply-laced root proof: the
    # projected external weights span alpha-perp and their doubled norms
    # exceed every orthogonal-root norm. The argument is not inferred from
    # the checked rank range.
    root_records = []
    for kind, n in [('A',n) for n in range(2,9)]+[('D',n) for n in range(4,9)]+[('E',6),('E',7),('E',8)]:
        cartan, roots = roots_from_cartan(kind, n)
        for j in range(n):
            projected = []
            for root in roots:
                pairing = sum(root[i]*int(cartan[i,j]) for i in range(n))
                if pairing == 1:
                    twice = sp.Matrix([2*root[i]-int(i==j) for i in range(n)])
                    assert (twice.T*cartan*twice)[0] == 6
                    projected.append(list(twice))
            assert sp.Matrix(projected).rank() == n-1
        root_records.append(dict(type=f'{kind}{n}', roots=len(roots), simple_roots_checked=n))

    # All adjoint-torus involutions, not a Levi enumeration. The general
    # proof uses the dual-Coxeter obstruction; this independently checks
    # the exceptional endpoints and the genuine E6 failure.
    involution_root_records = []
    for n in (6,7,8):
        cartan, roots = roots_from_cartan('E',n)
        labels = [sp.Matrix([sum(root[i]*int(cartan[i,j]) for i in range(n))
                             for j in range(n)]) for root in roots]
        rank_counts, centralizer_counts = Counter(), Counter()
        for bits in product((0,1), repeat=n):
            H = sp.zeros(n)
            centralizer_roots = 0
            for root, label in zip(roots, labels):
                sign = (-1)**sum(a*b for a,b in zip(root,bits))
                H += sign*label*label.T
                centralizer_roots += int(sign == 1)
            rank_counts[H.rank()] += 1
            centralizer_counts[centralizer_roots] += 1
        if n in (7,8):
            assert rank_counts == {n:2**n}
        else:
            assert rank_counts == {6:28, 1:36}
        involution_root_records.append(dict(type=f'E{n}', characters_checked=2**n,
            hessian_rank_counts=dict(rank_counts), centralizer_root_counts=dict(centralizer_counts)))

    # E8: verify the compact witness, then separately use the uneliminated
    # archived full character and freshly enumerated 240 roots.
    x, y, w, a, A, B, C = sp.symbols('x y w a A B C')
    p, q, r = x*y+w*a, x*w+y*a, x*a+y*w
    residual = x*x+y*y+w*w+a*a-p*p-q*q-r*r+x*y*w*a-A*q*r-B*p*r-C*p*q-B*C*p-A*C*q-A*B*r-A*B*C-4
    point = {x:0, y:0, w:0, a:0, A:2, B:2, C:2}
    hh = sp.hessian(residual, (x,y,w,a)).subs(point)
    derived = hh[:3,:3]
    normal = (hh[3,3]-(hh[3,:3]*derived.inv()*hh[:3,3])[0])
    assert derived.det() == -216 and normal == 10
    uu, vv, ss = sp.symbols('uu vv ss', nonzero=True)
    traces = {A:uu+1/uu, B:vv+1/vv, C:uu*vv+1/(uu*vv), a:ss+1/ss}
    pulled = residual.subs(traces, simultaneous=True)
    all_variables = (x,y,w,ss,uu,vv)
    operators = [lambda f,t=t: sp.diff(f,t) if t in (x,y,w) else t*sp.diff(f,t) for t in all_variables]
    at = {x:0,y:0,w:0,ss:sp.I,uu:1,vv:1}
    assert all(sp.simplify(op(pulled).subs(at)) == 0 for op in operators)
    hessian = sp.Matrix([[sp.simplify(op1(op2(pulled)).subs(at)) for op2 in operators] for op1 in operators])
    schur = hessian[3:,3:]-hessian[3:,:3]*hessian[:3,:3].inv()*hessian[:3,3:]
    assert schur == sp.diag(-40, sp.Matrix([[-16,-8],[-8,-16]]))
    assert schur.det() == -7680

    row = next(r for r in json.loads((ROOT/'e8-final-character-models.json').read_text()) if r['seed']==[1,2,3,5,7])
    variables = sp.symbols(' '.join(row['variables']))
    full = sp.Add(*(sp.Rational(t['coefficient'])*sp.prod(v**e for v,e in zip(variables,t['exponent']+t['charge'])) for t in row['terms']))
    full = full.subs(variables[5], ss+1/ss)
    full_variables = variables[:5]+(ss,)+variables[6:]
    full_point = dict(zip(full_variables, (-3,0,-3,0,0,sp.I,1,1)))
    full_ops = [lambda f,t=t,i=i: sp.diff(f,t) if i<5 else t*sp.diff(f,t) for i,t in enumerate(full_variables)]
    assert all(sp.simplify(op(full).subs(full_point)) == 0 for op in full_ops)
    full_hessian = sp.Matrix([[sp.simplify(op1(op2(full)).subs(full_point)) for op2 in full_ops] for op1 in full_ops])
    full_schur = full_hessian[5:,5:]-full_hessian[5:,:5]*full_hessian[:5,:5].inv()*full_hessian[:5,5:]
    assert full_hessian[:5,:5].det() == 216
    assert full_schur.det() == -69120 and full_hessian.det() == -14929920

    cartan, roots = roots_from_cartan('E',8)
    columns = row['columns']
    raw = Counter()
    for root in roots:
        labels = [sum(root[i]*int(cartan[i,j]) for i in range(8)) for j in range(8)]
        raw[tuple(sum(labels[j]*column[j] for j in range(8)) for column in columns)] += 1
    raw[(0,)*8] += 8
    fundamental = {}
    for component in row['components']:
        letters = []
        for k in range(len(component)+1):
            weight = [0]*8
            for j, index in enumerate(component):
                weight[index] = int(k==j)-int(k==j+1)
            letters.append(weight)
        for j, index in enumerate(component):
            fundamental[index] = Counter(tuple(sum(letters[k][v] for k in subset) for v in range(8))
                                         for subset in combinations(range(len(component)+1),j+1))
    rebuilt = Counter()
    for term in row['terms']:
        expanded = Counter({(0,)*6+tuple(term['charge']): term['coefficient']})
        for i, power in enumerate(term['exponent']):
            for _ in range(power):
                updated = Counter()
                for left, lc in expanded.items():
                    for right, rc in fundamental[i].items():
                        updated[tuple(a+b for a,b in zip(left,right))] += lc*rc
                expanded = updated
        rebuilt.update(expanded)
    assert {e:c for e,c in rebuilt.items() if c} == dict(raw)

    # All-rank type A argument: check identities with independent formal
    # trace sums, then their one-parameter restriction. Integer rank bounds
    # are not used for the written proof.
    X, AA, BB, zz, zi, zj, UU, VV = sp.symbols('X AA BB zz zi zj UU VV', nonzero=True)
    type_a_phi = (zz*X+AA)*(X/zz+BB)-1
    type_a_critical = -(zz*BB+AA/zz)/2
    type_a_value = -1-(zz*BB-AA/zz)**2/4
    assert sp.cancel(sp.diff(type_a_phi,X).subs(X,type_a_critical)) == 0
    assert sp.cancel(type_a_phi.subs(X,type_a_critical)-type_a_value) == 0
    assert sp.expand(zi*(zi*VV+UU/zi)-zj*(zj*VV+UU/zj)-(zi**2-zj**2)*VV) == 0
    uu, qrank, singles = sp.symbols('uu qrank singles', positive=True)
    qrank = sp.Symbol('qrank', integer=True, positive=True)
    singles = sp.Symbol('singles', integer=True, positive=True)
    character_restriction = type_a_value.subs({zz:uu**singles, AA:singles*uu**(-2*qrank), BB:singles*uu**(2*qrank)})
    expected_restriction = -1-singles**2*(uu**(2*qrank+singles)-uu**(-2*qrank-singles))**2/4
    assert sp.simplify(character_restriction-expected_restriction) == 0

    # Independent companion-matrix route for the uniform B/C/D four-block
    # theorem, including repeated-root parameters (no discriminant inverse).
    tau, LL, TT, number = sp.symbols('tau LL TT number')
    companion = sp.Matrix([[0,0,0,1],[1,0,0,-tau],[0,1,0,-LL],[0,0,1,tau]])
    inverse = companion.inv()
    assert sp.trace(companion) == sp.trace(inverse) == tau
    assert sp.expand(sp.trace(companion**2)) == tau**2-2*LL
    assert sp.expand(sp.trace(inverse**2)) == tau**2+2*LL
    orthogonal_numerator = tau*(companion-inverse)-(companion**2-inverse**2)
    symplectic_numerator = -tau*(companion-inverse)+(companion**2-inverse**2)
    assert sp.simplify(orthogonal_numerator-LL*sp.eye(4)) == sp.zeros(4)
    assert sp.simplify(symplectic_numerator+LL*sp.eye(4)) == sp.zeros(4)
    T_B = 1/(1-2*number)
    b_B = sp.cancel((TT**2-(2*number*TT**2+1)).subs(TT,T_B)/2)
    assert sp.cancel(b_B+number/(2*number-1)) == 0
    assert sp.cancel((2*number-1)*b_B+number) == 0
    assert ((TT**2-2*number*TT**2)/2).subs(TT,0) == 0
    assert ((TT**2+2*number*TT**2)/2).subs(TT,0) == 0
    direct_matrix_records = direct_adjoint_matrix_controls()

    report = dict(result='passed', rational_support_vs_log_hessian=rank_records,
        nonlinear_constraint_conormal_inverse=True, rank_drop_component_retained=True,
        E6_contact_excess=2, SL4_nonclean_point_control=True,
        componentwise_stabilizer_control=True,
        algebraic_critical_value_equation='t^2-4*(z+w-1)^3/9',
        E6_composition_original_ideal_target_power_certificates=composition_powers,
        SL3_nonconstant_critical_value_composition_control=True,
        type_A_composition_completed_square_and_unit_elimination=True,
        type_A_composition_all_rank_nonconstant_restriction=True,
        BCD_four_block_companion_identity=True,
        BCD_four_block_constant_values=['-q/(2*q-1)','0','0'],
        BCD_direct_adjoint_matrix_controls=direct_matrix_records,
        type_A_involution_hessian_checks=involution_records,
        exceptional_involution_root_checks=involution_root_records,
        simply_laced_root_checks=root_records,
        E8_full_roots_rebuilt=len(roots), E8_original_derived_hessian_determinant=216,
        E8_original_central_schur_determinant=-69120,
        E8_original_full_hessian_determinant=-14929920,
        E8_residual_central_schur_determinant=-7680,
        E8_last_generic_nonemptiness_witness=True,
        scope='Independent exact controls and full-root character witness; the general stabilizer criterion is a written proof, not mechanically verified.',
        literature_priority_claim=False)
    names = ['tools/audit_critical_value_symmetry.py','ADJOINT_NONEMPTINESS_CRITICAL_VALUES.md',
             'CHARACTER_COMPOSITION_FIBERS.md','e8-final-character-models.json']
    report['source_sha256'] = {name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in names}
    (ROOT/'critical-value-symmetry-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__ == '__main__':
    main()
