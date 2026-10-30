"""Direct root audit of the uniform nonlinear-coefficient theorems.

The proof is in UNIFORM_ADJOINT_PROOF_TARGET.md. Reconstructs cubic
constituents from highest roots, independently of stored character
polynomial inversion; checks no Springer-fiber smoothness claim.
"""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import gzip
import json

import sympy as sp
from audit_affine_branching_certificates import cartan, root_data

ROOT = Path(__file__).resolve().parents[1]


def nonlinear_coefficients(grouped, A, I, H, gram_inverse, roots, counts):
    """Check the simplex/cross-polytope proof on actual highest roots."""
    inverse = A.inv()
    central_roots = [sp.Matrix(b) for b in roots if not any(b[i] for i in I)]
    positive = {tuple(b) for b in central_roots if all(v >= 0 for v in inverse*b)}
    central_simple = [sp.Matrix(b) for b in sorted(positive)
                      if not any(tuple(x-y for x,y in zip(b,a)) in positive for a in positive)]
    coordinate_monomials = set()
    for exponent, charges in grouped.items():
        gammas = [A*H*gram_inverse*q for q in charges]
        norm = (charges[0].T*gram_inverse*charges[0])[0]
        m = 2-norm
        assert 1 <= m < 2
        assert all((q.T*gram_inverse*q)[0] == norm for q in charges)
        if m > 1:
            p = len(gammas)
            assert p*(m-1) <= 1
            factor = {tuple(a-b) for a in gammas for b in gammas if a != b}
            assert len(factor) == p*(p-1) and factor <= roots
            for a,b in combinations(gammas,2):
                assert (a.T*inverse*b)[0] == 1-m
            barycenter = sum(gammas,sp.zeros(A.rows,1))/p
            assert (barycenter.T*inverse*barycenter)[0] == 1-m+sp.Rational(1,p)
            counts[f'degree_{sum(exponent)}_simplex_{p}_weights'] += 1
        else:
            assert sum(exponent) == 2
            weights = {tuple(a) for a in gammas}
            assert {tuple(-a) for a in gammas} == weights
            axes = [a for a in gammas if tuple(a) > tuple(-a)]
            p = len(axes)
            assert len(gammas) == 2*p
            assert all((a.T*inverse*a)[0] == 1 for a in axes)
            assert all((a.T*inverse*b)[0] == 0 for a,b in combinations(axes,2))
            factor = {tuple(e*a+f*b) for a,b in combinations(axes,2)
                      for e,f in product((-1,1),repeat=2)}
            assert len(factor) == 2*p*(p-1) and factor <= roots
            counts[f'degree_2_cross_polytope_{p}_axes'] += 1
        # Verify that these are complete isolated factors of the centralizer,
        # not only root subsystems of the indicated abstract type.
        weights = {tuple(a) for a in gammas}
        for beta in central_roots:
            coroot = inverse*beta
            for gamma in gammas:
                pairing = (gamma.T*coroot)[0]
                assert pairing in (-1,0,1)
                assert tuple(gamma-pairing*beta) in weights
                if tuple(beta) not in factor:
                    assert pairing == 0
        # A single fixed simple system for H, common to ALL coefficients.
        highest = [g for g in gammas if all(tuple(g+b) not in weights for b in central_simple)]
        labels = [tuple((g.T*inverse*b)[0] for b in central_simple) for g in highest]
        assert all(v in (0,1) for label in labels for v in label)
        if factor:
            assert len(highest) == 1
            support = frozenset(i for i,v in enumerate(labels[0]) if v)
            assert len(support) == (2 if m == 1 and p == 2 else 1)
            if len(support) == 2:
                for i in support:
                    assert all((central_simple[i].T*inverse*b)[0] == 0
                               for j,b in enumerate(central_simple) if i != j)
            coordinate_monomials.add(support)
            counts[f'common_frame_{len(support)}_fundamental_factors'] += 1
        else:
            assert len(highest) == (1 if m > 1 else 2)
            assert all(not any(label) for label in labels)
            counts['common_frame_unit' if m > 1 else 'common_frame_central_binomial'] += 1
        counts['nonlinear_coefficient_weight_sets'] += 1
    # Independent coordinate-ideal check: reconstruct a squarefree monomial
    # ideal from its minimal coordinate-prime decomposition (blocker duality).
    vertices = set().union(*coordinate_monomials) if coordinate_monomials else set()
    def minimal(sets):
        return {a for a in sets if not any(b < a for b in sets)}
    def blocker(edges):
        return minimal({frozenset(c) for k in range(len(vertices)+1)
                        for c in combinations(sorted(vertices),k)
                        if all(set(c)&edge for edge in edges)})
    if coordinate_monomials:
        assert blocker(blocker(coordinate_monomials)) == minimal(coordinate_monomials)
    counts['common_centralizer_coordinate_frames_checked'] += 1


def check_coefficient_boundaries():
    """Exact derivatives and all sign strata through the ambient rank bound.

    The proof in the note works for arbitrary p. This finite check verifies
    its formulas and the D2 crossing independently of character inversion.
    """
    a,b = sp.symbols('a b',nonzero=True)
    trace = a+1/a+b+1/b
    assert sp.cancel(trace-(a+b)*(a*b+1)/(a*b)) == 0
    crossing_jacobian = sp.Matrix([a+b,a*b+1]).jacobian([a,b])
    assert crossing_jacobian.det().subs({a:1,b:-1}) == 2
    assert crossing_jacobian.det().subs({a:-1,b:1}) == -2
    for p in range(1,9):
        q = sp.symbols(f'q0:{p}',nonzero=True)
        standard = sum(q)
        for i in range(1,p):
            assert q[i]*sp.diff(standard,q[i])-q[0]*sp.diff(standard,q[0]) == q[i]-q[0]
        assert standard.subs({v:q[0] for v in q}) == p*q[0]
        f = sum(v+1/v for v in q)
        grad = [sp.expand(v*sp.diff(f,v)) for v in q]
        assert grad == [v-1/v for v in q]
        for signs in product((-1,1),repeat=p):
            if sum(signs):
                continue
            point = dict(zip(q,signs))
            assert f.subs(point) == 0
            hessian = sp.Matrix([[w*sp.diff(grad[i],w) for w in q] for i in range(p)])
            assert hessian.subs(point) == sp.diag(*[2*v for v in signs])
            if p >= 4:
                i,j = next((i,j) for i,j in combinations(range(p),2) if signs[i] == signs[j])
                assert (q[i]/q[j]-1).subs(point) == 0
            else:
                assert p == 2
                assert (q[0]/q[1]).subs(point) == -1
                assert (q[0]*q[1]).subs(point) == -1


def check_linear_resonance(rows):
    """Test the failed inference from coefficient boundaries to smoothness.

    Only the linear coefficient is changed in an archived adjoint model.
    This is an honest representation of its reductive Levi, NOT an
    adjoint representation of an ambient simple group.
    """
    z = sp.Symbol('z',nonzero=True)
    a,b,u,v,x = sp.symbols('a b u v x')
    S = z**3+z**-3
    left = z**-2+z*x
    right = z**2+z**-1*x
    original = a*b+u*v+x*x-2+S*x+left*a*v+right*b*u
    row = next(r['character'] for r in rows if r['character']['type']=='E6'
               and r['character']['nodes']==[1,2,3,5,6])
    variables = [a,x,b,u,v]
    archived = sum(t['coefficient']*z**t['charge'][0]*sp.prod(
        q**e for q,e in zip(variables,t['exponent'])) for t in row['terms'])
    assert sp.expand(original-archived) == 0
    changed = original+S*x
    difference = sp.Poly(sp.expand(changed-original),*variables)
    assert difference.total_degree() == 1
    assert sp.expand(changed.subs({a:3,b:3,u:3,v:3,x:2,z:1})) == 82
    def D(f): return z*sp.diff(f,z)
    assert D(D(original)).subs({a:3,b:3,u:3,v:3,x:2,z:1}) == 144
    assert D(D(changed)).subs({a:3,b:3,u:3,v:3,x:2,z:1}) == 180

    # Root-space reconstruction of the four-dimensional added module and
    # its trace form, independent of character differentiation.
    C = sp.Matrix(cartan('E6'))
    positive = root_data(tuple(map(tuple,C.tolist())))[1]
    roots = set(positive) | {tuple(-q for q in r) for r in positive}
    H = sp.Matrix(row['central_cocharacters'][0])
    retained = [i-1 for i in row['nodes']]
    added = []
    added_weights = set()
    root_norm = 0
    for alpha in roots:
        charge = (sp.Matrix(alpha).T*H)[0]
        root_norm += charge**2
        restricted = tuple(alpha[i] for i in retained)
        if restricted in ((0,1,0,0,0),(0,-1,0,0,0)) and abs(charge)==3:
            added.append(charge)
            added_weights.add(alpha)
    assert len(added) == 4 and sum(q*q for q in added) == 36
    assert root_norm == 144

    # The flat torus Laplacian is diagonal on weight monomials. Both
    # characters obey Delta(Phi)=2*(Phi-6), despite the altered critical
    # scheme being nonreduced. This is NOT the radial group Casimir.
    inverse = C.inv()
    assert all((sp.Matrix(alpha).T*inverse*sp.Matrix(alpha))[0] == 2
               for alpha in roots)
    altered_multiplicities = {alpha:1+int(alpha in added_weights) for alpha in roots}
    assert sum(altered_multiplicities.values())+6 == 82
    def reflected(alpha, i):
        return tuple(alpha[j]-alpha[i]*C[j,i] for j in range(6))
    assert all(reflected(alpha,i) in roots for alpha in roots for i in range(6))
    assert all(altered_multiplicities[reflected(alpha,i)] == altered_multiplicities[alpha]
               for alpha in roots for i in retained)
    witness = next((alpha,i,reflected(alpha,i)) for alpha in sorted(roots)
                   for i in range(6)
                   if altered_multiplicities[reflected(alpha,i)] != altered_multiplicities[alpha])

    # Torus invariant derivations and ordinary quotient derivatives have
    # different ideals at the discriminant, already for SL2 adjoint.
    q = sp.Symbol('q')
    t = sp.Symbol('t', nonzero=True)
    trace = t+1/t
    adjoint = q*q-1
    quotient_derivative = sp.diff(adjoint,q)
    invariant_derivative = (q*q-4)*quotient_derivative
    assert quotient_derivative.subs(q,2) == 4
    assert invariant_derivative.subs(q,2) == 0
    assert sp.simplify((t-1/t)*t*sp.diff(adjoint.subs(q,trace),t)
                       -invariant_derivative.subs(q,trace)) == 0
    unipotent = sp.Matrix([[1,1],[0,1]])
    lowering = sp.Matrix([[0,0],[1,0]])
    assert 2*sp.trace(unipotent)*sp.trace(unipotent*lowering) == 4

    # Exact completion, with no specialization of z.
    ap,bp,X = sp.symbols('ap bp X')
    completed = changed.subs({a:ap-right*u,b:bp-left*v}, simultaneous=True)
    expected = ap*bp+x*x+2*S*x-x*(x+S)*u*v-2
    assert sp.expand(completed-expected) == 0
    translated = sp.expand(expected.subs(x,X-S)-ap*bp)
    assert sp.expand(translated-(-S*S-2+X*X*(1-u*v)+S*X*u*v)) == 0

    # Localize at X-S (a unit on the branch X=u=v=0, S!=0),
    # and check BOTH ideal inclusions by elimination over Q(S).
    s,w = sp.symbols('s w')
    f = -s*s-2+X*X*(1-u*v)+s*X*u*v
    derivatives = [sp.diff(f,q) for q in (X,u,v)]
    field = sp.QQ.frac_field(s)
    localized = sp.groebner(derivatives+[w*(X-s)-1],w,X,u,v,domain=field)
    eliminated = [p.as_expr() for p in localized.polys if not p.as_expr().has(w)]
    claimed = sp.groebner([2*X+s*u*v,u*u*v,u*v*v],X,u,v,domain=field)
    actual = sp.groebner(eliminated,X,u,v,domain=field)
    assert actual == claimed
    assert claimed.reduce(u*v)[1] == u*v
    assert claimed.reduce((u*v)**2)[1] == 0
    assert claimed.reduce(sp.diff(f,s)+2*s)[1] == 0

    # Replay all original affine and central derivatives after the inverse
    # substitution in the nonreduced quotient, not just at its origin.
    nilpotent = sp.groebner([u*u*v,u*v*v],u,v,domain=sp.QQ.frac_field(z))
    value_x = -S-S*u*v/2
    substitution = {x:value_x, a:-right.subs(x,value_x)*u,
                    b:-left.subs(x,value_x)*v}
    for q in variables:
        residual = sp.expand(sp.diff(changed,q).subs(substitution,simultaneous=True))
        assert nilpotent.reduce(residual)[1] == 0
    mu = sp.expand(D(changed).subs(substitution,simultaneous=True))
    assert nilpotent.reduce(mu+2*S*D(S))[1] == 0
    assert sp.expand(-2*S*D(S)+6*(z**6-z**-6)) == 0
    assert sp.expand(D(-2*S*D(S))+36*(z**6+z**-6)) == 0
    c,t = sp.symbols('c t')
    assert sp.discriminant(6*t*t+c*t-6,t) == c*c+144

    # Exact source of the resonance: unchanged pencil, doubled linear term.
    pencil = sp.Matrix([[1,left],[right,1]])
    assert sp.expand(pencil.det()+x*(x+S)) == 0
    B = pencil.diff(x)
    assert B.det() == -1
    assert sp.expand((2*x+S).subs(x,-S)) == -S
    assert sp.expand((2*x+2*S).subs(x,-S)) == 0
    return dict(same_all_nonlinear_coefficients=True,same_constant_term=True,
        altered_linear_coefficient='S*x -> 2*S*x; S=z^3+z^(-3)',
        actual_reductive_Levi_representation=True,representation_dimension=82,
        ambient_simple_adjoint_representation=False,
        local_ideal='(u^2*v,u*v^2)',nonzero_nilpotent='u*v',
        both_local_ideal_inclusions_verified=True,
        all_original_derivatives_replayed_in_nonreduced_ring=True,
        generic_bad_target_map='-6*(z^6-z^(-6))',central_trace_norm=180,
        flat_root_laplacian_identity_preserved=True,
        full_ambient_Weyl_invariance_preserved=False,
        Weyl_invariance_failure={'weight':list(map(int,witness[0])),
                                'simple_node':witness[1]+1,
                                'reflected_weight':list(map(int,witness[2]))},
        invariant_derivations_not_ordinary_quotient_derivatives=True,
        sl2_adjoint_regular_unipotent_derivative=4,
        generic_adjoint_conjecture_disproved=False,
        coefficient_boundary_structure_alone_implies_rank_bound=False)


def main():
    archive = ROOT/'unit-constraint-cases.json.gz'
    rows = json.loads(gzip.decompress(archive.read_bytes()))
    counts = Counter()
    families = Counter()
    for record in rows:
        row = record['character']
        A = sp.Matrix(cartan(row['type']))
        n = A.rows
        I = [i-1 for i in row['nodes']]
        H = sp.Matrix(row['central_cocharacters']).T
        central_gram_inverse = (H.T*A*H).inv()
        positive = root_data(tuple(map(tuple,A.tolist())))[1]
        roots = set(positive) | {tuple(-v for v in a) for a in positive}
        simple = [tuple(A[:,i]) for i in I]
        direct = Counter()
        nonlinear = Counter()
        # Raising-operator test, without a character-polynomial expansion.
        for alpha in roots:
            if any(tuple(a+b for a,b in zip(alpha,beta)) in roots for beta in simple):
                continue
            charge = tuple(sum(alpha[i]*H[i,j] for i in range(n)) for j in range(H.cols))
            exponent = tuple(alpha[i] for i in I)
            if any(charge) and sum(exponent) in (2,3):
                assert all(v>=0 for v in exponent)
                nonlinear[exponent,charge] += 1
                if sum(exponent)==3:
                    direct[exponent,charge] += 1
        stored_nonlinear = Counter({(tuple(t['exponent']),tuple(t['charge'])):t['coefficient']
                          for t in row['terms'] if sum(t['exponent']) in (2,3) and any(t['charge'])})
        assert nonlinear == stored_nonlinear, (row['type'],I)
        all_grouped = defaultdict(list)
        for (exponent,charge),multiplicity in nonlinear.items():
            assert multiplicity == 1
            all_grouped[exponent].append(sp.Matrix(charge))
        nonlinear_coefficients(all_grouped,A,I,H,central_gram_inverse,roots,counts)
        stored = Counter({(tuple(t['exponent']),tuple(t['charge'])):t['coefficient']
                          for t in row['terms'] if sum(t['exponent'])==3})
        assert direct == stored, (row['type'],I)
        grouped = defaultdict(list)
        for (exponent,charge),multiplicity in direct.items():
            assert multiplicity == 1
            grouped[exponent].append(sp.Matrix(charge))
        singleton = {p[0] for p in row['factor_orders'] if len(p)==1}
        pure = []
        mixed = []
        for exponent,charges in grouped.items():
            support = {i for i,v in enumerate(exponent) if v}
            assert len(support)==3 and all(exponent[i]==1 for i in support)
            if not support<=singleton:
                assert len(charges)==1
                mixed.append(charges[0])
                counts['mixed_cubic_monomials'] += 1
                continue
            assert len(charges)==2 and charges[0] == -charges[1]
            q = charges[0]
            assert (q.T*central_gram_inverse*q)[0] == sp.Rational(1,2)
            gamma = A*H*central_gram_inverse*q
            delta = 2*gamma
            assert tuple(delta) in roots and all(delta[i]==0 for i in I)
            # Coroot is an actual integral cocharacter of simply connected G.
            delta_coroot = A.inv()*delta
            assert all(v.q==1 for v in delta_coroot)
            assert (gamma.T*delta_coroot)[0] == 1
            betas = [A[:,I[i]] for i in sorted(support)]
            axes = betas+[delta]
            d4 = {tuple(sign*b) for b in axes for sign in (-1,1)}
            d4 |= {tuple(sum((sign*b for sign,b in zip(signs,axes)),sp.zeros(n,1))/2)
                   for signs in product((-1,1),repeat=4)}
            assert len(d4)==24 and d4<=roots
            # Delta is an isolated A1 factor in the centralizer root system.
            for beta in roots:
                pairing = (sp.Matrix(beta).T*delta_coroot)[0]
                assert pairing in (-2,-1,0,1,2)
                if abs(pairing)==2:
                    assert tuple(beta) in (tuple(delta),tuple(-delta))
                reflected = sp.Matrix(beta)-pairing*delta
                assert tuple(reflected) in roots
                if any(beta[i] for i in I):
                    continue
                assert pairing==0 or tuple(beta) in (tuple(delta),tuple(-delta))
            pure.append((support,q,gamma,delta_coroot))
            counts['pure_A1_cubic_pairs'] += 1
            counts['D4_subsystems_checked'] += 1
        for (support,q,gamma,coroot),(other,p,eta,dual) in combinations(pure,2):
            assert len(support & other) in (0,2)
            assert (q.T*central_gram_inverse*p)[0] == 0
            assert (gamma.T*dual)[0] == (eta.T*coroot)[0] == 0
            counts['orthogonal_charge_pairs'] += 1
        assert len(pure)<=H.cols
        for mixed_charge in mixed:
            for _,charge,_,_ in pure:
                assert (mixed_charge.T*central_gram_inverse*(2*charge))[0] == 0
        families[row['type']] += 1
        counts['maximum_simultaneous_binomials'] = max(
            counts['maximum_simultaneous_binomials'],len(pure))

    # Exact coordinate identities, retaining both sheets of the torus map.
    q,u,v=sp.symbols('q u v',nonzero=True)
    a=q+1/q
    assert sp.cancel(q*sp.diff(a,q)-(q-1/q)) == 0
    assert sp.cancel((q-1/q)**2-(a*a-4)) == 0
    assert sp.cancel((q-1/q)-(q*q-1)/q) == 0
    assert sp.expand((q*q+1).subs(q,sp.I)) == 0
    assert sp.simplify((q-1/q).subs(q,sp.I)) == 2*sp.I
    cover = sp.groebner([q*q-u*q+1,v*q-1],v,q,u)
    assert cover.reduce((q-v)**2-(u*u-4))[1] == 0
    assert sp.discriminant(q*q-u*q+1,q) == u*u-4
    d,c = sp.symbols('d c',nonzero=True)
    equation = (u*u-4)*d*d-c*c
    reconstructed_q = (u+c/d)/2
    reconstructed_inverse = (u-c/d)/2
    assert sp.cancel(reconstructed_q*reconstructed_inverse-1-equation/(4*d*d)) == 0
    assert sp.cancel(reconstructed_q+reconstructed_inverse-u) == 0
    assert sp.cancel((reconstructed_q-reconstructed_inverse)*d-c) == 0
    assert dict(families)=={'E6':11,'E7':32,'E8':93}
    assert counts['mixed_cubic_monomials']==230
    assert counts['pure_A1_cubic_pairs']==78
    assert counts['orthogonal_charge_pairs']==48
    assert counts['nonlinear_coefficient_weight_sets']==1036
    assert sum(v for k,v in counts.items() if k.startswith('degree_2_')) == 728
    assert sum(v for k,v in counts.items() if k.startswith('degree_3_')) == 308
    assert counts['common_centralizer_coordinate_frames_checked'] == 136
    assert counts['common_frame_unit'] == 376
    assert counts['common_frame_central_binomial'] == 58
    assert counts['common_frame_1_fundamental_factors'] == 514
    assert counts['common_frame_2_fundamental_factors'] == 88
    check_coefficient_boundaries()
    resonance = check_linear_resonance(rows)
    report = dict(result='passed',uses_Sage=False,
        theorem_scope='Simply laced simple G, adjoint representation, proper all-A Levi; proof is uniform, checks are finite.',
        finite_cubic_character_records_checked=len(rows),records_by_type=dict(families),
        root_checks=dict(counts),
        independent_of_stored_character_polynomial_inversion=True,
        unit_or_binomial_cubic_coefficients_verified=True,
        D4_and_centralizer_A1_structure_verified=True,
        integral_coroot_character_pairings_verified=True,
        exact_etale_coordinate_identities_verified=True,
        full_character_unit_quadratic_coefficients_verified=True,
        scheme_elimination_inverse_verified=True,
        all_nonlinear_coefficient_simplex_or_cross_polytope_structure_verified=True,
        isolated_centralizer_factors_verified=True,
        individual_coefficient_boundary_formulas_verified=True,
        common_fundamental_character_coordinate_frames_verified=True,
        squarefree_coordinate_ideal_decompositions_verified=True,
        simultaneous_coefficient_zero_schemes_reduced_proved=True,
        linear_resonance_control=resonance,
        simultaneous_boundary_transversality_proved=False,
        generic_adjoint_smoothness_proved=False,new_Levi_coverage=False,
        literature_priority_claim=False,
        proof_note='UNIFORM_ADJOINT_PROOF_TARGET.md')
    sources = [Path(__file__),ROOT/report['proof_note'],
               ROOT/'tools/audit_affine_branching_certificates.py']
    report['source_sha256'] = {p.relative_to(ROOT).as_posix():
        sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in sources}
    report['raw_sha256'] = {archive.name:sha256(archive.read_bytes()).hexdigest()}
    (ROOT/'cubic-root-coefficients-audit.json').write_text(
        json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('source_sha256','raw_sha256')},indent=2))


if __name__=='__main__':
    main()
