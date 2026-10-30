"""Exact G2 vertical-rank audit, including independent 7-dimensional matrices.

Run from borel with PYTHONPATH=src:tests sage -python tools/verify_g2_borel_rank.py.
No sampling is used to establish the stated rank strata: the proof uses
centralizer reduction and the exact commutator matrices checked here.
"""
from hashlib import sha256
import json
from pathlib import Path
import sys
from sage.all import (QQ, PolynomialRing, NumberField, matrix, vector,
                      diagonal_matrix, identity_matrix, gcd, prod)
from sage.env import SAGE_VERSION
from sage_borel import root_borel_map

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tests'))
from test_root_borel import adjoint_data, exponential_nilpotent, trace_pair


def minimal_data(api, torus):
    # The seven weights form the lowering chain 1,2,1,1,2,1.
    # On each sl2 string F has entries 1 and E has k(m-k+1).
    weights = [(1,0),(-1,1),(2,-1),(0,0),(-2,1),(1,-1),(-1,0)]
    E = [matrix(QQ, 7, sparse=True) for _ in range(2)]
    F = [matrix(QQ, 7, sparse=True) for _ in range(2)]
    for node, entries in enumerate([[(0,1),(2,2),(3,2),(5,1)],[(1,1),(4,1)]]):
        for j, c in entries:
            E[node][j,j+1] = c
            F[node][j+1,j] = 1
    H = [e*f-f*e for e,f in zip(E,F)]
    assert H == [diagonal_matrix([w[i] for w in weights]) for i in range(2)]
    for i in range(2):
        for j in range(2):
            assert E[i]*F[j]-F[j]*E[i] == (H[i] if i == j else matrix(QQ,7))
            assert H[i]*E[j]-E[j]*H[i] == api.root_weights[j][i]*E[j]
            assert H[i]*F[j]-F[j]*H[i] == -api.root_weights[j][i]*F[j]
            if i != j:
                for generators in [E,F]:
                    out = generators[j]
                    for _ in range(1-api.root_weights[j][i]):
                        out = generators[i]*out-out*generators[i]
                    assert out.is_zero()
    positive, negative = list(E), list(F)
    for alpha in api.positive_roots[2:]:
        i,j = next((i,j) for i,j in api.brackets
                   if tuple(a+b for a,b in zip(api.positive_roots[i],api.positive_roots[j])) == alpha)
        positive.append((positive[i]*positive[j]-positive[j]*positive[i])/api.brackets[i,j])
        negative.append(negative[i]*negative[j]-negative[j]*negative[i])
    index = {r:i for i,r in enumerate(api.positive_roots)}
    for i in range(6):
        for j in range(i+1,6):
            k = index.get(tuple(a+b for a,b in zip(api.positive_roots[i],api.positive_roots[j])))
            assert positive[i]*positive[j]-positive[j]*positive[i] == (
                positive[k]*api.brackets[i,j] if k is not None else matrix(QQ,7))
    assert matrix(QQ, [m.list() for m in H+positive+negative]).rank() == 14
    T = diagonal_matrix([prod(t**a for t,a in zip(torus,w)) for w in weights])
    return T, positive, negative, H


def full_trace_check(api, minimal):
    forms = api.symbolic_coordinates()
    x,y = api.torus_parameters
    us = api.unipotent_parameters
    field = api.symbolic_ring.fraction_field()
    T,E,F,H = (minimal_data if minimal else adjoint_data)(api,(x,y))
    g = T.change_ring(field)
    for e,u in zip(E,us):
        g *= exponential_nilpotent(e.change_ring(field)*u)
    actual = tuple(trace_pair(g,f)/trace_pair(e,f) for e,f in zip(E,F))
    assert actual == forms
    gram = matrix(QQ, [[trace_pair(h,k) for k in H] for h in H])
    observed = gram.change_ring(field).solve_right(vector(field,[trace_pair(g,h) for h in H]))
    assert tuple(observed) == tuple(field(h(x,y)) for h in api.toral_coordinates())
    assert all(trace_pair(g,e) == 0 for e in E)
    return forms


def specialize(api, forms, point, values=None):
    # Construct a homomorphism to avoid Python integer negative powers/floats.
    ground = point[0].parent()
    ring = PolynomialRing(ground, names=['a','b','c','d','e','f'])
    us = ring.gens()
    hom = api.symbolic_ring.hom(list(point)+list(us),ring)
    out = [ring(hom(f.numerator())/hom(f.denominator())) for f in forms]
    jac = matrix([[f.derivative(u) for u in us] for f in out])
    if values is not None:
        sub = dict(zip(us,values))
        return tuple(f.subs(sub) for f in out), jac.apply_map(lambda f:f.subs(sub))
    return out,jac


def commutator_matrix(api, coefficients):
    ring = coefficients[0].parent()
    M = matrix(ring,6)
    index = {r:i for i,r in enumerate(api.positive_roots)}
    for (i,j),c in api.brackets.items():
        k = index[tuple(a+b for a,b in zip(api.positive_roots[i],api.positive_roots[j]))]
        M[k,j] += coefficients[i]*c
        M[k,i] -= coefficients[j]*c
    return M


def main():
    univar = PolynomialRing(QQ,'s'); s = univar.gen(); K = univar.fraction_field()
    curves = [(s,s**2),(s**2,s**3),(s,s),(K(1),s),(s,s**3),(s,K(1))]
    involutions = [(-1,1),(1,-1),(-1,-1)]
    central_indices = [(0,5),(3,1),(2,4)] # short, long, zero based
    records = []
    point_rank_checks = 0
    for minimal in [True,False]:
        api = root_borel_map(('G',2),(1,0) if minimal else (0,1))
        forms = full_trace_check(api,minimal)
        print(('minimal' if minimal else 'adjoint')+': full symbolic trace agrees',flush=True)
        qshort = (s*s+4*s+1)/(6*s) if minimal else (5*s*s+2*s+5)/(12*s)
        qlong = (s*s+1)/(2*s) if minimal else (s+1)**2*(s**4-2*s**3+4*s*s-2*s+1)/(8*s**3)
        hs = (s-1)*(s+1) if minimal else (s-1)*(s+1)**3
        hl = (s-1)*(s+1)*(s*s+s+1)
        if not minimal:
            hl *= s*s-QQ(2)/3*s+1
        curve_records = []
        for i,pt in enumerate(curves):
            pt = tuple(K(c) for c in pt)
            h = [K(f(*pt)) for f in api.toral_coordinates()]
            q = K(api.root_factors()[i](*pt))
            short = i in (0,2,3)
            assert q == (qshort if short else qlong)
            common = gcd([f.numerator() for f in h]).monic()
            assert common == (hs if short else hl).monic()
            # Exhaust all intersections with the other five hypertori.
            intersection = univar(1)
            for j,weight in enumerate(api.root_weights):
                if j != i:
                    equation = K(prod(t**a for t,a in zip(pt,weight))-1).numerator()
                    intersection *= equation
            squarefree = prod(p for p,multiplicity in intersection.factor()
                              if p != s).monic()
            expected_intersection = (s-1)*(s+1)
            if not short:
                expected_intersection *= s*s+s+1
            assert squarefree == expected_intersection
            curve_records.append(dict(root=list(api.positive_roots[i]),parameter=list(map(str,pt)),
                                      central_factor=str(q.factor()),toral_zero_polynomial=str(common.factor())))
        # Exact symbolic rank at every centralizer-normal-form point of each
        # involution. The generic matrices are also checked on all zero strata.
        involution_records = []
        P = PolynomialRing(QQ,names=['A','B']); A,B = P.gens()
        for pt,(short,long) in zip(involutions,central_indices):
            pt = tuple(QQ(v) for v in pt)
            q = api.root_factors(pt)
            assert q[short] == (QQ(1)/3 if minimal else -QQ(2)/3)
            assert q[long] == (-1 if minimal else 0)
            assert all(q[i] == 0 for i in range(6) if i not in (short,long))
            coefficients = [P(0)]*6
            coefficients[short] = q[short]*A
            coefficients[long] = q[long]*B
            columns = [short,long] if minimal else [short]
            M = commutator_matrix(api,coefficients).augment(identity_matrix(P,6).matrix_from_columns(columns))
            ranks = {}
            for av,bv in [(0,0),(0,1),(1,0),(1,1)]:
                values = [QQ(0)]*6; values[short]=QQ(av); values[long]=QQ(bv)
                _,J = specialize(api,forms,pt,values)
                rank = M.apply_map(lambda p:p(av,bv)).rank()
                assert J.rank() == rank
                ranks[str((av,bv))] = rank
                point_rank_checks += 1
            expected = ([2,2,5,5] if short==0 else [2,4,4,4] if short==3 else [2,3,4,4]) if minimal else (
                [1,1,4,4] if short==0 else [1,1,3,3])
            assert list(ranks.values()) == expected
            # The nonzero entries are monomials in A,B. Verify all minors
            # generating the relevant rank conditions, beyond point samples.
            for rank in range(1,7):
                minors = [v for v in M.minors(rank) if v]
                ideal = P.ideal(minors)
                if rank > M.rank():
                    assert ideal.is_zero()
                elif rank <= (2 if minimal else 1):
                    assert ideal == P.ideal(1)
                else:
                    expected_radical = P.ideal(A) if short==0 or not minimal or (short==2 and rank==4) else P.ideal(A,B)
                    assert ideal.radical() == expected_radical
            involution_records.append(dict(torus=list(map(str,pt)),short_root=short+1,long_root=long+1,
                                           ranks=ranks,determinantal_ideals_verified=True))
        # Long A2 centralizer at a primitive cube root, both embeddings covered
        # by a polynomial identity in the quadratic number field.
        L = NumberField(s*s+s+1,'zeta'); zeta = L.gen(); pt=(zeta,L(1))
        assert all(f(*pt)==0 for f in api.toral_coordinates())
        assert [f(*pt) for f in api.root_factors()] == [0, -QQ(1)/2 if minimal else QQ(5)/8,0,0,
                                                     -QQ(1)/2 if minimal else QQ(5)/8,
                                                     -QQ(1)/2 if minimal else QQ(5)/8]
        out,J = specialize(api,forms,pt)
        aa,bb,cc,dd,ee,ff = J.base_ring().gens()
        sub = {aa:0,cc:0,dd:0}
        N = J.apply_map(lambda f:f.subs(sub))
        # All 5-minors vanish identically; the 4-minor ideal has radical (b).
        P3=PolynomialRing(L,names=['b','e','f']); b,e,f=P3.gens()
        hom=J.base_ring().hom([0,b,0,0,e,f],P3)
        N=N.apply_map(hom)
        assert N.rank()==4
        assert P3.ideal(N.minors(4)).radical()==P3.ideal(b)
        assert P3.ideal(N.minors(3))==P3.ideal(1)
        extra=[]
        if not minimal:
            L=NumberField(3*s*s-2*s+3,'z'); z=L.gen()
            for root,pt in [(1,(z*z,z**3)),(4,(z,z**3)),(5,(z,L(1)))]:
                assert all(f(*pt)==0 for f in api.toral_coordinates())
                q=[f(*pt) for f in api.root_factors()]
                assert q[root]==QQ(10)/27 and all(q[i]==0 for i in range(6) if i!=root)
                out,J=specialize(api,forms,pt)
                us=J.base_ring().gens(); param=us[root]
                N=J.apply_map(lambda f:f.subs({u:0 for i,u in enumerate(us) if i!=root}))
                top={1:3,4:2,5:1}[root]
                assert N.rank()==top
                if top>1:
                    assert J.base_ring().ideal(N.minors(top)).radical()==J.base_ring().ideal(param)
                assert J.base_ring().ideal(N.minors(1))==J.base_ring().ideal(1)
                extra.append(dict(root=root+1,nonzero_parameter_rank=top,zero_parameter_rank=1))
        records.append(dict(representation='minimal' if minimal else 'adjoint',curves=curve_records,
                            involutions=involution_records,order_three_rank='4 if b != 0; 3 otherwise',
                            extra_long_A1=extra))
    paths=['G2_BOREL_RANK.md','tools/verify_g2_borel_rank.py','src/sage_borel/root_borel.py',
           'src/sage_borel/characters.py','src/sage_borel/borel_reconstruction.py',
           'src/sage_borel/geometry.py','src/sage_borel/__init__.py','tests/test_root_borel.py']
    report=dict(result='passed',sage_version=SAGE_VERSION,coefficient_field='characteristic zero; exact QQ and quadratic fields',
                full_symbolic_trace_representations=[7,14],centralizer_curves_checked=12,
                rational_normal_form_rank_checks=point_rank_checks,
                determinantal_strata_verified=True,root_hypertorus_intersections_exhausted=True,records=records,
                theorem_scope='All G2 vertical ranks, using torus evaluation off the root hypertori and centralizer normal forms on them',
                general_deduction='rank D_U theta(tz) = dim([u,theta(tz)] + d(theta|tU_t)_z(T_z U_t))',
                independence='Minimal matrices are built from explicit sl2 weight strings and checked against Chevalley-Serre relations. Adjoint trace uses full Lie matrices. Both avoid the fixed-vector recurrence; sign alignment and the adjoint route share Sage root constants. Not a separate CAS or external review.',
                limitations='No full B/G differential, scheme-theoretic fiber or target branch-locus classification; no novelty claim.',
                source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'g2-borel-rank-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('source_sha256','records')},indent=2))


if __name__=='__main__':
    main()
