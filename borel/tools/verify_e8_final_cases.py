"""Exact full-character and boundary checks for three further E8 families.

The other two archived models are reconstructed, not declared proved.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
from hashlib import sha256
import json, gzip
from sage.all import QQ, GF, PolynomialRing, LaurentPolynomialRing, RootSystem, matrix, vector, prod
from linear_character_elimination import as_laurent

ROOT=Path(__file__).resolve().parents[1]
def encode(f):return [[list(map(int,e)),str(c)] for e,c in f.dict().items()]
def conv(a,b):
    out=Counter()
    for x,c in a.items():
        for y,d in b.items():out[tuple(u+v for u,v in zip(x,y))]+=c*d
    return out

def check_models():
    rows=json.loads((ROOT/'e8-final-character-models.json').read_text())
    lattice=RootSystem('E8').root_lattice();A=lattice.cartan_type().cartan_matrix()
    roots=[vector(QQ,[r[i+1] for i in range(8)]) for r in lattice.roots()]
    pairs={'1257':None,'12356':([0],[2]),'12357':([0],[2]),
           '12457':([1],[3]),'123567':([3,4],[0,5])}
    for row in rows:
        H=matrix(QQ,row['columns']).transpose();assert H.det()!=0
        k=row['derived']+row['central_semisimple'];n=row['derived']
        for j in range(row['derived'],8):
            assert all((A*H.column(j))[i-1]==0 for i in row['seed'])
        for j in range(k):assert H.column(j) in roots
        M=H.transpose()*A*H
        for comp in row['components']:
            for i,p in enumerate(comp):
                for j,q in enumerate(comp):assert M[p,q]==(2 if i==j else -int(abs(i-j)==1))
        fundamental={}
        for path in row['components']:
            N=len(path)+1;letters=[]
            for j in range(N):
                w=[0]*8
                for p,i in enumerate(path):w[i]=int(j==p)-int(j==p+1)
                letters.append(w)
            for j,i in enumerate(path):
                fundamental[i]=Counter(tuple(sum(letters[t][v] for t in subset) for v in range(8))
                                       for subset in combinations(range(N),j+1))
        actual=Counter({tuple(r*A*H):1 for r in roots});actual[(0,)*8]=8
        rebuilt=Counter()
        for term in row['terms']:
            start=(0,)*k+tuple(term['charge']);weights={start:term['coefficient']}
            for i,e in enumerate(term['exponent']):
                for _ in range(e):weights=conv(weights,fundamental[i])
            rebuilt.update(weights)
        assert {w:c for w,c in rebuilt.items() if c}==dict(actual)
        tor=row['variables'][k:];L=LaurentPolynomialRing(QQ,tor if tor else ['dummy'])
        P=PolynomialRing(L.fraction_field(),names=row['variables'][:k]);v=P.gens()
        phi=sum(t['coefficient']*prod(w**e for w,e in zip(v,t['exponent']))*
                prod(w**e for w,e in zip(L.gens(),t['charge'])) for t in row['terms'])
        pair=pairs[''.join(map(str,row['seed']))];removed=[];graph={}
        if pair:
            taus=[v[i] for i in pair[0]];aas=[v[i] for i in pair[1]]
            constraints=[phi.derivative(t) for t in taus]
            C=matrix(P,[[f.derivative(a) for a in aas] for f in constraints]);det=C.det()
            assert as_laurent(L,P.base_ring()(det)).is_unit()
            assert all(c.derivative(a)==0 for c in C.list() for a in aas+taus)
            rhs=vector(P,[sum(C[i,j]*aas[j] for j in range(len(aas)))-constraints[i] for i in range(len(taus))])
            assert all(f.derivative(a)==0 for f in rhs for a in aas+taus)
            Ci=C.adjugate()/det;graph=dict(zip(aas,Ci*rhs));G=phi.subs({t:0 for t in taus})
            graph.update(dict(zip(taus,-Ci.transpose()*vector(P,[G.derivative(a).subs(graph) for a in aas]))))
            assert all(phi.derivative(a).subs(graph)==0 for a in aas+taus)
            removed=pair[0]+pair[1]
        residual=P(phi.subs(graph));retained=[i for i in range(k) if i not in removed]
        expected=sum(QQ(t['coefficient'])*prod(v[i]**e for i,e in zip(retained,t['exponent']))*
                     prod(w**e for w,e in zip(L.gens(),t['charge'])) for t in row['residual_terms'])
        assert residual==expected
        for i in retained:assert phi.derivative(v[i]).subs(graph)==residual.derivative(v[i])
        for z in L.gens():
            D=lambda f:P({ex:z*c.derivative(z) for ex,c in P(f).dict().items()})
            assert D(phi).subs(graph)==D(residual)
        row['exact_root_reconstruction_and_all_chain_rules']=True
    return rows

def binary(row):
    P=PolynomialRing(QQ,names='x,y,a,t');x,y,a,t=P.gens();K=P.fraction_field()
    f=K(0)
    for q in row['residual_terms']:
        e=q['exponent'];power=q['charge'][0]-2*e[1];assert power%6==0
        f+=QQ(q['coefficient'])*x**e[0]*y**e[1]*a**e[2]*t**(power//6)
    aa=x*x-(t+1/t)*a*x+a*a-1
    bb=a*x**3-2*t*a*a*x*x+a**3*x+(t-1/t)*x*x-a*x+(t-1/t)*(a*a-1)
    cc=-t*a**3*x**3+t*a*x**3+a*a*x*x+t*a**3*x-(2*t+1/t)*a*x-3
    assert f==aa*y*y+bb*y+cc
    N=P(t**5*(4*aa**2*cc.derivative(x)-2*aa*bb*bb.derivative(x)+bb**2*aa.derivative(x)))
    assert N.degree(x)==7 and N.coefficient({x:7})==-4*a*a*t**5
    assert gcd_list([N.coefficient({x:i}) for i in range(8)])==t*t
    assert P(aa*t).gcd(P(bb*t))==1
    U=PolynomialRing(QQ,'X');X=U.gen();special=U(N(x=X,y=0,a=3,t=2))
    assert special.degree()==7 and special.gcd(special.derivative())==1
    x0=(t+1/t)*a/2
    F=P(aa.subs({x:x0}).numerator());G=P(bb.subs({x:x0}).numerator())
    res=F.resultant(G,a)
    assert res==65536*(t-1)**2*(t+1)**2*t**4*(t*t+1)**4*(t**4+1)**2
    # On every possible t, F is either a nonzero constant or a nonzero
    # quadratic in a. Thus the resultant bound has no positive-dimensional exception.
    assert F.subs({t:1})!=0 and F.subs({t:-1})!=0
    return dict(seed='123567',binary_degree=7,content='t^2 (Laurent unit)',
                squarefree_specialization=encode(P(special(x))),exceptional_boundary_dimension_bound=0,
                full_boundary_dimension_bound=1,regular_fiber_dimension_if_nonempty=20)

def gcd_list(fs):
    g=fs[0]
    for f in fs[1:]:g=g.gcd(f)
    return g

def pair(row):
    P=PolynomialRing(QQ,names='x,y,z,a,b,t');x,y,z,a,b,t=P.gens();K=P.fraction_field();f=K(0)
    for q in row['residual_terms']:
        e=q['exponent'];power=q['charge'][0]+e[1]-e[2]+e[3]-e[4];assert power%3==0
        f+=QQ(q['coefficient'])*prod(w**i for w,i in zip((x,y,z,a,b),e))*t**(power//3)
    H=matrix(K,[[f.derivative(v).derivative(w) for v in (y,z)] for w in (y,z)])
    ell=vector(K,[f.derivative(v).subs({y:0,z:0}) for v in (y,z)]);Q=f.subs({y:0,z:0});D=H.det()
    N=Q*D-ell*H.adjugate()*ell/2
    J=P(t**6*(N.derivative(x)*D-N*D.derivative(x)))
    assert J.degree(x)==9 and J.coefficient({x:9})==2*t**6
    U=PolynomialRing(QQ,'X');X=U.gen();sp=U(J(x=X,y=0,z=0,a=2,b=3,t=2))
    assert sp.degree()==9 and sp.gcd(sp.derivative())==1
    R=PolynomialRing(QQ,names='x,a,b,t,u',order='degrevlex');X,A,B,T,U=R.gens()
    convert=lambda z:R({(e[0],e[3],e[4],e[5],0):c for e,c in P(z.numerator()).dict().items()})
    certificates=[]
    for j in (0,1):
        C=H[j,j]*ell[1-j]-H[j,1-j]*ell[j];unit=convert(H[j,j])*T
        gens=[convert(D),convert(C),U*unit-1]
        for extra,dimension in (([],2),([convert(D.derivative(x))],1)):
            I=R.ideal(gens+extra);basis=list(I.groebner_basis(algorithm='libsingular:std'))
            assert I.dimension()==dimension
            certificates.append(dict(chart=j,dimension=dimension,variables=list(R.variable_names()),
                                     generators=[encode(g) for g in gens+extra],basis=[encode(g) for g in basis]))
    zero=R.ideal([convert(z) for z in H.list()+list(ell)]+[U*T-1])
    target=R.ideal([A*A*B-A,A*B*B-B,A*U*U-A,B*U*U-B,
                   X*X+A*B+X*T+X*U,X*A+A*U,X*B+B*U,A*T-A*U,B*T-B*U,T*U-1])
    assert zero==target and zero.dimension()==1
    assert f.derivative(x).subs({a:0,b:0})==(2*x+t+1/t)*(1-y*z)
    assert H[0,0].derivative(x)==-2*b and H[1,1].derivative(x)==-2*a
    # Rank-one symmetric 2x2 matrices always have a nonzero diagonal.
    # Rank zero: a=b=0, x=0 or -t-t^-1; or ab=1,t=+/-1,x=-t.
    return dict(seed='12356',schur_numerator_degree=9,leading_coefficient='2t^6 (unit)',
                squarefree_specialization=True,rank_one_chart_dimensions=[2,2],
                rank_one_identically_quadratic_exception_bounds=[1,1],rank_zero_base_dimension=1,
                full_boundary_dimension_bound=2,regular_fiber_dimension_if_nonempty=14,
                boundary_certificates=certificates)

def multiplier(row):
    P=PolynomialRing(QQ,names='x,y,v,a,b,c');x,y,v,a,b,c=P.gens()
    f=-x*x*y*y-c*x*y*y*a-y*y*a*a-c*x*x*y*b-3*x*y*a*b-c*y*a*a*b-x*x*b*b-c*x*a*b*b-a*a*b*b+v*(a*x+b*y+c)+x*x+y*y+c*x*a+a*a+c*y*b+b*b-4
    # Check the entire residual against the original central torus; no denominator boundary is removed.
    L=LaurentPolynomialRing(QQ,'z');z=L.gen();S=PolynomialRing(L,names='x,y,v,a,b');xx,yy,vv,aa,bb=S.gens()
    expected=S(sum(QQ(q['coefficient'])*prod(w**e for w,e in zip((xx,vv,yy,aa,bb),q['exponent']))*z**q['charge'][0] for q in row['residual_terms']))
    assert P.hom([xx,yy,vv,aa,bb,S(z*z+z**-2)],S)(f)==expected
    g=[f.derivative(w).subs({c:-a*x-b*y}) for w in (x,y)]
    Q=PolynomialRing(QQ,names='x,y,v,a,b,h');h=Q.gen(5)
    gs=[Q({tuple(e[:5])+(0,):co for e,co in z.dict().items()}) for z in g]
    J=matrix(Q,[[z.derivative(w) for w in Q.gens()[:-1]] for z in gs])
    original=[z.homogenize(h) for z in gs+J.minors(2)]
    reports=[]
    for prime,algorithm in ((32003,'libsingular:slimgb'),(32009,'libsingular:std')):
        F=Q.change_ring(GF(prime));I=F.ideal([F(z) for z in original]);G=list(I.groebner_basis(algorithm=algorithm))
        assert I.dimension()==3
        reports.append(dict(prime=prime,algorithm=algorithm,basis_size=len(G),homogeneous_dimension=3))
    certificate=dict(variables=list(Q.variable_names()),original_generators=[encode(z) for z in original],
                     prime=32003,dimension=3)
    return dict(seed='12457',homogeneous_dimension_bounds=reports,actual_singular_dimension_bound=2,
                regular_fiber_dimension_if_nonempty=16,certificate=certificate)

def last_model(row):
    L=LaurentPolynomialRing(QQ,names='U,V');U,V=L.gens()
    P=PolynomialRing(L,names='x,y,w,a');x,y,w,a=P.gens()
    A=U+1/U;B=V+1/V;C=U*V+1/(U*V)
    assert A*A+B*B+C*C-A*B*C-4==0
    f=x*x+y*y+w*w+a*a-x*x*y*y-x*x*w*w-y*y*w*w-x*x*a*a-y*y*a*a-w*w*a*a-C*(x*x*y*w+x*y*y*a+x*w*w*a+y*w*a*a)-B*(x*y*y*w+x*x*y*a+y*w*w*a+x*w*a*a)-A*(x*y*w*w+x*x*w*a+y*y*w*a+x*y*a*a)-5*x*y*w*a-B*C*(x*y+w*a)-A*C*(x*w+y*a)-A*B*(y*w+x*a)-A*B*C-4
    actual=P(0)
    for q in row['residual_terms']:
        r,s=q['charge'];assert r%3==0
        actual+=QQ(q['coefficient'])*prod(v**e for v,e in zip(P.gens(),q['exponent']))*U**(-r//3-2*s)*V**(r//3+s)
    assert f==actual
    assert matrix(QQ,[[3,-1],[6,-1]]).det()==3
    for u in (-1,1):
        for v in (-1,1):
            assert all((z*q.derivative(z))(U=u,V=v)==0 for q in (A,B,C) for z in (U,V))
    return dict(full_residual_verified=True,torus_isogeny_degree=3,
                four_trace_quotient_branch_points_have_two_zero_central_derivatives=True,
                smoothness_proved=False)

def main():
    rows=check_models();lookup={''.join(map(str,r['seed'])):r for r in rows}
    cases=[binary(lookup['123567']),pair(lookup['12356']),multiplier(lookup['12457'])]
    four=[]
    for prime in (32003,32009):
        r=json.loads((ROOT/f'e8-a1-four-{prime}-verification.json').read_text())
        assert r['result']=='passed' and r['characteristic_zero_affine_dimension_upper_bound']==3
        four.append(r)
    cases.append(dict(seed='1257',affine_singular_dimension_bound=3,regular_fiber_dimension_if_nonempty=8,
                      homogeneous_certificate_primes=[r['prime'] for r in four]))
    report=dict(result='passed',full_characters_and_scheme_eliminations_checked=5,cases=cases,
                newly_completed_classes=['A3+A2+A1','A2^2+A1','A3+A1^2','A1^4'],
                unresolved_classes=['A2+A1^3'],covered_proper_E8_standard_subsets=247,
                remaining_standard_subsets=8,remaining_Weyl_classes=1,
                remaining_model=last_model(lookup['12357']),
                all_E8_Levis_resolved=False,nonemptiness_claim=False,literature_priority_claim=False)
    names=['E8_FINAL_CASES.md','e8-final-character-models.json','tools/verify_e8_final_cases.py',
           'BINARY_BRANCHING_FIBERS.md','e8-weyl-transfers-verification.json',
           'e8-a1-four-32003-verification.json','e8-a1-four-32009-verification.json']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'e8-final-cases-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('cases','source_sha256')},indent=2))
if __name__=='__main__':main()

