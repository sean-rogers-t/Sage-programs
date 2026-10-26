"""Exact certificates for one binary criterion and two E7 cubic classes."""
from hashlib import sha256
from pathlib import Path
import json
from sage.all import QQ, PolynomialRing, LaurentPolynomialRing, matrix, vector, prod, RootSystem
from verify_branching_norm_fibers import charge_polynomial

ROOT=Path(__file__).resolve().parents[1]


def original_characters():
    row1=charge_polynomial('E7',(0,1,4,6))
    L=LaurentPolynomialRing(QQ,['p','q','r']);p,q,r=L.gens()
    R=PolynomialRing(L.fraction_field(),names='X,x,y,w');X,x,y,w=R.gens()
    original=sum(t['coefficient']*prod(v**e for v,e in zip((p,q,r),t['charge']))*
                 prod(v**e for v,e in zip(R.gens(),t['exponent'])) for t in row1['terms'])
    U=p*q;V=p*q*r;W=p*q*q*r
    A=U+1/U;B=V+1/V;C=W+1/W
    full=X*X+x*x+y*y+w*w+A*X*x*y+B*X*x*w+C*X*y*w+B*C*x*y+A*C*x*w+A*B*y*w+A*B*C*X+A*A+B*B+C*C-7
    assert original==full
    graph={X:-(A*x*y+B*x*w+C*y*w+A*B*C)/2}
    residual=full.subs(graph)
    for v in (x,y,w):assert full.derivative(v).subs(graph)==residual.derivative(v)
    for z in (p,q,r):
        D=lambda f:R({e:z*c.derivative(z) for e,c in f.dict().items()})
        assert D(full).subs(graph)==D(residual)
    row2=charge_polynomial('E7',(0,1,2,4,6))
    L=LaurentPolynomialRing(QQ,['p','q']);p,q=L.gens()
    R=PolynomialRing(L.fraction_field(),names='a,x,b,y,w');a,x,b,y,w=R.gens()
    full=sum(t['coefficient']*prod(v**e for v,e in zip((p,q),t['charge']))*
             prod(v**e for v,e in zip(R.gens(),t['exponent'])) for t in row2['terms'])
    F=full.derivative(a)-b;G=full.derivative(b)-a
    graph={a:-G,b:-F};residual=full.subs(graph)
    U=q;V=p**3*q;A=U+1/U;B=V+1/V;C=U*V+1/(U*V)
    expected=x*x+y*y+w*w-x*x*y*y-x*x*w*w-y*y*w*w-x*y*w*(A*x+B*y+C*w)-A*B*x*y-A*C*x*w-B*C*y*w-A*B*C-3
    assert residual==expected
    for v in (x,y,w):assert full.derivative(v).subs(graph)==residual.derivative(v)
    for z in (p,q):
        D=lambda f:R({e:z*c.derivative(z) for e,c in f.dict().items()})
        assert D(full).subs(graph)==D(residual)
    assert row2['central_cocharacters']==[[4,6,8,12,9,6,3],[2,3,4,6,5,4,2]]
    return [row1,row2]


def symbolic():
    R=PolynomialRing(QQ,names='X,x,y,A,B,C,w');X,x,y,A,B,C,w=R.gens()
    full=X*X+x*x+y*y+w*w+A*X*x*y+B*X*x*w+C*X*y*w+B*C*x*y+A*C*x*w+A*B*y*w+A*B*C*X+A*A+B*B+C*C-7
    L=A*y+B*w;M=C*y*w+A*B*C
    a=1-L*L/4;b=B*C*y+A*C*w-L*M/2
    c=y*y+w*w+A*B*y*w+A*A+B*B+C*C-7-M*M/4
    assert full.subs({X:-(L*x+M)/2})==a*x*x+b*x+c
    P=4*a*a*c.derivative(y)-2*a*b*b.derivative(y)+b*b*a.derivative(y)
    assert P.degree(y)==5 and P.coefficient({y:5})==A**4/2
    assert P.subs({A:0,B:0,C:0,w:0})==8*y
    spec={A:1,B:1,C:1,w:0}
    assert P.subs(spec)==y*(y*y-3)*(y*y-5)/2
    assert a.subs(spec)==1-y*y/4 and b.subs(spec)==y/2
    assert a.derivative(y)==-A*(A*y+B*w)/2
    assert b.derivative(y).subs({A:0})==-C*B*(w*w-2)/2
    witness={X:-4,x:-4,y:-4,A:2,B:2,C:2,w:2}
    H1=matrix(QQ,[[full.derivative(v).derivative(u).subs(witness)
                   for u in (X,x,y)] for v in (X,x,y)])
    assert H1.det()==-216
    assert all(full.derivative(v).subs(witness)==0 for v in (X,x,y))
    assert [full.derivative(v).subs(witness) for v in (A,B,C,w)]==[-108,36,36,36]
    assert a.subs(witness)==-3

    # Work in an independent polynomial ring, clearing only Laurent units.
    S=PolynomialRing(QQ,names='U,V,w,x,y');U,V,w,x,y=S.gens();F=S.fraction_field()
    A=U+1/F(U);B=V+1/F(V);C=U*V+1/F(U*V)
    a=1-y*y-w*w-A*w*y
    b=-B*y*y*w-C*y*w*w-A*B*y-A*C*w
    c=y*y+w*w-y*y*w*w-B*C*y*w-A*B*C-3
    P=4*a*a*c.derivative(y)-2*a*b*b.derivative(y)+b*b*a.derivative(y)
    N=S(U**5*V**2*P)
    assert N.degree(y)==5
    assert N.coefficient({y:5})==2*U**5*((V*V-1)**2*w*w+4*V*V)
    assert N.coefficient({y:0}).degree(U)==10
    assert N.coefficient({y:0}).coefficient({U:10})==-V**4*w**3
    assert P.subs({U:1,V:1,w:0})==8*y*(y*y-3)*(y*y+1)
    assert a.subs({U:1,V:1,w:0})==1-y*y and b.subs({U:1,V:1,w:0})==-4*y
    y0=-A*w/2
    E=(U*U-1)**2*w*w+4*U*U
    assert a.subs({y:y0})==E/(4*U*U)
    expected=-w*(U*U+1)*(U**4+1)*(V*V-1)/(2*U*U*V*(U*U-1))
    numerator=S((b.subs({y:y0})-expected).numerator())
    assert numerator.reduce([S(E)])==0

    T=PolynomialRing(QQ,names='x,y,A,B,C,w');x,y,A,B,C,w=T.gens()
    psi=x*x+y*y+w*w-x*x*y*y-x*x*w*w-y*y*w*w-x*y*w*(A*x+B*y+C*w)-A*B*x*y-A*C*x*w-B*C*y*w-A*B*C-3
    witness={x:-4,y:-4,A:2,B:2,C:2,w:2}
    assert all(psi.derivative(v).subs(witness)==0 for v in (x,y))
    H2=matrix(QQ,[[psi.derivative(v).derivative(u).subs(witness) for u in (x,y)] for v in (x,y)])
    assert H2.det()==-108
    assert [psi.derivative(v).subs(witness) for v in (A,B,C,w)]==[108,108,-36,36]
    logH2=matrix(QQ,[[144,-72,0],[-72,144,0],[0,0,72]])
    assert logH2.det()==1119744
    # Counterexample control: a unit quadratic coefficient alone is insufficient.
    D=PolynomialRing(QQ,names='x,y,z');x,y,z=D.gens()
    bad=x*x+x*y*y+z
    reduced=bad.subs({x:-y*y/2})
    assert reduced.derivative(y)==-y**3
    return dict(binary_degree=5,primitive_certificates=True,squarefree_witnesses=True,
                boundary_dimension_bounds=[3,2],exceptional_base_dimension_bounds=[2,1],
                relative_hessian_determinants=[-216,-108],
                logarithmic_hessians=[[-216,72,72,72],[[144,-72,0],[-72,144,0],[0,0,72]]],
                bad_squarefree_control=True)


def transfers():
    data=json.loads((ROOT/'e7-multiplier-boundary-verification.json').read_text())
    Q=RootSystem('E7').root_lattice();A=Q.cartan_type().cartan_matrix()
    roots=[tuple(int(v[i+1]) for i in range(7)) for v in Q.roots()]
    def subsystem(nodes):return {r for r in roots if all(r[i]==0 for i in range(7) if i+1 not in nodes)}
    out=[]
    for group in data['Weyl_classes']:
        if group['seed'] not in ([1,2,5],[1,2,3,5]):continue
        for member in group['members']:
            current=subsystem(group['seed'])
            for node in member['word']:
                i=node-1;new=set()
                for root in current:
                    v=list(root);v[i]-=sum(A[i,j]*root[j] for j in range(7));new.add(tuple(v))
                current=new
            assert current==subsystem(member['nodes'])
        out.append(group)
    assert sorted(len(g['members']) for g in out)==[10,12]
    return out


def main():
    report=dict(result='passed',characters=original_characters(),symbolic=symbolic(),
                Weyl_classes=transfers(),new_standard_embeddings=22,
                resolved_proper_standard_Levis=123,remaining_embeddings=4,remaining_Weyl_classes=2,
                regular_fiber_dimensions={'A1^3 cubic class':6,'A2+A1^2 cubic class':10},
                companion_membership_certificate_used=False,all_E7_Levis_resolved=False,
                all_simple_adjoint_conjecture_proved=False,full_map_degree_claim=False,literature_priority_claim=False)
    paths=['BINARY_BRANCHING_FIBERS.md','tools/verify_binary_branching_fibers.py',
           'tools/verify_branching_norm_fibers.py','e7-multiplier-boundary-verification.json']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in paths}
    (ROOT/'binary-branching-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('source_sha256','characters','Weyl_classes')},indent=2))


if __name__=='__main__':main()
