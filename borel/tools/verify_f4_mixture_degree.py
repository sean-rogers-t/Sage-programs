"""Exact F4 mixture geometry, normal jets and independently constructed fibers."""
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json
import argparse
from sage.all import (QQ, ZZ, GF, matrix, vector, prod, factorial, gcd,
                      Polyhedron, PolynomialRing, LaurentPolynomialRing,
                      WeylCharacterRing, RootSystem)
from sage.env import SAGE_VERSION
from monomial_staircase import staircase_length

ROOT=Path(__file__).resolve().parents[1]


def degree(a,b):
    assert a>=0 and b>=0 and a+b>0
    return 96 if b==0 else 240 if a==b else 336 if a==5*b else 384


def geometry():
    E=list(matrix.identity(QQ,4).rows())
    short=[s*v for v in E for s in [-1,1]]
    short += [vector(QQ,s)/2 for s in product([-1,1],repeat=4)]
    long=[s*E[i]+t*E[j] for i,j in combinations(range(4),2)
          for s,t in product([-1,1],repeat=2)]
    basis=matrix(QQ,E[:3]+[vector(QQ,[1]*4)/2]).transpose()
    to_M=lambda w:tuple(basis.inverse()*w)
    Wshort=list(map(to_M,short)); Wlong=list(map(to_M,long))
    chars=WeylCharacterRing('F4',style='coroots')
    for highest,weights,zero,dim in [((0,0,0,1),short,2,26),((1,0,0,0),short+long,4,52)]:
        ch=chars(*highest)
        mult={tuple(w.to_vector()):int(m) for w,m in ch.weight_multiplicities().items()}
        assert mult.pop((0,0,0,0))==zero and ch.degree()==dim
        assert mult=={tuple(w):1 for w in weights}
    assert sum((w.column()*w.row() for w in short),matrix(QQ,4,4))==6*matrix.identity(QQ,4)
    assert sum((w.column()*w.row() for w in long),matrix(QQ,4,4))==12*matrix.identity(QQ,4)
    assert abs(basis.det())==QQ(1)/2
    assert [d for d in matrix(ZZ,Wshort+Wlong).elementary_divisors() if d]==[1]*4
    P=Polyhedron(vertices=Wlong)
    assert P.interior_contains(vector(QQ,[0]*4))
    volumes=[]; centers=[]
    for face in P.facets():
        V=[vector(QQ,v) for v in face.vertices()]; assert len(V)==6
        center=sum(V)/6; centers.append(tuple(center)); assert tuple(center) in Wshort
        h=face.ambient_Hrepresentation()[0]
        assert h[0]==2 and gcd(list(h)[1:])==1
        on_face=[w for w in Wshort+Wlong if h[0]+sum(x*y for x,y in zip(h[1:],w))==0]
        assert set(on_face)=={tuple(v) for v in V}|{tuple(center)}
        remaining={tuple(v) for v in V}; pairs=[]
        while remaining:
            a=vector(QQ,min(remaining)); b=2*center-a
            assert tuple(b) in remaining and a!=b
            remaining.remove(tuple(a)); remaining.remove(tuple(b)); pairs.append((a,b))
        dets=[abs(matrix(list(pairs[0])+list(c)).det()) for c in product(*pairs[1:])]
        assert dets==[4]*4; volumes.append(int(sum(dets)))
    assert len(volumes)==len(set(centers))==24
    for d in range(3):
        for face in P.faces(d):
            assert len(face.vertices())==d+1
            assert not any(face.as_polyhedron().contains(w) for w in Wshort)
    lattice=RootSystem(['F',4]).weight_lattice()
    independent=Polyhedron(vertices=[tuple(w.to_vector()) for w in lattice.fundamental_weights()[1].orbit()])
    assert sum(volumes)==factorial(4)*P.volume()==factorial(4)*independent.volume()==384
    # A primitive normal chart: s=e^-delta; x,y,z=e^e2,e^e3,e^e4.
    # This is an integral basis, not an unrecorded finite cover.
    chart=matrix(QQ,[-vector(QQ,[1]*4)/2,E[1],E[2],E[3]]).transpose()
    assert abs(chart.det())==abs(basis.det())
    LR=LaurentPolynomialRing(QQ,['s','x','y','z']); s,x,y,z=LR.gens(); tang=[x,y,z]
    def mon(w): return prod(t**k for t,k in zip(LR.gens(),chart.inverse()*w))
    S=sum(mon(w) for w in short); L=sum(mon(w) for w in long)
    H=sum(t+1/t for t in tang); K=sum((u+1/u)*(v+1/v) for u,v in combinations(tang,2))
    A=prod(1+t for t in tang); D=prod(tang)
    assert S==s**-2/D+H+s*A+s**-1*A/D+s*s*D
    assert L==(s**-2/D+s*s*D)*H+K
    rows=[]
    for a,b in [(0,1),(1,1),(5,1),(2,1),(1,2),(2,2),(10,2)]:
        c=a+b; F=s*s*(c*S+b*L)
        f0=(c+b*H)/D; f1=c*A/D
        assert F==f0+s*f1+s*s*(c*H+b*K)+s**3*c*A+s**4*D*(c+b*H)
        points=[]
        for signs in product([-1,1],repeat=3):
            subs=dict(zip(tang,signs))
            if c+2*b*sum(signs): continue
            # Exact logarithmic Hessian and first intermediate-layer jet.
            H0=matrix(QQ,[[ (u*(v*f0.derivative(v)).derivative(u)).subs(subs) for v in tang] for u in tang])
            assert H0.det()!=0 and f0.subs(subs)==0
            assert all((u*f0.derivative(u)).subs(subs)==0 for u in tang)
            assert f1.subs(subs)==0
            assert all((u*f1.derivative(u)).subs(subs)==0 for u in tang)
            points.append(dict(point=list(signs),log_hessian_determinant=str(H0.det()),
                              first_layer_value=0,first_layer_gradient=[0,0,0],normal_length=2))
        assert len(points)==(3 if a==b else 1 if a==5*b else 0)
        correction=24*2*len(points)
        assert 384-correction==degree(a,b)
        rows.append(dict(a=a,b=b,facet_singular_points=points,total_boundary_length=correction,degree=degree(a,b)))
    # Elimination of the critical values is symbolic in r=c/b; no sampled ratio inference.
    R=PolynomialRing(QQ,'r'); r=R.gen()
    critical_product=prod(r+2*sum(e) for e in product([-1,1],repeat=3))
    assert critical_product==(r-6)*(r-2)**3*(r+2)**3*(r+6)
    return Wshort,Wlong,dict(actual_lattice='P(F4)=Q(F4)=Z^4+Z*delta',central_cover_degree=1,
        character='(a+b)*S_short+b*S_long+(2*a+4*b)',trace_moment='(6*a+18*b)*I4',
        independent_highest_weight_characters_agree=True,normalized_volume=384,
        independent_fundamental_basis_volume=384,facet_count=24,primitive_facet_height=2,
        facet_cone_volumes=volumes,lower_faces_simplex_without_extra_weights=True,
        primitive_normal_chart_verified=True,critical_value_product=str(critical_product),
        nonnegative_exceptional_ratios=['a=b','a=5*b'],boundary_checks=rows)


def fiber(a,b,prime,model,Wshort,Wlong):
    A=GF(prime); target=[1,2,3,5]; scale=6*a+18*b; c=a+b
    if model=='weight_lattice':
        W=Wshort+(Wlong if b else []); mult=[c]*len(Wshort)+([b]*len(Wlong) if b else [])
        R=PolynomialRing(A,['z1','z2','z3','u','v'],order='degrevlex')
        z=R.gens()[:4]; inv=R.gen(4)
        minimum=[min(0,min(w[i] for w in W)) for i in range(4)]
        den=prod(t**(-m) for t,m in zip(z,minimum))
        mon=[prod(t**(e-m) for t,e,m in zip(z,w,minimum)) for w in W]
        rhs=[scale*(target[i]-target[3]) for i in range(3)]+[2*scale*target[3]]
        equations=[sum(A(m*w[i])*u for w,m,u in zip(W,mult,mon))-rhs[i]*den for i in range(4)]
        equations.append(inv*prod(z)-1); ell=z[3]+2*z[0]+3*z[1]
    else:
        R=PolynomialRing(A,['z1','z2','z3','z4','u','v'],order='degrevlex')
        z=R.gens()[:4]; u,inv=R.gens()[4:]; D=prod(z); F=prod(1+t for t in z)
        equations=[]
        for i,t in enumerate(z):
            # Differentiate factored S_short and pair-product S_long in epsilon coordinates.
            others=[z[j] for j in range(4) if j!=i]; rest=prod(others)
            longsum=sum((q*q+1)*prod(v for v in others if v!=q) for q in others)
            eq=c*(2*u*(t*t-1)+t*(2*t*F.derivative(t)-F))*rest
            eq+=2*b*u*(t*t-1)*longsum-2*scale*target[i]*u*D
            equations.append(eq)
        equations += [u*u-D,inv*u*D-1]; ell=u+2*z[0]+3*z[1]
    I=R.ideal(equations); gb=I.groebner_basis(algorithm='libsingular:slimgb')
    leading=[tuple(g.lm().exponents()[0]) for g in gb]
    length=staircase_length(leading); normal=list(I.normal_basis())
    assert length==len(normal)==degree(a,b),(a,b,prime,model,length)
    positions={tuple(v.exponents()[0]):i for i,v in enumerate(normal)}; entries={}
    for j,v in enumerate(normal):
        for e,co in (ell*v).reduce(gb).dict().items(): entries[positions[tuple(e)],j]=co
    p=matrix(A,length,length,entries).charpoly('q')
    assert p.gcd(p.derivative()).degree()==0,(a,b,prime,model,'not reduced')
    return dict(a=a,b=b,prime=prime,model=model,normalized_epsilon_target=target,
        saturated_fiber_length=length,reduced_by_squarefree_multiplication=True,
        multiplication_function='u+2*z1+3*z2',characteristic_polynomial_coefficients=[int(v) for v in p.list()])


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--geometry-only',action='store_true')
    args=parser.parse_args(); Ws,Wl,geo=geometry()
    if args.geometry_only: print(json.dumps(geo,indent=2)); return
    fibers=[]
    for prime in [32003,32009]:
        for a,b in [(0,1),(1,1),(5,1),(2,1),(1,0)]:
            pair=[]
            for model in ['weight_lattice','product_character']:
                row=fiber(a,b,prime,model,Ws,Wl); pair.append(row)
                print('Verified',a,b,prime,model,'length',row['saturated_fiber_length'],flush=True)
            assert pair[0]['characteristic_polynomial_coefficients']==pair[1]['characteristic_polynomial_coefficients']
            fibers+=pair
    report=dict(result='passed',date='2026-10-22',sage_version=SAGE_VERSION,geometry=geo,
        independent_fibers=fibers,same_target_multiplication_polynomials_agree=True,
        formula={'b=0':96,'b>0,a=b':240,'b>0,a=5b':336,'b>0,other':384},
        characteristic_zero_proof='F4_MIXTURE_DEGREES.md',
        finite_fibers_used_as_proof_of_generic_degree=False,novelty_claim=False,externally_reviewed=False,
        runtime_API_changed=False,independence='Two equation constructions and geometry; shared Sage/Singular arithmetic.')
    paths=['F4_MIXTURE_DEGREES.md','tools/verify_f4_mixture_degree.py','tools/monomial_staircase.py']
    report['source_sha256']={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths}
    (ROOT/'f4-mixture-degree-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['independent_fibers','source_sha256']},indent=2))


if __name__=='__main__': main()
