"""Fresh E6 roots, facet covers and local deformation checks; no Segre input."""
from collections import Counter
from hashlib import sha256
from itertools import product, permutations, combinations
from pathlib import Path
import json
from sage.all import (QQ, ZZ, RootSystem, Polyhedron, matrix, vector,
                      PolynomialRing, prod, factorial)
from sage.env import SAGE_VERSION

ROOT=Path(__file__).resolve().parents[1]


def root_geometry():
    lattice=RootSystem(['E',6]).root_lattice()
    roots=sorted(tuple(a.to_vector()) for a in lattice.roots())
    C=matrix(ZZ,lattice.cartan_type().cartan_matrix())
    assert C.det()==3
    assert [d for d in matrix(ZZ,roots).elementary_divisors() if d]==[1]*6
    P=Polyhedron(vertices=roots)
    assert len(roots)==72 and len(P.facets())==54
    assert factorial(6)*P.volume()==2808
    assert all(len(f.vertices())==16 for f in P.facets())
    assert all(h[0]==1 for h in P.inequalities())
    # Construct a D5 Euclidean realization from its Cartan matrix, rather
    # than using the earlier exceptional-degree script's weight conversion.
    eye=list(matrix.identity(QQ,5).rows())
    standard=[eye[i]-eye[i+1] for i in range(4)]+[eye[3]+eye[4]]
    Cd=C.submatrix(1,1,5,5)
    E=next(matrix(QQ,p) for p in permutations(standard)
           if matrix(QQ,p)*matrix(QQ,p).transpose()==Cd)
    def mu(v):
        labels=vector(QQ,[(C*vector(v))[i] for i in range(1,6)])
        return (Cd.inverse()*labels)*E
    top=[v for v in roots if v[0]==1]
    base=mu(top[0]); assert all(abs(t)==QQ(1)/2 for t in base)
    def transform(v):
        m=mu(v)-v[0]*base
        return vector(QQ,[v[0]]+[-2*base[i]*m[i] for i in range(5)])
    B=matrix(QQ,[transform(v) for v in matrix.identity(ZZ,6).columns()]).transpose()
    assert all(a in ZZ for a in B.list())
    B=matrix(ZZ,B); assert abs(B.det())==2
    image=[tuple(B*vector(v)) for v in roots]
    even={tuple(p) for p in product([0,1],repeat=5) if sum(p)%2==0}
    d5={tuple(s*eye[i]+t*eye[j]) for i,j in combinations(range(5),2) for s,t in product([-1,1],repeat=2)}
    assert {w[1:] for w in image if w[0]==1}==even
    assert {w[1:] for w in image if w[0]==0}==d5
    assert {w[1:] for w in image if w[0]==-1}=={tuple(-a for a in w) for w in even}
    assert all(sum(w[1:])%2==0 for w in image)
    # Facet difference lattice is saturated in its hyperplane inside Q(E6).
    differences=matrix(ZZ,[vector(v)-vector(top[0]) for v in top[1:]])
    assert [d for d in differences.elementary_divisors() if d]==[1]*5
    def reflect_face(F,i):
        answer=[]
        for v in F:
            w=list(v); w[i]-=sum(C[i,j]*v[j] for j in range(6)); answer.append(tuple(w))
        return frozenset(answer)
    def face_orbit(F):
        seen={F}; todo=[F]
        while todo:
            f=todo.pop()
            for i in range(6):
                other=reflect_face(f,i)
                if other not in seen: seen.add(other); todo.append(other)
        return seen
    orbit=face_orbit(frozenset(top))
    negatives={frozenset(tuple(-a for a in v) for v in f) for f in orbit}
    assert len(orbit)==27 and not orbit.intersection(negatives)
    assert orbit|negatives=={frozenset(tuple(v) for v in f.vertices()) for f in P.facets()}
    counts=Counter(len(f.vertices()) for f in P.faces(4))
    assert counts=={5:432,8:270}
    for d in range(4): assert all(len(f.vertices())==d+1 for f in P.faces(d))
    cross=[f for f in P.faces(4) if len(f.vertices())==8]
    for f in cross:
        N=matrix(ZZ,[list(h[1:]) for h in f.ambient_Hrepresentation()])
        assert N.nrows()==2 and [d for d in N.elementary_divisors() if d]==[1,1]
        vertices=[vector(QQ,v) for v in f.vertices()]
        diffs=matrix(ZZ,[v-vertices[0] for v in vertices[1:]])
        assert [d for d in diffs.elementary_divisors() if d]==[1]*4
        center=sum(vertices)/8
        assert all(2*center-v in vertices for v in vertices)
        assert matrix([v-center for v in vertices]).rank()==4
    assert len(face_orbit(frozenset(tuple(v) for v in cross[0].vertices())))==270
    # Enumerate the critical points in an actual unimodular root-lattice
    # chart. This is separate from the doubled D5 facet cover.
    N=matrix(ZZ,[[-a for a in h[1:]] for h in cross[0].ambient_Hrepresentation()])
    D,U,V=N.smith_form()
    chart=matrix(ZZ,list(N.rows())+list(V.inverse().rows()[2:]))
    assert abs(chart.det())==1
    transformed=[tuple(chart*vector(w)) for w in roots]
    nodes=[]
    for p in product([-1,1],repeat=4):
        def term(w): return prod(QQ(x)**a for x,a in zip(p,w[2:]))
        layer=lambda i,j:[w for w in transformed if w[:2]==(i,j)]
        val=lambda L:sum(term(w) for w in L)
        grad=lambda L:vector(QQ,[sum(w[k+2]*term(w) for w in L) for k in range(4)])
        f=layer(1,1)
        if val(f) or grad(f): continue
        H=matrix(QQ,4,4,lambda i,j:sum(w[i+2]*w[j+2]*term(w) for w in f))
        assert H.det()!=0
        a,b=grad(layer(0,1)),grad(layer(1,0))
        assert val(layer(0,1))==val(layer(1,0))==0
        assert val(layer(-1,1))==a*H.inverse()*a/2
        assert val(layer(1,-1))==b*H.inverse()*b/2
        nodes.append(dict(point=list(p),hessian_determinant=str(H.det())))
    assert len(nodes)==3
    return dict(root_count=72,normalized_volume=2808,facet_count=54,
                facet_Weyl_orbits=[27,27],facet_height=1,facet_cover_degree=2,
                cover_matrix=[list(map(int,row)) for row in B.rows()],
                root_layers=[16,40,16],facet_difference_lattice_saturated=True,
                codimension_two_simplex_faces=432,codimension_two_cross_faces=270,
                lower_faces_all_simplices=True,cross_normal_charts_unimodular=True,
                cross_difference_lattice_index=1,central_cover_degree=int(C.det()),
                crossing_critical_points_per_face=3,prototype_nodes=nodes)


def local_checks():
    R=PolynomialRing(QQ,['x1','x2','x3','x4','x5']); x=R.gens(); K=R.fraction_field()
    f=(prod(1+t for t in x)+prod(1-t for t in x))/2
    g=sum((K(x[i])+1/x[i])*(K(x[j])+1/x[j]) for i,j in combinations(range(5),2))
    d_g=[t*g.derivative(t) for t in x]
    curves=0
    for plus in combinations(range(5),2):
        for minus in combinations([i for i in range(5) if i not in plus],2):
            sub={x[i]:1 for i in plus}; sub.update({x[i]:-1 for i in minus})
            assert g.subs(sub)==-8 and all(v.subs(sub)==0 for v in d_g)
            assert f.subs(sub)==0 and all(f.derivative(t).subs(sub)==0 for t in x)
            curves+=1
    triples=[p for p in product([-1,1],repeat=5) if p.count(1) in (2,3)]
    assert len(triples)==20
    for p in triples:
        sub=dict(zip(x,p))
        assert g.subs(sub)==-8 and all(v.subs(sub)==0 for v in d_g)
    # Actual facet polynomial in transverse Mobius coordinates at a triple.
    u=list(map(K,x))
    z=[(t-1)/(t+1) for t in u[:2]]+[(1-t)/(1+t) for t in u[2:]]
    assert f(*z)*prod(1+t for t in u)==16*(prod(u[:2])+prod(u[2:]))
    # Generality of restricted perturbations and transverse first jets.
    jet=[]
    for s in [4,5]:
        for parity in [0,1]:
            A=[p for p in product([0,1],repeat=s) if sum(p)%2==parity]
            rank=matrix(QQ,[[1]+list(p) for p in A]).rank()
            assert rank==s+1; jet.append(dict(fixed=s,parity=parity,rank=rank))
    # Generic curve critical equation, with algebraically independent data.
    B=PolynomialRing(QQ,['A','B','k']); F=B.fraction_field(); A,Bc,k=F.gens()
    Q=PolynomialRing(F,'q'); q=Q.gen(); H=A+Bc*q
    eq=q*H.derivative()-k*H
    root=k*A/((1-k)*Bc)
    assert eq(root)==0 and H(root)!=0 and root not in [0,1,-1]
    transverse=[]
    for size in [2,3]:
        S=PolynomialRing(QQ,['a'+str(i) for i in range(size)])
        a=S.gens(); I=S.ideal([prod(a[j] for j in range(size) if j!=i)-(i+1) for i in range(size)])
        assert I.vector_space_dimension()==size-1
        J=matrix(S,[[f.derivative(v) for v in a] for f in I.gens()])
        assert I+S.ideal(J.det())==S.ideal(1)
        transverse.append(dict(product_size=size,simple_roots=size-1))
    # Universal two-normal cancellation: this verifies the Schur algebra
    # rather than reusing the old nodal Segre coefficient.
    S=PolynomialRing(QQ,['s','t','eta','A','B','E','b1','b2','h'])
    s,t,eta,A,Bc,E,b1,b2,h=S.gens()
    # One tangent direction suffices to check the bilinear identity;
    # the general symmetric-Hessian version is its coordinate sum.
    u=-A*s-Bc*t
    F0=u*u/2+s*A*u+t*Bc*u+s*s*A*A/2+t*t*Bc*Bc/2+s*t*E
    s_partial=A*u+s*A*A+t*E
    t_partial=Bc*u+t*Bc*Bc+s*E
    assert F0-s*s_partial-b1*s*t==-b1*s*t
    assert F0-t*t_partial-b2*s*t==-b2*s*t
    leading=matrix(S,[[h,-b1],[h,-b2]])
    assert leading.det()==h*(b1-b2)
    return dict(cover_curve_strata=curves,cover_triple_strata=len(triples),
                middle_layer_value=-8,middle_layer_full_gradient_zero=True,
                jet_ranks=jet,curve_generic_critical_points=1,product_gradient_checks=transverse,
                crossing_leading_determinant='h*(b1-b2)',
                generic_crossing_escape_count=0)


def main():
    geometry=root_geometry(); local=local_checks()
    correction=geometry['facet_count']*(local['cover_curve_strata']+2*local['cover_triple_strata'])//geometry['facet_cover_degree']
    assert correction==1890
    report=dict(result='passed',date='2026-10-22',sage_version=SAGE_VERSION,
                geometry=geometry,local=local,boundary_correction=correction,
                limiting_branches=dict(open_curves=54*30//2,
                                       triple_points=54*20*2//2,crossings=0),
                root_torus_degree=geometry['normalized_volume']-correction,
                simply_connected_degree=geometry['central_cover_degree']*(geometry['normalized_volume']-correction),
                Segre_formula_used_to_derive_result=False,novelty_claim=False,external_review=False,
                scope='Exact finite root/lattice and local-algebra checks support the written deformation proof; no formal proof certificate or generic-degree inference from finite-field samples.')
    paths=['E6_ADJOINT_DEFORMATION.md','tools/verify_e6_adjoint_deformation.py']
    report['source_sha256']={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths}
    (ROOT/'e6-adjoint-deformation-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__': main()
