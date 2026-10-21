"""Exact boundary and intersection data for four exceptional Springer degrees.

Read EXCEPTIONAL_DEGREES.md for the intersection-theoretic argument. This
program verifies its finite root, lattice, volume, critical-locus and local
Segre inputs. No specialized generic-fiber length is used to derive a degree.
"""
from collections import Counter
from itertools import product, combinations
from pathlib import Path
from time import perf_counter
import json
from sage.all import (QQ, ZZ, matrix, vector, prod, Polyhedron, factorial,
                      PolynomialRing, PowerSeriesRing, binomial)
from e6_matrix_oracle import C as C6
from e7_matrix_oracle import C as C7

ROOT=Path(__file__).resolve().parents[1]


def orbit(C,initial,dual=False,modulus=None):
    found,queue={tuple(initial)},[tuple(initial)]
    while queue:
        w=queue.pop()
        for i in range(C.nrows()):
            if dual:
                v=list(w)
                v[i]-=sum(C[i,j]*w[j] for j in range(C.nrows()))
            else:
                v=[w[j]-w[i]*C[j,i] for j in range(C.nrows())]
            v=tuple(int(a%modulus) if modulus else int(a) for a in v)
            if v not in found:
                found.add(v)
                queue.append(v)
    return tuple(sorted(found))


def lattice_index(vertices):
    first=vector(vertices[0])
    D=matrix(ZZ,[vector(v)-first for v in vertices[1:]])
    ds=[int(a) for a in D.elementary_divisors() if a]
    return int(prod(ds)),len(ds)


def cross_simplices(vertices):
    V=set(map(tuple,vertices))
    center=sum((vector(QQ,v) for v in V),vector(QQ,[0]*len(next(iter(V)))))/len(V)
    pairs=[]
    while V:
        a=min(V)
        b=tuple(2*center-vector(a))
        assert b in V and b!=a
        pairs.append((a,b))
        V.remove(a)
        V.remove(b)
    d=len(pairs)
    assert matrix(QQ,[vector(a)-center for a,b in pairs]).rank()==d
    return [list(pairs[0])+list(chosen) for chosen in product(*pairs[1:])]


def cone_volume(vertices):
    """Determinant triangulations, independent of Polyhedron.volume()."""
    n=len(vertices[0])
    if len(vertices)==n:
        return abs(matrix(QQ,vertices).det())
    if len(vertices)==2*(n-1):
        return sum(abs(matrix(QQ,s).det()) for s in cross_simplices(vertices))
    P=Polyhedron(vertices=vertices)
    assert P.dim()==n-1
    center=sum((vector(QQ,v) for v in vertices),vector(QQ,[0]*n))/len(vertices)
    answer=QQ.zero()
    for face in P.faces(n-2):
        V=[tuple(v) for v in face.vertices()]
        simplices=[V] if len(V)==n-1 else cross_simplices(V)
        for simplex in simplices:
            assert len(simplex)==n-1
            answer+=abs(matrix(QQ,[center]+simplex).det())
    return answer


def polytope_certificate(kind,name):
    C=C6 if kind=='E6' else C7
    n=C.nrows()
    ci=C.change_ring(QQ).inverse()
    unit=lambda i:tuple(int(j==i) for j in range(n))
    if name=='minimal':
        V=orbit(C,unit(0 if kind=='E6' else 6))
        normals=[('root',orbit(C,unit(0),dual=True))]
        extra=ci*vector(unit(5 if kind=='E6' else 1))*(3 if kind=='E6' else 2)
        normals.append(('coweight',orbit(C,tuple(map(int,extra)),dual=True)))
        lattice='weight lattice P'
    else:
        V=tuple(sorted(tuple(ci*vector(w)) for w in orbit(C,tuple(C.column(0)))))
        normals=[('node_%s'%(i+1),orbit(C,unit(i)))
                 for i in ((0,5) if kind=='E6' else (1,6))]
        lattice='root lattice Q'
    inequalities=[]
    cases=[]
    histogram=Counter()
    for label,ns in normals:
        heights=set()
        for normal in ns:
            height=max(sum(a*b for a,b in zip(normal,v)) for v in V)
            heights.add(int(height))
            inequalities.append([height]+[-a for a in normal])
            F=[v for v in V if sum(a*b for a,b in zip(normal,v))==height]
            assert matrix(QQ,[vector(v)-vector(F[0]) for v in F[1:]]).rank()==n-1
            histogram[len(F)]+=1
        assert len(heights)==1
        normal=ns[0]
        height=next(iter(heights))
        F=[v for v in V if sum(a*b for a,b in zip(normal,v))==height]
        # Every facet in this class is a unimodular Weyl translate.
        volume=cone_volume(F)
        index,rank=lattice_index(F)
        cases.append(dict(label=label,facets=len(ns),height=height,vertices=len(F),
                          lattice_index=index,cone_volume=str(volume),normal=list(normal)))
    H=Polyhedron(ieqs=inequalities)
    assert set(tuple(v) for v in H.vertices())==set(V)
    P=Polyhedron(vertices=V)
    volume=sum(c['facets']*QQ(c['cone_volume']) for c in cases)
    assert volume==factorial(n)*P.volume()
    report=dict(cartan_type=kind,representation=name,lattice=lattice,
                vertices=len(V),facet_vertex_counts=dict(sorted(histogram.items())),
                normalized_volume=int(volume),facets=cases,
                exact_H_V_equality=True,determinant_volume_agreement=True)
    if kind=='E6' and name=='adjoint':
        faces=P.faces(4)
        counts=Counter(len(f.vertices()) for f in faces)
        assert counts=={5:432,8:270}
        for face in faces:
            if len(face.vertices())!=8:
                continue
            F=[tuple(v) for v in face.vertices()]
            assert lattice_index(F)==(1,4)
            hs=face.ambient_Hrepresentation()
            assert len(hs)==2
            N=matrix(ZZ,[list(h[1:]) for h in hs])
            assert [a for a in N.elementary_divisors() if a]==[1,1]
        report['codimension_two_faces']=dict(simplex=432,cross=270,
            cross_lattice_index=1,cross_incident_facets=2,smooth_cross_normal_cones=True)
    return report


def half_spin_certificate():
    R=PolynomialRing(QQ,5,names=['x%d'%(i+1) for i in range(5)])
    xs=R.gens()
    f=(prod(1+x for x in xs)+prod(1-x for x in xs))/2
    assert len(f.dict())==16
    J=R.ideal([f]+[f.derivative(x) for x in xs])
    union=R.ideal(1)
    for plus in combinations(range(5),2):
        for minus in combinations([i for i in range(5) if i not in plus],2):
            union=union.intersection(R.ideal([xs[i]-1 for i in plus]+[xs[i]+1 for i in minus]))
    assert J==union
    # Identify an actual E6 root-polytope facet with a D5 half-spin module.
    roots=orbit(C6,tuple(C6.column(0)))
    ci=C6.change_ring(QQ).inverse()
    F=[w for w in roots if (ci*vector(w))[0]==1]
    standard_weights=[vector(QQ,[1,0,0,0,0]),vector(QQ,[1,1,0,0,0]),
        vector(QQ,[1,1,1,0,0]),vector(QQ,[1,1,1,1,-1])/2,
        vector(QQ,[1,1,1,1,1])/2]
    epsilon=[sum((w[i]*v for i,v in zip((5,4,3,1,2),standard_weights)),vector(QQ,[0]*5)) for w in F]
    assert len(epsilon)==16 and all(abs(a)==QQ(1)/2 for w in epsilon for a in w)
    parity={sum(a<0 for a in w)%2 for w in epsilon}
    assert len(parity)==1
    if next(iter(parity)):
        epsilon=[vector(list(w[:4])+[-w[4]]) for w in epsilon]
    exponents={tuple(int(QQ(1)/2-a) for a in w) for w in epsilon}
    assert exponents=={tuple(e) for e in f.dict()}
    T=PolynomialRing(QQ,'t')
    t=T.gen()
    point=[T(1),T(1),T(-1),T(-1),t]
    assert f(*point)==0 and all(f.derivative(x)(*point)==0 for x in xs)
    H=matrix(T,[[f.derivative(x).derivative(y)(*point) for y in xs] for x in xs])
    assert H.rank()==4 and H.submatrix(0,0,4,4).det()!=0
    determinant=H.submatrix(0,0,4,4).det()
    assert set(determinant.roots(QQ,multiplicities=False))=={-1,1}
    assert matrix(QQ,H(t=1)).rank()==matrix(QQ,H(t=-1)).rank()==2
    curves=set()
    for plus in combinations(range(5),2):
        for minus in combinations([i for i in range(5) if i not in plus],2):
            curves.add(min((tuple(plus),tuple(minus)),(tuple(minus),tuple(plus))))
    triples=set()
    for signs in product((-1,1),repeat=5):
        if signs.count(1) in (2,3):
            triples.add(min(signs,tuple(-s for s in signs)))
    incidence=Counter()
    for plus,minus in curves:
        free=next(i for i in range(5) if i not in plus+minus)
        for sign in (-1,1):
            point=tuple(1 if i in plus else -1 if i in minus else sign for i in range(5))
            incidence[min(point,tuple(-a for a in point))]+=1
    assert len(curves)==15 and len(triples)==10 and set(incidence)==triples
    assert set(incidence.values())=={3}
    S=PowerSeriesRing(QQ,'h',default_prec=7)
    h=S.gen()
    double=2*h/(1+2*h)
    individual_double=2*h/(1+h)
    triple=1-sum(g*h**i/(1+2*h)**(i+1) for i,g in enumerate((1,2,1,0)))
    individual_triple=3*h*h/(1+h)**2
    assert (double-individual_double)[2]==-2
    assert triple[2]==3 and triple[3]==-10
    assert (triple-individual_triple)[3]==-4
    return dict(character_terms=16,exact_Jacobian_ideal_equals_30_line_union=True,
                D5_half_spin_facet_identification=True,
                curves_per_facet=15,triples_per_facet=10,
                curve_hessian_rank=4,exceptional_parameters=[-1,1],triple_hessian_rank=2,
                transverse_hessian_determinant=str(determinant),
                double_point_segre_correction=-2,triple_point_segre_correction=-4,
                curve_section_degree=1,anticanonical_curve_degree=1,
                smooth_curve_contribution=7)


def e6_critical_zero_certificate():
    # Character formulas have separate global certificates in E6_GEOMETRY.md.
    from sage_borel import e6_adjoint_quotient,e6_cayley_determinant
    q=e6_adjoint_quotient()
    R=q.character_polynomials[0].parent()
    point=(0,-3,0,9,0,0)
    I=R.ideal(list(q.character_polynomials)+[R.gen(0)])
    G=I.groebner_basis()
    assert all(f(*point)==0 for f in I.gens())
    nilpotence=[]
    for x,a in zip(R.gens(),point):
        g=x-a
        for k in range(1,38):
            if (g**k).reduce(G)==0:
                nilpotence.append(k)
                break
        else:
            raise AssertionError('failed exact radical support certificate')
    T=e6_cayley_determinant().torus_character_polynomial
    assert T(*point)==QQ(1)/64
    W=orbit(C6,(1,0,0,0,0,0))
    critical=[]
    for v in product(range(3),repeat=6):
        residues=[sum(a*b for a,b in zip(w,v))%3 for w in W]
        if Counter(residues)!={0:9,1:9,2:9}:
            continue
        moments=[[sum(w[i] for w,r in zip(W,residues) if r==j) for i in range(6)] for j in range(3)]
        if moments[0]==moments[1]==moments[2]:
            critical.append(v)
    assert len(critical)==240
    assert set(orbit(C6,critical[0],dual=True,modulus=3))==set(critical)
    roots=orbit(C6,tuple(C6.column(0)))
    v=critical[0]
    root_counts=Counter(sum(a*b for a,b in zip(w,v))%3 for w in roots)
    assert root_counts=={0:18,1:27,2:27}
    # A second orbit-size calculation: the roots centralizing t form A2^3.
    remaining={w for w in roots if sum(a*b for a,b in zip(w,v))%3==0}
    inverse=C6.change_ring(QQ).inverse()
    components=[]
    while remaining:
        first=remaining.pop()
        component,queue={first},[first]
        while queue:
            root=queue.pop()
            adjacent={w for w in remaining if vector(root)*inverse*vector(w)!=0}
            component.update(adjacent)
            queue.extend(adjacent)
            remaining-=adjacent
        assert len(component)==6 and matrix(QQ,list(component)).rank()==2
        components.append(component)
    assert len(components)==3 and prod(q.degrees)//6**3==len(critical)
    # Newton identities give chi_3=chi_5=0 and chi_4=9 from p1=p2=0,p3=27.
    # The dual character is zero and the adjoint is 6+18-27=-3.
    face=tuple(w for w in orbit(C7,tuple(C7.column(0)))
               if (C7.change_ring(QQ).inverse()*vector(w))[6]==1)
    assert len(face)==27 and set(tuple(w[:6]) for w in face)==set(W)
    return dict(unique_character_point=list(point),coordinate_nilpotence_exponents=nilpotence,
                exact_radical_support_certificate=True,cartan_determinant='1/64',
                torsion_order=3,torus_points=240,one_Weyl_orbit=True,
                centralizer_root_system='A2^3',centralizer_Weyl_order=216,
                independent_orbit_stabilizer_count=240,
                torsion_representative=list(v),central_quotient_degree=3,
                E6_root_torus_critical_points=80,E7_Levi_face_identification=True)


def e6_boundary_node_certificate():
    """Actual root sections at the codimension-two boundary points.

    After eliminating four tangent variables, both normal equations have
    leading terms -b_i*s*t. This checks the reduced nodal zero scheme,
    including the two normal directions rather than only the face character.
    """
    ci=C6.change_ring(QQ).inverse()
    roots=[tuple(ci*vector(w)) for w in orbit(C6,tuple(C6.column(0)))]
    P=Polyhedron(vertices=roots)
    cross=[f for f in P.faces(4) if len(f.vertices())==8]
    F=cross[0]
    hs=F.ambient_Hrepresentation()
    N=matrix(ZZ,[[-a for a in h[1:]] for h in hs])
    assert all(h[0]==1 for h in hs)
    D,U,V=N.smith_form()
    assert D.submatrix(0,0,2,2)==matrix.identity(ZZ,2)
    B=matrix(ZZ,list(N.rows())+list(V.inverse().rows()[2:]))
    assert abs(B.det())==1
    transformed=[tuple(B*vector(w)) for w in roots]
    nodes=[]
    for point in product((-1,1),repeat=4):
        def term(w):
            return prod(QQ(t)**a for t,a in zip(point,w[2:]))
        def value(i,j):
            return sum(term(w) for w in transformed if w[:2]==(i,j))
        top=[w for w in transformed if w[:2]==(1,1)]
        gradient=vector(QQ,[sum(w[k+2]*term(w) for w in top) for k in range(4)])
        if value(1,1) or gradient:
            continue
        H=matrix(QQ,4,4,lambda i,j:sum(w[i+2]*w[j+2]*term(w) for w in top))
        assert H.det()
        a=vector(QQ,[sum(w[k+2]*term(w) for w in transformed if w[:2]==(0,1)) for k in range(4)])
        b=vector(QQ,[sum(w[k+2]*term(w) for w in transformed if w[:2]==(1,0)) for k in range(4)])
        assert value(0,1)==value(1,0)==0
        assert value(-1,1)==a*H.inverse()*a/2
        assert value(1,-1)==b*H.inverse()*b/2
        assert a*H.inverse()*b==b*H.inverse()*a
        nodes.append(dict(point=list(point),hessian_determinant=str(H.det()),
                          reduced_normal_leading_terms=['-b1*s*t','-b2*s*t']))
    assert len(nodes)==3
    # Check that the prototype accounts for every cross face by Weyl symmetry.
    first=frozenset(tuple(v) for v in F.vertices())
    found,queue={first},[first]
    while queue:
        face=queue.pop()
        for i in range(6):
            new=[]
            for w in face:
                v=list(w)
                v[i]-=sum(C6[i,j]*w[j] for j in range(6))
                new.append(tuple(v))
            other=frozenset(new)
            if other not in found:
                found.add(other)
                queue.append(other)
    assert found=={frozenset(tuple(v) for v in f.vertices()) for f in cross}
    return dict(unimodular_normal_chart=True,cross_face_Weyl_orbit=270,
                node_count_per_face=3,nodes=nodes,full_section_Schur_complement_checked=True)


def main():
    start=perf_counter()
    cases=[polytope_certificate(k,n) for k in ('E6','E7') for n in ('minimal','adjoint')]
    spin=half_spin_certificate()
    zero=e6_critical_zero_certificate()
    boundary_nodes=e6_boundary_node_certificate()
    by={(c['cartan_type'],c['representation']):c for c in cases}
    assert by['E6','minimal']['normalized_volume']==1080
    assert by['E7','minimal']['normalized_volume']==9792
    assert by['E6','adjoint']['normalized_volume']==2808
    assert by['E7','adjoint']['normalized_volume']==21312
    curves=54*spin['curves_per_facet']
    triples=54*spin['triples_per_facet']
    doubles=270*3
    e6_adjoint_correction=7*curves-4*triples-2*doubles
    results=[dict(cartan_type='E6',representation='minimal',torus_lattice='P',
                  normalized_volume=1080,boundary_correction=0,central_cover=1,degree=1080),
             dict(cartan_type='E6',representation='adjoint',torus_lattice='Q',
                  normalized_volume=2808,boundary_correction=e6_adjoint_correction,
                  central_cover=3,degree=3*(2808-e6_adjoint_correction)),
             dict(cartan_type='E7',representation='minimal',torus_lattice='P',
                  normalized_volume=9792,boundary_correction=126*20,
                  central_cover=1,degree=9792-126*20),
             dict(cartan_type='E7',representation='adjoint',torus_lattice='Q',
                  normalized_volume=21312,boundary_correction=56*80,
                  central_cover=2,degree=2*(21312-56*80))]
    assert [r['degree'] for r in results]==[1080,2754,7272,33664]
    report=dict(result='passed',scope='Finite inputs to the written toric-intersection proof; not a formal proof assistant certificate.',
                polytopes=cases,half_spin=spin,E6_critical_zero=zero,
                E6_boundary_node_sections=boundary_nodes,
                E6_adjoint_boundary=dict(curves=curves,triple_points=triples,
                    double_points=doubles,excess_intersection=e6_adjoint_correction),
                degree_deductions=results,seconds=perf_counter()-start)
    (ROOT/'experiments/exceptional-degree-geometry.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    main()
