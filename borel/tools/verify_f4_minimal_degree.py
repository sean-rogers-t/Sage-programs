"""F4 minimal degree: independent weights, determinant volumes and fibers."""
from collections import Counter
from hashlib import sha256
from itertools import product, combinations
from pathlib import Path
import json
from sage.all import (QQ, ZZ, GF, matrix, vector, prod, factorial,
                      Polyhedron, PolynomialRing, WeylCharacterRing, RootSystem)
from sage.env import SAGE_VERSION
from monomial_staircase import staircase_length

ROOT=Path(__file__).resolve().parents[1]


def weights_and_geometry():
    E=list(matrix.identity(QQ,4).rows())
    explicit=[sign*v for v in E for sign in [-1,1]]
    explicit += [vector(QQ,s)/2 for s in product([-1,1],repeat=4)]
    # A separate highest-weight character computation includes the zero weights.
    character=WeylCharacterRing('F4',style='coroots')(0,0,0,1)
    multiplicities={tuple(w.to_vector()):int(m) for w,m in character.weight_multiplicities().items()}
    assert character.degree()==26 and multiplicities.pop((0,0,0,0))==2
    assert multiplicities=={tuple(w):1 for w in explicit}
    C=RootSystem(['F',4]).cartan_type().cartan_matrix()
    assert C.det()==1
    # Actual lattice basis e1,e2,e3,delta, with delta=(1/2,...,1/2).
    basis=matrix(QQ,E[:3]+[vector(QQ,[1,1,1,1])/2]).transpose()
    assert abs(basis.det())==QQ(1)/2
    W=[tuple(basis.inverse()*w) for w in explicit]
    assert all(a in ZZ for w in W for a in w)
    assert [d for d in matrix(ZZ,W).elementary_divisors() if d]==[1]*4
    moment=sum((w.column()*w.row() for w in explicit),matrix(QQ,4,4))
    assert moment==6*matrix.identity(QQ,4)
    # A noncentral involution has zero image; finiteness of the torus
    # restriction must not be transferred to the full group morphism.
    evaluate=lambda w:QQ(-1)**ZZ(sum((basis.inverse()*w)[:3]))
    assert sum((evaluate(w)*w for w in explicit),vector(QQ,[0]*4))==0
    long_roots=[s*E[i]+t*E[j] for i,j in combinations(range(4),2) for s,t in product([-1,1],repeat=2)]
    centralizing=sum(evaluate(w)==1 for w in explicit+long_roots)
    assert centralizing==32
    # Independent polytope in the fundamental-weight basis has the same volume.
    P=RootSystem(['F',4]).weight_lattice()
    orbit=[tuple(w.to_vector()) for w in P.fundamental_weights()[4].orbit()]
    Q=Polyhedron(vertices=W); independent=Polyhedron(vertices=orbit)
    assert Q.interior_contains(vector(QQ,[0]*4))
    volumes=[]; heights=[]
    for face in Q.facets():
        V=[vector(QQ,v) for v in face.vertices()]
        assert len(V)==6
        center=sum(V)/6; remaining={tuple(v) for v in V}; pairs=[]
        while remaining:
            a=vector(QQ,min(remaining)); b=2*center-a
            assert tuple(b) in remaining and a!=b
            remaining.remove(tuple(a)); remaining.remove(tuple(b)); pairs.append((a,b))
        assert len(pairs)==3 and matrix([a-center for a,b in pairs]).rank()==3
        determinants=[abs(matrix(list(pairs[0])+list(c)).det()) for c in product(*pairs[1:])]
        assert determinants==[1]*4
        volumes.append(int(sum(determinants)))
        h=face.ambient_Hrepresentation()[0]
        from sage.all import gcd
        assert all(a in ZZ for a in h) and gcd(list(h)[1:])==1
        heights.append(int(h[0]))
    assert len(volumes)==24 and set(heights)=={1}
    assert sum(volumes)==factorial(4)*Q.volume()==factorial(4)*independent.volume()==96
    for d in range(3): assert all(len(f.vertices())==d+1 for f in Q.faces(d))
    critical_values=sorted({2*sum(s) for s in product([-1,1],repeat=3)})
    assert critical_values==[-6,-2,2,6]
    # Exact unit-ideal certificate for f=df=0 on the three-variable cover.
    R=PolynomialRing(QQ,['x','y','z']); x,y,z=R.gens()
    assert R.ideal(x*x-1,y*y-1,z*z-1,x+y+z)==R.ideal(1)
    # Unequal coefficients show why odd crosspolytope geometry alone is insufficient.
    assert 1+1+2*(-1)==0
    return W,dict(dimension=26,nonzero_weights=24,zero_weight_multiplicity=2,
                  actual_lattice='Q(F4)=P(F4)=Z^4+Z*(1/2,1/2,1/2,1/2)',
                  euclidean_lattice_covolume='1/2',weight_lattice_index=1,
                  central_cover_degree=1,trace_moment_matrix='6*I4',
                  independent_highest_weight_character_agrees=True,
                  facet_count=24,facet_type='three-dimensional crosspolytope',
                  primitive_height=1,normalized_cone_volumes=volumes,
                  f_vector=list(map(int,Q.f_vector())),
                  normalized_volume=96,independent_fundamental_basis_volume=96,
                  cross_cover_critical_values=critical_values,
                  face_singular_zero_ideal_is_unit=True,boundary_correction=0,
                  noncentral_zero_witness=dict(epsilon_characters=[-1]*4,
                    delta_character=1,centralizing_roots=centralizing,
                    conjugacy_class_dimension=52-(4+centralizing)),
                  group_torus_and_quotient_degree=96)


def fiber(prime,model,W):
    A=GF(prime); target=[1,2,3,5]
    if model=='weight_lattice':
        R=PolynomialRing(A,['z1','z2','z3','b','v'],order='degrevlex')
        z=R.gens()[:4]; inv=R.gen(4)
        minimum=[min(0,min(w[i] for w in W)) for i in range(4)]
        den=prod(t**(-m) for t,m in zip(z,minimum))
        mon=[prod(t**(a-m) for t,a,m in zip(z,w,minimum)) for w in W]
        # Gradient coordinates in M, obtained from the epsilon target 6*y.
        btarget=[6*(target[i]-target[3]) for i in range(3)]+[12*target[3]]
        equations=[sum(A(w[i])*u for w,u in zip(W,mon))-btarget[i]*den for i in range(4)]
        equations.append(inv*prod(z)-1)
        ell=z[3]+2*z[0]+3*z[1]
    else:
        R=PolynomialRing(A,['z1','z2','z3','z4','b','v'],order='degrevlex')
        z=R.gens()[:4]; b,inv=R.gens()[4:]; D=prod(z); F=prod(1+t for t in z)
        equations=[2*b*(t*t-1)+t*(2*t*F.derivative(t)-F)-12*y*b*t for t,y in zip(z,target)]
        equations += [b*b-D,inv*b*D-1]
        ell=b+2*z[0]+3*z[1]
    I=R.ideal(equations); gb=I.groebner_basis(algorithm='libsingular:slimgb')
    leading=[tuple(g.lm().exponents()[0]) for g in gb]
    length=staircase_length(leading); normal=list(I.normal_basis())
    assert length==len(normal)==96
    positions={tuple(v.exponents()[0]):i for i,v in enumerate(normal)}
    entries={}
    for j,v in enumerate(normal):
        for e,c in (ell*v).reduce(gb).dict().items(): entries[positions[tuple(e)],j]=c
    operator=matrix(A,length,length,entries); p=operator.charpoly('q')
    assert p.gcd(p.derivative()).degree()==0
    return dict(prime=prime,model=model,epsilon_normalized_target=target,
                saturated_fiber_length=length,reduced_by_squarefree_multiplication=True,
                linear_form='b+2*z1+3*z2',
                characteristic_polynomial_coefficients=[int(a) for a in p.list()],
                leading_monomials=[list(map(int,e)) for e in leading])


def main():
    W,geometry=weights_and_geometry(); fibers=[]
    for prime in [32003,32009]:
        rows=[fiber(prime,model,W) for model in ['weight_lattice','product_character']]
        assert rows[0]['characteristic_polynomial_coefficients']==rows[1]['characteristic_polynomial_coefficients']
        fibers+=rows; print('Verified both F4 fiber models at',prime,flush=True)
    report=dict(result='passed',date='2026-10-22',sage_version=SAGE_VERSION,
                geometry=geometry,independent_fibers=fibers,
                same_target_multiplication_polynomials_agree=True,
                characteristic_zero_degree_from_boundary_proof=True,
                finite_flat_torus_and_quotient_written_corollary=True,
                full_group_map_not_finite=True,
                finite_fibers_used_as_proof_of_generic_degree=False,
                novelty_claim=False,externally_reviewed=False,runtime_API_changed=False,
                scope='Characteristic-zero written boundary proof; independent weight, lattice and determinant checks plus four reduced specialized fibers in two constructions, sharing Sage/Singular arithmetic.')
    paths=['F4_MINIMAL_DEGREE.md','tools/verify_f4_minimal_degree.py','tools/monomial_staircase.py']
    report['source_sha256']={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths}
    (ROOT/'f4-minimal-degree-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['source_sha256','independent_fibers']},indent=2))


if __name__=='__main__': main()
