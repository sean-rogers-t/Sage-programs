"""Audit the half-spin rank law and independently replay Spin_12 fibers.

Finite computations corroborate the written characteristic-zero deformation
argument; the numerical rank checks do not prove the all-rank theorem.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json
from random import Random
from sage.all import QQ, ZZ, GF, PolynomialRing, Polyhedron, matrix, prod, factorial
from sage.env import SAGE_VERSION
from verify_half_spin_degree import supports
from monomial_staircase import staircase_length

ROOT=Path(__file__).resolve().parents[1]

def degree(r):
    return 2**(r-1)*(r*r-3*r-2)+4*r

def combinatorial_checks():
    rows=[]
    for r in range(3,41):
        m=r-1
        strata=[(k,l,m-k-l) for k in range(2,m+1) for l in range(2,m-k+1)]
        facet=sum(factorial(m)//(factorial(k)*factorial(l))*(k-1)*(l-1)
                  for k,l,f in strata)
        closed=factorial(m)-(m-2)*2**(m-1)-2
        assert facet==closed
        total=2*factorial(r)-2**r
        assert total-2*r*facet==degree(r)
        if r%2==0: assert degree(r)%2==0
        rows.append(dict(rank=r,facet_correction=int(facet),degree=int(degree(r))))
    # Symbolic generating-function identity, independent of the finite sums.
    # ((x-1)e^x+1)^2/(1-x)=(1-x)e^(2x)-2e^x+1/(1-x).
    R=PolynomialRing(QQ,['x','e']); x,e=R.gens()
    assert ((x-1)*e+1)**2==(1-x)**2*e**2-2*(1-x)*e+1
    return rows

def polytope_checks():
    rows=[]
    for sigma in [1,-1]:
        W=supports(6,sigma); P=Polyhedron(vertices=W,base_ring=QQ)
        assert factorial(6)*P.volume()==1376
        assert Counter(len(f.vertices()) for f in P.facets())=={6:32,16:12}
        assert Counter(len(f.vertices()) for f in P.faces(4))=={5:192,8:60}
        for d in range(4):
            for f in P.faces(d):
                assert len(f.vertices())==d+1
                assert matrix(QQ,[list(v) for v in f.vertices()]).rank()==d+1
        assert prod(d for d in matrix(ZZ,W).elementary_divisors() if d)==2
        rows.append(dict(rank=6,chirality=sigma,volume=1376,
            coordinate_facets=12,simplex_facets=32,
            nonsimplex_codimension_two_faces=60,weight_sublattice_index=2))
    return rows

def local_checks():
    # Enumerate facet strata directly as sign/free patterns.
    patterns=list(product([-1,0,1],repeat=5))
    curves=[p for p in patterns if p.count(-1)==2 and p.count(1)==2]
    triples=[p for p in patterns if sorted([p.count(-1),p.count(1)])==[2,3]]
    assert len(curves)==30 and len(triples)==20
    assert all(sum(all(a==0 or a==b for a,b in zip(c,t)) for c in curves)==3
               for t in triples)
    # The square-free product gradient has k-1 solutions for nonzero constants.
    product_rows=[]
    for k in range(2,8):
        R=PolynomialRing(QQ,['a'+str(i) for i in range(k)]+['v'],order='degrevlex')
        a=R.gens()[:k]; v=R.gen(k); P=prod(a)
        eq=[P.derivative(x)-(i+2) for i,x in enumerate(a)]+[v*P-1]
        gb=R.ideal(eq).groebner_basis()
        length=staircase_length([tuple(g.lm().exponents()[0]) for g in gb])
        assert length==k-1
        product_rows.append(dict(number_of_variables=k,length=int(length)))
    # Independent exact free-coordinate critical equations on full cubes.
    cube_rows=[]
    rng=Random(687013)
    for f in range(1,5):
        R=PolynomialRing(GF(32003),['q'+str(i) for i in range(f)]+['v'],order='degrevlex')
        q=R.gens()[:f]; v=R.gen(f)
        coefficients=[rng.randrange(1,32003) for _ in range(2**f)]
        H=sum(c*prod(x**a for x,a in zip(q,e))
              for c,e in zip(coefficients,product([0,1],repeat=f)))
        eq=[x*H.derivative(x)-(i+3)*H for i,x in enumerate(q)]+[v*prod(q)*H-1]
        gb=R.ideal(eq).groebner_basis()
        length=staircase_length([tuple(g.lm().exponents()[0]) for g in gb])
        assert length==factorial(f)
        cube_rows.append(dict(free_coordinates=f,length=int(length),coefficients=coefficients))
    # Uniform first-jet freedom: on each fixed-sign block, parity weights
    # affinely span. Value and every transverse derivative can be prescribed.
    jets=[]
    for s in range(4,10):
        for parity in [0,1]:
            E=[e for e in product([0,1],repeat=s) if sum(e)%2==parity]
            assert matrix(QQ,[[1]+list(e) for e in E]).rank()==s+1
            jets.append(dict(fixed_coordinates=s,parity=parity,jet_rank=s+1))
    # Exact local product identities in Mobius coordinates.
    mobius=[]
    for k,l in [(2,2),(2,3),(3,2),(2,4),(3,3),(4,2)]:
        R=PolynomialRing(QQ,['u'+str(i) for i in range(k+l)]); K=R.fraction_field()
        u=list(map(K,R.gens())); z=[(x-1)/(x+1) for x in u[:k]]+[(1-x)/(1+x) for x in u[k:]]
        den=prod(1+x for x in u)
        for sigma in [1,-1]:
            F=prod(1+x for x in z)+sigma*prod(1-x for x in z)
            assert F*den/2**(k+l)==prod(u[:k])+sigma*prod(u[k:])
        mobius.append(dict(minus=k,plus=l))
    return dict(Spin12_curves_per_facet=30,Spin12_triple_points_per_facet=20,
        each_triple_has_three_curve_branches=True,
        Spin12_curve_contribution=360,Spin12_triple_contribution=480,
        product_gradient_checks=product_rows,free_cube_checks=cube_rows,
        transverse_jet_checks=jets,Mobius_product_identities=mobius)

def fiber(p,sigma,model):
    r=6; A=GF(p)
    if model=='weight_lattice':
        R=PolynomialRing(A,['z1','z2','z3','z4','z5','b','v'],order='degrevlex')
        vs=R.gens()[:r]; inv=R.gen(r); W=supports(r,sigma)
        mins=[min(0,min(w[i] for w in W)) for i in range(r)]
        den=prod(v**(-m) for v,m in zip(vs,mins))
        mon=[prod(v**(e-m) for v,e,m in zip(vs,w,mins)) for w in W]
        target=[2,3,5,7,11,17]
        eq=[sum(A(w[i])*term for w,term in zip(W,mon))-A(target[i])*den for i in range(r)]
        eq.append(inv*prod(vs)-1)
    else:
        R=PolynomialRing(A,['z1','z2','z3','z4','z5','z6','b','v'],order='degrevlex')
        z=R.gens()[:r]; b,v=R.gens()[r:]; D=prod(z)
        F=(prod(1+x for x in z)+sigma*prod(1-x for x in z))/2
        target=[1,2,3,5,7,11]
        eq=[2*x*F.derivative(x)-F-16*b*y for x,y in zip(z,target)]
        eq += [b*b-D,v*b*D-1]
        vs=list(z)+[b]
    gb=R.ideal(eq).groebner_basis(algorithm='libsingular:slimgb')
    lms=[tuple(g.lm().exponents()[0]) for g in gb]
    length=staircase_length(lms); assert length==536
    # Independent Sage normal basis confirms the native staircase count.
    basis=list(R.ideal(gb).normal_basis()); assert len(basis)==536
    # One linear multiplication operator certifies geometric reducedness
    # when its characteristic polynomial is squarefree.
    positions={tuple(b.exponents()[0]):i for i,b in enumerate(basis)}
    ell=vs[-1]+2*vs[0]+3*vs[1]
    entries={}
    for j,b in enumerate(basis):
        for e,c in (b*ell).reduce(gb).dict().items(): entries[positions[tuple(e)],j]=c
    M=matrix(A,length,length,entries); f=M.charpoly('q')
    gcd_degree=f.gcd(f.derivative()).degree()
    # A nonseparating linear form is reported honestly; length remains exact.
    return dict(prime=p,chirality=sigma,model=model,target=target,length=int(length),
        squarefree_multiplication_characteristic_polynomial=(gcd_degree==0),
        characteristic_polynomial_gcd_degree=int(gcd_degree),linear_form=str(ell),
        leading_monomials=[list(e) for e in lms],
        characteristic_polynomial_coefficients=[int(c) for c in f.list()])

def main():
    counts=combinatorial_checks(); polytopes=polytope_checks(); local=local_checks()
    fibers=[]
    for p,sigma in [(32003,1),(32009,-1)]:
        for model in ['weight_lattice','product_character']:
            fibers.append(fiber(p,sigma,model))
            print('Verified Spin12 fiber',p,sigma,model,flush=True)
    paths=['HALF_SPIN_RANK_DEGREES.md','tools/verify_half_spin_rank_degrees.py',
           'tools/verify_half_spin_degree.py','tools/monomial_staircase.py',
           'tools/explore_half_spin_rank6_product.py','experiments/half-spin-rank7-fiber.json']
    rank7=json.loads((ROOT/'experiments/half-spin-rank7-fiber.json').read_text())
    assert rank7['fiber_length']==degree(7)==1692
    for p,digest in rank7['source_sha256'].items():
        assert sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest()==digest
    report=dict(result='passed',date='2026-10-22',sage_version=SAGE_VERSION,
        all_rank_characteristic_zero_written_proof=True,rank_minimum=3,
        degree_formula='2^(r-1)*(r^2-3*r-2)+4*r',
        Spin12_degree=536,Spin12_total_intersection=1376,Spin12_boundary_correction=840,
        deeper_face_noncontribution_written_proof=True,
        actual_Spin_character_lattice=True,even_rank_representation_kernel_order=2,
        combinatorial_checks=counts,polytope_checks=polytopes,local_checks=local,
        independent_fibers=fibers,finite_fiber_checks=len(fibers),
        additional_completed_rank7_fiber=rank7,
        squarefree_fiber_certificates=sum(r['squarefree_multiplication_characteristic_polynomial'] for r in fibers),
        novelty_claim=False,externally_reviewed=False,runtime_API_changed=False,
        scope='All-rank degree is a written deformation deduction. Finite polytopes, local models and specialized original fibers corroborate it, not a formal proof or priority certificate.',
        source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'half-spin-rank-degrees-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in
        ('combinatorial_checks','polytope_checks','local_checks','independent_fibers','source_sha256')},indent=2))

if __name__=='__main__': main()
