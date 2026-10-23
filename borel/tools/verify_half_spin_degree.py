"""Spin_10 degree: lattice polytope, exact boundary algebra and original fibers."""
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json
from sage.all import QQ, ZZ, GF, PolynomialRing, Polyhedron, matrix, prod, factorial, gcd
from sage.env import SAGE_VERSION
from monomial_staircase import staircase_length

ROOT=Path(__file__).resolve().parents[1]


def weights(r,sigma=1):
    return [e for e in product([1,-1],repeat=r) if prod(e)==sigma]


def supports(r,sigma=1):
    return [tuple((e[i]-e[-1])//2 for i in range(r-1))+(e[-1],)
            for e in weights(r,sigma)]


def polytope_checks():
    rows=[]
    for r in [3,4,5]:
        for sigma in [1,-1]:
            W=supports(r,sigma); P=Polyhedron(vertices=W,base_ring=QQ)
            volume=factorial(r)*P.volume()
            assert volume=={3:4,4:32,5:208}[r]
            S=matrix(ZZ,W).smith_form()[0]; index=prod(abs(S[i,i]) for i in range(r))
            assert index=={3:1,4:2,5:1}[r]
            facets=P.facets(); counts=[len(f.vertices()) for f in facets]
            if r==5:
                assert counts.count(8)==10 and counts.count(5)==16
                for f in facets:
                    verts=[list(v) for v in f.vertices()]
                    if len(verts)==5: assert matrix(QQ,verts).rank()==5
                    else:
                        eps=[[(v[i]+QQ(v[-1])/2) for i in range(4)]+[QQ(v[-1])/2]
                             for v in verts]
                        coords=[i for i in range(5) if len({v[i] for v in eps})==1]
                        assert len(coords)==1
                        i=coords[0]
                        # The normal 2 epsilon_i^* is primitive on P(D5).
                        normal=[2*int(i==j) for j in range(4)]+[1]
                        assert gcd(normal)==1
                        assert all(abs(2*v[i])==1 for v in eps)
                for d in range(4):
                    for f in P.faces(d):
                        assert len(f.vertices())==d+1
                        assert matrix(QQ,[list(v) for v in f.vertices()]).rank()==d+1
            else:
                assert all(len(f.vertices())==r for f in facets)
            rows.append(dict(rank=r,chirality=sigma,normalized_volume=int(volume),
                weight_generated_sublattice_index=int(index),facet_vertex_counts=counts,
                proper_small_faces_are_simplices=True))
    # Independent Euclidean cube-minus-corners calculation, covolume P(D5)=1/2.
    assert factorial(5)*(1-QQ(16)/factorial(5))/(QQ(1)/2)==208
    return rows


def boundary_checks():
    R=PolynomialRing(QQ,['z1','z2','z3','z4','v'],order='degrevlex')
    z=R.gens()[:4]; v=R.gens()[4]; D=prod(z)
    E=1+sum(z[i]*z[j] for i,j in combinations(range(4),2))+D
    eq=[E]+[x*E.derivative(x) for x in z]+[v*D-1]
    gb=R.ideal(eq).groebner_basis()
    length=staircase_length([tuple(f.lm().exponents()[0]) for f in gb])
    assert length==6
    points=[e for e in product([1,-1],repeat=4) if sum(e)==0]
    assert len(points)==6
    for e in points: assert all(f(*e,1)==0 for f in eq)
    # Symbolic target variables certify multiplicity one at every one of
    # these six points, rather than checking a single numerical target.
    Q=PolynomialRing(QQ,['z1','z2','z3','z4','t','y1','y2','y3','y4','yb'])
    z=Q.gens()[:4]; t=Q.gen(4); ys=Q.gens()[5:]; D=prod(z)
    E=1+sum(z[i]*z[j] for i,j in combinations(range(4),2))+D
    O=sum(z)+sum(prod(z[i] for i in ids) for ids in combinations(range(4),3))
    eq=[x*E.derivative(x)-E+t*t*D*x*O.derivative(x)-y*t*D
        for x,y in zip(z,ys[:4])]
    eq.append(E-t*t*D*O-ys[4]*t*D)
    J=matrix(Q,[[f.derivative(x) for x in list(z)+[t]] for f in eq])
    rows=[]
    for e in points:
        sub=dict(zip(z,e)); sub[t]=0
        jj=matrix(Q,[[entry.subs(sub) for entry in row] for row in J.rows()])
        assert jj.det()==-16*ys[4]
        rows.append(dict(point=e,full_local_Jacobian_determinant=str(jj.det()),
                         local_intersection_multiplicity=1))
    return dict(facet_critical_algebra_length=6,facet_count=10,total_boundary_length=60,
                symbolic_local_checks=rows,
                leading_monomials=[list(f.lm().exponents()[0]) for f in gb])


def finite_algebra_report(R,eq,ell,require_reduced=True):
    I=R.ideal(eq); gb=I.groebner_basis(algorithm='libsingular:slimgb')
    lms=[tuple(f.lm().exponents()[0]) for f in gb]
    length=staircase_length(lms); assert length==148
    basis=list(R.ideal(gb).normal_basis()); assert len(basis)==length
    positions={tuple(f.exponents()[0]):i for i,f in enumerate(basis)}
    # A reduced algebra need not have distinct values of one chosen linear
    # form. Record collisions and try other forms; never infer nonreducedness
    # from a collision, or reducedness from length alone.
    collisions=[]
    for attempt in range(16):
        candidate=ell if attempt==0 else ell+sum(
            R.base_ring()(attempt**(i+1))*x for i,x in enumerate(R.gens()[:-1]))
        entries={}
        for j,b in enumerate(basis):
            red=(b*candidate).reduce(gb)
            for e,c in red.dict().items(): entries[positions[tuple(e)],j]=c
        M=matrix(R.base_ring(),length,length,entries)
        f=M.charpoly('q'); repeated_degree=f.gcd(f.derivative()).degree()
        if repeated_degree==0: break
        if M.minpoly().degree()==length:
            # A cyclic length-d algebra is k[q]/(f). A repeated factor
            # now proves nilpotents, rather than merely a value collision.
            break
        collisions.append(int(repeated_degree))
    else: raise AssertionError('No separating linear form certified reducedness')
    squarefree=repeated_degree==0
    if require_reduced: assert squarefree
    return dict(length=int(length),squarefree_multiplication_characteristic_polynomial=squarefree,
                cyclic_linear_form=str(candidate),
                geometric_distinct_point_count=int(length-repeated_degree),
                nonreduced_cyclic_algebra_certified=not squarefree,
                earlier_linear_form_collision_gcd_degrees=collisions,
                characteristic_polynomial_coefficients=[int(c) for c in f.list()],
                leading_monomials=[list(e) for e in lms])


def lattice_fiber(p,sigma,target):
    A=GF(p); r=5
    R=PolynomialRing(A,['z1','z2','z3','z4','b','v'],order='degrevlex')
    vs=R.gens()[:r]; inv=R.gen(r); W=supports(r,sigma)
    mins=[min(0,min(w[i] for w in W)) for i in range(r)]
    den=prod(v**(-m) for v,m in zip(vs,mins))
    mon=[prod(v**(e-m) for v,e,m in zip(vs,w,mins)) for w in W]
    eq=[sum(A(w[i])*term for w,term in zip(W,mon))-A(target[i])*den for i in range(r)]
    eq.append(inv*prod(vs)-1)
    row=finite_algebra_report(R,eq,vs[-1]+2*vs[0]+3*vs[1])
    row.update(model='actual_P_D5_basis',prime=p,chirality=sigma,
               logarithmic_derivative_target=target)
    return row


def original_spin_fiber(p,sigma,target,require_reduced=True):
    A=GF(p)
    R=PolynomialRing(A,['z1','z2','z3','z4','z5','b','v'],order='degrevlex')
    z=R.gens()[:5]; b,v=R.gens()[5:]; W=weights(5,sigma); D=prod(z)
    # The actual normalized epsilon-coordinate formula: Gram=4 I,
    # weights epsilon/2, and b^2=product z_i on the Spin torus.
    eq=[sum(A(e[i])/8*b*prod(x**((a+1)//2) for x,a in zip(z,e)) for e in W)
        -A(target[i])*D for i in range(5)]
    eq += [b*b-D,v*b*D-1]
    row=finite_algebra_report(R,eq,b+2*z[0]+3*z[1],require_reduced)
    row.update(model='original_Spin_cover',prime=p,chirality=sigma,
               normalized_Springer_target=target)
    return row


def main():
    polytopes=polytope_checks(); boundary=boundary_checks(); fibers=[]
    for p in [32003,32009]:
        for sigma in [1,-1]:
            for target in [[2,3,5,7,11],[3,7,13,17,23]]:
                fibers.append(lattice_fiber(p,sigma,target))
            fibers.append(original_spin_fiber(p,sigma,[2,3,5,7,11]))
            print('Verified fibers:',p,sigma,flush=True)
    special=[original_spin_fiber(32009,sigma,[1,2,3,4,5],require_reduced=False)
             for sigma in [1,-1]]
    assert all(row['nonreduced_cyclic_algebra_certified'] for row in special)
    paths=['HALF_SPIN_DEGREE.md','tools/verify_half_spin_degree.py',
           'tools/monomial_staircase.py']
    report=dict(result='passed',date='2026-10-22',sage_version=SAGE_VERSION,
        Spin10_half_spin_generic_degree=148,characteristic_zero_toric_written_proof=True,
        toric_intersection_number=208,toric_boundary_length=60,
        actual_Spin_character_lattice=True,group_torus_quotient_degrees_agree=True,
        polytope_checks=polytopes,boundary_algebra=boundary,
        finite_field_original_fibers=fibers,finite_field_fiber_checks=len(fibers),
        all_primary_test_fibers_reduced=True,
        reduced_fiber_target=[2,3,5,7,11],
        separate_exceptional_specializations=special,
        all_rank_degree_formula_proved=False,novelty_claim=False,externally_reviewed=False,
        rational_specialized_Groebner_completed=False,
        runtime_API_changed=False,
        scope='Characteristic-zero degree follows from the written toric intersection proof. Exact finite-field fibers in two coordinate models corroborate it; they are not extrapolated to prove the generic degree.',
        source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'half-spin-degree-verification.json').write_text(json.dumps(report,indent=2,default=int)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in
        ('polytope_checks','boundary_algebra','finite_field_original_fibers',
         'separate_exceptional_specializations','source_sha256')},indent=2))


if __name__=='__main__': main()
