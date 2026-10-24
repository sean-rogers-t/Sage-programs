"""Classical defining/adjoint degrees: root geometry and independent fibers."""
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import argparse
import json
from sage.all import (QQ, ZZ, GF, matrix, vector, prod, factorial, binomial,
                      Polyhedron, PolynomialRing, LaurentPolynomialRing, WeylCharacterRing, gcd)
from sage.env import SAGE_VERSION
from monomial_staircase import staircase_length

ROOT=Path(__file__).resolve().parents[1]


def exceptional_ratios(n):
    return [(4*k-2*n+1,k,int(2*n*binomial(n-1,k)))
            for k in range(n) if 4*k-2*n+1>=0]


def degree(kind,n,a,b):
    assert kind in ['B','C'] and n>=1 and a>=0 and b>=0 and a+b>0
    if b==0: return 2**n
    if kind=='C': return 4**n
    correction=sum(v for ratio,k,v in exceptional_ratios(n) if a==ratio*b)
    return 2**n*(2**n-n)-correction


def weights(kind,n,a,b):
    E=list(matrix.identity(ZZ,n).rows())
    short=[tuple(s*v) for v in E for s in [-1,1]]
    pairs=[tuple(s*E[i]+t*E[j]) for i,j in combinations(range(n),2)
           for s,t in product([-1,1],repeat=2)]
    doubled=[tuple(2*s*v) for v in E for s in [-1,1]]
    if kind=='B': out={w:a+b for w in short}; out.update({w:b for w in pairs})
    else: out={w:a for w in short}; out.update({w:b for w in pairs+doubled})
    return {w:m for w,m in out.items() if m}


def geometry():
    rows=[]
    for kind in ['B','C']:
        for n in range(2,7):
            chars=WeylCharacterRing([kind,n],style='coroots')
            defining=[1]+[0]*(n-1)
            # In B2 the vector module is omega_1, adjoint is 2*omega_2.
            adjoint=([0,2] if n==2 else [0,1]+[0]*(n-2)) if kind=='B' else [2]+[0]*(n-1)
            for high,a,b,zero,dim in [(defining,1,0,int(kind=='B'),2*n+int(kind=='B')),
                                      (adjoint,0,1,n,n*(2*n+1))]:
                ch=chars(*high)
                mult={tuple(w.to_vector()):int(m) for w,m in ch.weight_multiplicities().items()}
                assert ch.degree()==dim
                assert mult.pop((0,)*n,0)==zero
                assert mult==weights(kind,n,a,b),(kind,n,high)
            W=weights(kind,n,1,1)
            moment=sum((m*vector(QQ,w).column()*vector(QQ,w).row() for w,m in W.items()),matrix(QQ,n,n))
            scale=2+(4*n-2 if kind=='B' else 4*(n+1))
            assert moment==scale*matrix.identity(QQ,n)
            P=Polyhedron(vertices=list(W))
            expected=2**n*(2**n-n) if kind=='B' else 4**n
            volume=factorial(n)*P.volume(); assert volume==expected,(kind,n,volume)
            counts=Counter(); nonvertex_weights=set(W)-{tuple(v) for v in P.vertices()}
            for face in P.facets():
                h=face.ambient_Hrepresentation()[0]
                assert gcd(list(h)[1:])==1
                normal=list(h)[1:]
                face_support=[w for w in W if h[0]+sum(x*y for x,y in zip(normal,w))==0]
                if kind=='B' and sum(v!=0 for v in normal)==1:
                    assert h[0]==1 and len(face.vertices())==2*(n-1)
                    assert len(face_support)==2*(n-1)+1; counts['coordinate_cross']+=1
                elif kind=='B':
                    assert h[0]==2 and all(abs(v)==1 for v in normal)
                    assert len(face_support)==binomial(n,2); counts['hypersimplex']+=1
                else:
                    assert h[0]==2 and all(abs(v)==1 for v in normal)
                    assert len(face.vertices())==n and len(face_support)==n*(n+1)//2
                    counts['doubled_simplex']+=1
            if kind=='B':
                assert counts['coordinate_cross']==2*n
                assert counts['hypersimplex']==(2**n if n>=3 else 0)
                # The short-root centers are on no lower-dimensional face.
                for d in range(n-1):
                    for face in P.faces(d):
                        assert not any(face.as_polyhedron().contains(w) for w in nonvertex_weights)
            else: assert counts=={'doubled_simplex':2**n}
            # For pure adjoint C the weight-generated lattice has index 2.
            snf=[abs(int(d)) for d in matrix(ZZ,list(weights(kind,n,0,1))).elementary_divisors() if d]
            assert prod(snf)==(1 if kind=='B' else 2)
            rows.append(dict(kind=kind,rank=n,normalized_volume=int(volume),facet_types=dict(counts),
                             independent_highest_weight_characters_agree=True,trace_moment=scale,
                             adjoint_support_lattice_index=int(prod(snf))))
    # Independent sign enumeration and exact local tangent Hessians, not extrapolated point counts.
    sign_checks=[]
    LR=LaurentPolynomialRing(QQ,'q'); q=LR.gen()
    for n in range(2,11):
        counts=Counter(-2*sum(s)-1 for s in product([-1,1],repeat=n-1))
        positive={int(r):int(c) for r,c in counts.items() if r>=0}
        assert positive=={r:int(binomial(n-1,k)) for r,k,c in exceptional_ratios(n)}
        for r,k,correction in exceptional_ratios(n):
            signs=[-1]*k+[1]*(n-1-k)
            assert r+1+2*sum(signs)==0
            assert prod(2*s for s in signs)!=0
            # An actual torus curve maps to zero, not just a compactification base point.
            z=[q]+[LR(s) for s in signs]
            W=weights('B',n,r,1)
            assert all(sum(m*w[i]*prod(t**e for t,e in zip(z,w)) for w,m in W.items())==0 for i in range(n))
            x=[t+1/t for t in z]; differences=[t-1/t for t in z]
            assert all(differences[i]*(r+1+sum(x)-x[i])==0 for i in range(n))
        sign_checks.append(dict(rank=n,points_per_facet=positive,
                                corrections={r:c for r,k,c in exceptional_ratios(n)},
                                exceptional_zero_curves_in_both_constructions=True))
    # Quadratic face-gradient matrices: e2 has J-I; h2 has J+I.
    # Their invertibility excludes zeros even before passing to face quotients.
    matrix_checks=[]
    for m in range(1,13):
        J=matrix(QQ,m,m,[1]*(m*m)); Id=matrix.identity(QQ,m)
        assert (J+Id).det()==m+1
        if m>=2: assert (J-Id).det()==(-1)**(m-1)*(m-1)
        matrix_checks.append(dict(variables=m,h2_determinant=int((J+Id).det()),
                                  e2_determinant=int((J-Id).det()) if m>=2 else None))
    return dict(polytope_and_character_checks=rows,sign_checks=sign_checks,
                face_gradient_matrix_checks=matrix_checks,
                arbitrary_rank_proof_in_note=True)


def fiber(kind,n,a,b,prime,model,spin_cover=False,adjoint_form=False):
    A=GF(prime); target=[1,3,7,13,23][:n]
    scale=2*a+(4*n-2 if kind=='B' else 4*(n+1))*b
    assert not (spin_cover and adjoint_form)
    if model=='weight_lattice':
        W=weights(kind,n,a,b)
        if adjoint_form:
            assert kind=='C' and a==0
            # Q(C_n) basis e_i-e_n (i<n), 2e_n.
            W={tuple(list(w[:-1])+[ZZ(sum(w)/2)]):m for w,m in W.items()}
            rhs=[scale*y for y in target[:-1]]+[A(scale)*sum(target)/2]
        else: rhs=[scale*y for y in target]
        names=['z'+str(i+1) for i in range(n)]+(['u'] if spin_cover else [])+['v']
        R=PolynomialRing(A,names,order='degrevlex'); z=R.gens()[:n]; inv=R.gens()[-1]
        mins=[min(0,min(w[i] for w in W)) for i in range(n)]
        den=prod(t**(-e) for t,e in zip(z,mins))
        mon={w:prod(t**(e-mi) for t,e,mi in zip(z,w,mins)) for w in W}
        equations=[sum(A(m*w[i])*mon[w] for w,m in W.items())-rhs[i]*den for i in range(n)]
        equations.append(inv*prod(z)-1)
        ell=sum((i+1)*z[i] for i in range(n))
        if spin_cover:
            u=R.gen(n); equations.append(u*u-prod(z)); ell+=u
    else:
        assert not spin_cover and not adjoint_form
        R=PolynomialRing(A,['x'+str(i+1) for i in range(n)]+['w'+str(i+1) for i in range(n)],order='degrevlex')
        x=R.gens()[:n]; w=R.gens()[n:]; total=sum(x)
        equations=[w[i]*w[i]-x[i]*x[i]+4 for i in range(n)]
        for i in range(n):
            factor=a+b+b*(total-x[i]) if kind=='B' else a+b*(total+x[i])
            equations.append(w[i]*factor-scale*target[i])
        # z_i=(x_i+w_i)/2, automatically a unit since x_i^2-w_i^2=4.
        ell=sum((i+1)*(x[i]+w[i])/2 for i in range(n))
    I=R.ideal(equations); gb=I.groebner_basis(algorithm='libsingular:slimgb')
    leading=[tuple(g.lm().exponents()[0]) for g in gb]
    length=staircase_length(leading); normal=list(I.normal_basis())
    expected=degree(kind,n,a,b)*(2 if spin_cover else 1)//(2 if adjoint_form else 1)
    assert length==len(normal)==expected,(kind,n,a,b,prime,model,length,expected)
    positions={tuple(v.exponents()[0]):i for i,v in enumerate(normal)}; entries={}
    for j,v in enumerate(normal):
        for e,c in (ell*v).reduce(gb).dict().items(): entries[positions[tuple(e)],j]=c
    p=matrix(A,length,length,entries).charpoly('q')
    assert p.gcd(p.derivative()).degree()==0,(kind,n,a,b,prime,model,'not separated/reduced')
    return dict(kind=kind,rank=n,a=a,b=b,prime=prime,model=model,
                normalized_epsilon_target=target,spin_cover=spin_cover,adjoint_form=adjoint_form,
                saturated_fiber_length=length,reduced_by_squarefree_multiplication=True,
                characteristic_polynomial_coefficients=[int(v) for v in p.list()])


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--geometry-only',action='store_true')
    args=parser.parse_args(); geo=geometry()
    if args.geometry_only: print(json.dumps(geo,indent=2)); return
    fibers=[]
    cases=[('B',1,0,1),('B',2,0,1),('B',2,1,1),('B',3,0,1),('B',3,3,1),
           ('B',4,0,1),('B',4,1,1),('B',4,5,1),('B',4,2,1),('B',3,1,0),
           ('C',1,0,1),('C',2,0,1),('C',2,1,1),('C',3,0,1),('C',3,3,1),
           ('C',4,1,1),('C',3,1,0)]
    for prime in [32003,32009]:
        for kind,n,a,b in cases:
            pair=[]
            for model in ['weight_lattice','trace_coordinates']:
                row=fiber(kind,n,a,b,prime,model); pair.append(row)
                print('Verified',kind,n,a,b,prime,model,'length',row['saturated_fiber_length'],flush=True)
            assert pair[0]['characteristic_polynomial_coefficients']==pair[1]['characteristic_polynomial_coefficients']
            fibers+=pair
    controls=[fiber('B',3,3,1,32003,'weight_lattice',spin_cover=True),
              fiber('C',2,0,1,32003,'weight_lattice',adjoint_form=True)]
    report=dict(result='passed',date='2026-10-22',sage_version=SAGE_VERSION,geometry=geo,
                independent_fibers=fibers,isogeny_controls=controls,
                characteristic_zero_proof='BC_MIXTURE_DEGREES.md',
                formula_B='2^n for b=0; 2^n*(2^n-n)-2*n*binomial(n-1,k) if a=(4*k-2*n+1)*b with 0<=k<=n-1; otherwise 2^n*(2^n-n)',
                formula_C='2^n for b=0; 4^n for b>0',
                spin_degree_multiplier=2,pure_C_adjoint_form_degree_divisor=2,
                same_target_multiplication_polynomials_agree=True,
                exceptional_B_torus_and_quotient_nonfinite_by_zero_curves=True,
                finite_fibers_used_as_proof_of_generic_degree=False,
                novelty_claim=False,externally_reviewed=False,runtime_API_changed=False,
                independence='Root-weight differentiation vs trace/difference coordinates; separate geometry; shared Sage/Singular arithmetic.')
    paths=['BC_MIXTURE_DEGREES.md','tools/verify_bc_mixture_degrees.py','tools/monomial_staircase.py']
    report['source_sha256']={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths}
    (ROOT/'bc-mixture-degrees-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['independent_fibers','source_sha256','geometry']},indent=2))


if __name__=='__main__': main()
