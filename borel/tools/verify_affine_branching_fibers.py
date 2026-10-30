"""Symbolic and independent weight checks of the affine branching criterion."""
from hashlib import sha256
from pathlib import Path
from itertools import combinations
from math import comb
import json

from sage.all import (QQ, PolynomialRing, CyclotomicField, WeylCharacterRing,
                      matrix, vector, prod)
from sage.env import SAGE_VERSION

ROOT=Path(__file__).resolve().parents[1]


def coefficient_checks():
    arows=[]; crows=[]
    for N in range(4,13):
        for k in range(1,N):
            for m in range(1,min(k,N-k)+1):
                q=N-m
                R=PolynomialRing(QQ,['a'+str(i) for i in range(1,m)]+['b'+str(i) for i in range(1,q+1)])
                gs=R.gens(); aa=gs[:m-1]; bb=gs[m-1:]; F=R.fraction_field()
                A=[F(1)]+list(aa)+[1/F(bb[-1])]
                B=[F(1)]+list(bb)
                chi=sum(A[i]*B[k-i] for i in range(m+1) if 0<=k-i<=q)
                actual=[chi.derivative(x) for x in aa]
                expected=[F(bb[k-i-1]) for i in range(1,m)]
                assert actual==expected
                assert len(set(actual))==m-1
                assert all(1<=k-i<=q-1 for i in range(1,m))
                arows.append(dict(N=N,k=k,m=m,codimension=m-1,excess=0))
    for n in range(2,14):
        for p in range(1,n):
            q=n-p
            R=PolynomialRing(QQ,['a'+str(i) for i in range(1,p+1)]+['b'+str(i) for i in range(1,q+1)])
            aa=R.gens()[:p]; bb=R.gens()[p:]; T=PolynomialRing(R,'t'); t=T.gen()
            def recip(rank,coords):
                return 1+t**(2*rank)+sum(coords[i-1]*(t**i+t**(2*rank-i)) for i in range(1,rank))+coords[-1]*t**rank
            chi=((1-t*t)*recip(p,aa)*recip(q,bb))[n]
            ideal=R.ideal([chi.derivative(a) for a in aa])
            if p>q:
                assert ideal==R.ideal(1)
                crows.append(dict(n=n,p=p,q=q,empty_regular_fiber=True))
            else:
                def b(j): return R(0) if j<0 else (R(1) if j==0 else bb[j-1])
                expected=R.ideal([b(j)-b(j-2) for j in range(q-p+1,q+1)])
                assert ideal==expected
                J=matrix(QQ,[[f.derivative(x) for x in bb] for f in expected.gens()])
                assert J.rank()==p
                crows.append(dict(n=n,p=p,q=q,codimension=p,excess=0))
    return arows,crows


def algebraic_controls():
    C=PolynomialRing(QQ,['c1','c2']); c1,c2=C.gens(); F=C.fraction_field()
    R=PolynomialRing(F,['x','z','w']); x,z,w=R.gens()
    # The singular coefficient curve (z-1)(w-1)=0 has only target (0,0)
    # above its crossing. Generic c1,c2 avoid both smooth components too.
    f=(z-1)*(w-1)
    I=R.ideal(f,z*x*(w-1)-c1,w*x*(z-1)-c2)
    assert I==R.ideal(1)
    shifted=R.ideal(f,z*(x*(w-1)+1)-c1,w*(x*(z-1)+1)-c2)
    assert shifted==shifted.radical() and shifted.vector_space_dimension()==2
    assert shifted+R.ideal(z-1,w-1)==R.ideal(1)
    # A nonreduced coefficient zero cannot dominate the one-parameter target.
    J=R.ideal((z-1)**2,2*x*z*(z-1)-c1)
    assert J==R.ideal(1)
    # Dropping affinity permits persistent double critical structure.
    K=R.ideal(3*x*x,z-c1,w-c2)
    assert K!=K.radical() and K.radical()==R.ideal(x,z-c1,w-c2)
    # This is W=x^3+z+w-c1 log z-c2 log w.
    # An empty coefficient list does not imply nonempty generic fibers:
    # W=-c1 log z has critical equation c1=0.
    assert R.ideal(F(c1))==R.ideal(1)
    return dict(singular_coefficient_crossing_avoided=True,
                shifted_singular_curve_has_two_reduced_generic_critical_points=True,
                nonreduced_coefficient_component_avoided=True,
                nonlinear_potential_has_persistent_double_structure=True,
                generic_nonemptiness_not_automatic=True)


def type_c_weights():
    rows=[]
    for n in range(2,8):
        W=WeylCharacterRing('C'+str(n),style='coroots')
        weights=W(*([0]*(n-1)+[1])).weight_multiplicities()
        for p in range(1,n//2+1):
            q=n-p; A=CyclotomicField(2*q+2); zeta=A.gen()
            ts=[A(2+j) for j in range(p)]+[zeta**j for j in range(1,q+1)]
            moment=vector(A,n); Hess=matrix(A,n); Gram=matrix(QQ,n)
            for weight,mult in weights.items():
                v=vector(QQ,list(weight.to_vector())[:n])
                term=prod(ts[j]**int(v[j]) for j in range(n))
                moment+=mult*term*v; Hess+=mult*term*v.column()*v.row()
                Gram+=mult*v.column()*v.row()
            image=Gram.inverse()*moment; torusJ=Gram.inverse()*Hess
            assert all(image[j]==0 for j in range(p))
            assert all(image[j]!=0 for j in range(p,n))
            assert all(image[i]!=s*image[j] for i,j in combinations(range(p,n),2) for s in [-1,1])
            # Source regularity and all root-direction differential ranks.
            roots=[]
            for i in range(n):
                roots.append((ts[i]**2,2*image[i]))
            for i,j in combinations(range(n),2):
                roots.extend([(ts[i]/ts[j],image[i]-image[j]),(ts[i]*ts[j],image[i]+image[j])])
            assert all(value!=1 for value,_ in roots)
            fullrank=torusJ.rank()+2*sum(hval!=0 for _,hval in roots)
            assert torusJ.rank()==n,(n,p,'torus Hessian',torusJ.rank())
            assert n*(2*n+1)-fullrank==2*p*p
            rows.append(dict(n=n,p=p,q=q,representation_dimension=int(sum(weights.values())),
                             torus_rank=int(torusJ.rank()),full_differential_rank=int(fullrank),
                             fiber_tangent_dimension=2*p*p,target_stratum_verified=True))
    return rows


def type_a_weights():
    rows=[]
    for N,k,m in [(5,2,2),(6,2,2),(6,3,3),(7,3,3),(8,4,3)]:
        q=N-m
        # Root-of-unity external spectrum sets all intermediate elementary
        # coefficients to zero, in particular the required contiguous block.
        A=CyclotomicField(q); zeta=A.gen()
        bs=[A(2)*zeta**j for j in range(q)]
        first=[A(3+j) for j in range(m-1)]
        first.append(1/prod(first+bs)); ts=first+bs
        assert prod(ts)==1 and len(set(ts))==N
        H=matrix(QQ,N,N-1,lambda i,j:1 if i==j else (-1 if i==N-1 else 0))
        moment=vector(A,N-1); Hess=matrix(A,N-1); Gram=matrix(QQ,N-1)
        for subset in combinations(range(N),k):
            weight=vector(QQ,[int(i in subset) for i in range(N)])
            v=H.transpose()*weight; term=prod(ts[i] for i in subset)
            moment+=term*v; Hess+=term*v.column()*v.row(); Gram+=v.column()*v.row()
        image=H*(Gram.inverse()*moment)
        assert all(image[j]==image[0] for j in range(m))
        assert all(image[j]!=image[0] for j in range(m,N))
        assert len(set(image[m:]))==q
        torusrank=Hess.rank()
        fullrank=torusrank+2*sum(image[i]!=image[j] for i,j in combinations(range(N),2))
        assert torusrank==N-1,(N,k,m,torusrank)
        assert N*N-1-fullrank==m*(m-1)
        rows.append(dict(N=N,k=k,m=m,representation_dimension=comb(N,k),
                         torus_rank=int(torusrank),full_differential_rank=int(fullrank),
                         fiber_tangent_dimension=m*(m-1),target_stratum_verified=True))
    return rows


def main():
    aa,cc=coefficient_checks(); controls=algebraic_controls()
    print('Coefficient identities and singular/nonlinear controls passed',flush=True)
    cweights=type_c_weights(); aweights=type_a_weights()
    report=dict(result='passed',date='2026-10-23',sage_version=SAGE_VERSION,
                type_A_coefficient_checks=aa,type_C_coefficient_checks=cc,
                algebraic_controls=controls,type_A_weight_differentials=aweights,
                type_C_Weyl_weight_differentials=cweights,
                general_proof_is_written=True,all_fundamental_branchings_affine=False,
                nonempty_generic_fibers_proved_for_all_admissible_cases=False,
                exceptional_target_values_classified=False)
    sources=['AFFINE_BRANCHING_FIBERS.md','tools/verify_affine_branching_fibers.py',
             'GENERAL_CAYLEY_FIBER_CRITICAL.md','SYMPLECTIC_TARGET_STRATA.md']
    report['source_sha256']={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in sources}
    (ROOT/'affine-branching-fibers-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(result='passed',type_A_coefficients=len(aa),type_C_coefficients=len(cc),
                         type_A_weight_cases=len(aweights),type_C_weight_cases=len(cweights)),indent=2))


if __name__=='__main__': main()
