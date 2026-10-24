"""Exact independent checks for polynomial-determinant Springer degrees.

The proof is uniform; these finite computations corroborate it. Characters
are built from exterior weights and tensor convolution. Small fibers are
also computed from their original Laurent moment equations by Groebner
bases, with only torus-coordinate saturation (no Cayley substitutions).
"""
from collections import defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json

from sage.all import QQ, ZZ, PolynomialRing, matrix, prod, vector
from sage.env import SAGE_VERSION
from monomial_staircase import staircase_length
from verify_mixed_regularity import sp_basis, exterior_derivative

ROOT=Path(__file__).resolve().parents[1]
Z=PolynomialRing(ZZ,'d')
d=Z.gen()


def examples():
    return [('power_'+str(q),d**q) for q in range(1,6)]+[
        ('sum_1_2',d+d**2),('sum_1_3',2*d+3*d**3+5),
        ('removed_trivial_square',(d-1)**2),
        ('removed_trivial_cube',(d-1)**3),
        ('rank_one_standard_square',(d-2)**2)]


def convolve(a,b):
    out=defaultdict(int)
    for u,c in a.items():
        for v,e in b.items():
            out[tuple(x+y for x,y in zip(u,v))]+=c*e
    return dict(out)


@lru_cache(None)
def exterior_weights(n):
    zero=(0,)*n
    weights={zero:1}
    # Each defining weight can be chosen or omitted in the exterior algebra.
    for j in range(n):
        for sign in [1,-1]:
            weight=tuple(sign if i==j else 0 for i in range(n))
            weights=convolve(weights,{zero:1,weight:1})
    return weights


def character_weights(n,F):
    zero=(0,)*n
    power={zero:1}
    out=defaultdict(int)
    for k,a in enumerate(F.list()):
        for weight,c in power.items(): out[weight]+=int(a)*c
        if k<F.degree(): power=convolve(power,exterior_weights(n))
    out={w:c for w,c in out.items() if c}
    assert all(c>0 for c in out.values())
    return out


def moment_gram(n,weights):
    return matrix(QQ,n,n,lambda i,j:sum(c*w[i]*w[j] for w,c in weights.items()))


def moment_value(weights,coordinates):
    n=len(coordinates)
    values=[coordinates[0].parent().zero() for _ in range(n)]
    # Cache positive and negative powers so quotient-algebra checks stay small.
    powers=[{a:z**a for a in {w[i] for w in weights}} for i,z in enumerate(coordinates)]
    for w,c in weights.items():
        value=c*prod(powers[i][w[i]] for i in range(n))
        for i in range(n): values[i]+=w[i]*value
    return values


def fiber_polynomial(n,F,ys):
    R=PolynomialRing(QQ,'s')
    s=R.gen()
    q=int(F.degree()); K=4**n
    delta=prod(1-s*s*y*y for y in ys)
    H=sum(k*F[k]*K**(k-1)*delta**(q-k) for k in range(1,q+1))
    P=F.derivative()(K)*delta**q-2*s*H
    return P,delta,H


def scalar_checks():
    rows=[]
    for n in range(1,9):
        for label,F in examples():
            ys=[QQ(j+2) for j in range(n)]
            P,delta,H=fiber_polynomial(n,F,ys)
            degree=2*n*int(F.degree())
            assert P.degree()==degree
            assert P.gcd(P.derivative()).degree()==0
            assert P.gcd(delta).degree()==0
            assert P.gcd(P.parent().gen()).degree()==0
            assert H.gcd(delta).degree()==0
            assert P[0]==F.derivative()(4**n)>0
            rows.append(dict(rank=n,character=label,polynomial=str(F),
                degree=degree,reduced_fiber=True,no_denominator_roots=True))
    return rows


def trace_checks():
    rows=[]
    selected=[(label,F) for label,F in examples() if F.degree()<=3]
    for n in range(1,4):
        for label,F in selected:
            weights=character_weights(n,F)
            G=moment_gram(n,weights)
            K=4**n
            assert sum(weights.values())==F(K)
            assert G==matrix.identity(QQ,n)*(K/2)*F.derivative()(K)
            zs=[QQ(j+2) for j in range(n)]
            moment=vector(QQ,moment_value(weights,zs))
            D=prod(2+z+1/z for z in zs)
            scalar=2*D*F.derivative()(D)/(K*F.derivative()(K))
            assert G.inverse()*moment==vector(QQ,[scalar*(z-1)/(z+1) for z in zs])
            # Verify every reconstructed fiber root simultaneously in a finite algebra.
            ys=[QQ(j+2) for j in range(n)]
            P,delta,H=fiber_polynomial(n,F,ys)
            A=P.parent().quotient(P,'sbar')
            s=A.gen()
            zbar=[(1+s*y)/(1-s*y) for y in ys]
            values=moment_value(weights,zbar)
            image=G.change_ring(A).inverse()*vector(A,values)
            assert image==vector(A,ys)
            rows.append(dict(rank=n,character=label,actual_dimension=int(F(K)),
                distinct_weights=len(weights),fiber_algebra_dimension=int(P.degree()),
                independent_weight_moments=True))
        print('Independent weight/fiber checks passed at rank',n,flush=True)
    return rows


def original_fiber_checks():
    rows=[]
    selected=[('power_1',d),('power_2',d**2),('sum_1_2',d+d**2),
              ('removed_trivial_square',(d-1)**2)]
    cases=[(n,label,F) for n in [1,2] for label,F in selected]
    cases += [(1,'power_3',d**3),(3,'power_1',d)]
    for n,label,F in cases:
        q=int(F.degree())
        weights=character_weights(n,F)
        G=moment_gram(n,weights)
        assert G.is_diagonal()
        R=PolynomialRing(QQ,['z'+str(i+1) for i in range(n)]+['u'],order='degrevlex')
        zs=R.gens()[:n]; u=R.gen(n)
        monomial=prod(z**q for z in zs)
        equations=[]
        for i in range(n):
            p=sum(c*w[i]*prod(z**(a+q) for z,a in zip(zs,w))
                  for w,c in weights.items())
            equations.append(p-G[i,i]*(i+2)*monomial)
        equations.append(u*prod(zs)-1)
        ideal=R.ideal(equations)
        gb=ideal.groebner_basis()
        dimension=staircase_length([tuple(p.lm().exponents()[0]) for p in gb])
        assert dimension==2*n*q
        rows.append(dict(rank=n,character=label,torus_fiber_length=dimension,
            saturation='product(z_i), no Cayley denominators',coefficient_field='QQ',
            Groebner_basis_size=len(gb)))
        print('Original Laurent fiber passed:',n,label,dimension,flush=True)
    return rows


def singular_target_checks():
    rows=[]
    # Invertibility is needed for constant fiber length, not for surjectivity.
    for n in range(1,7):
        for m in range(1,5):
            for r in range(n+1):
                ys=[QQ(j+2) if j<r else QQ(0) for j in range(n)]
                P,delta,H=fiber_polynomial(n,d**m,ys)
                assert P.degree()==(2*r*m if r else 1)
                assert P.gcd(delta).degree()==0
                rows.append(dict(rank=n,tensor_power=m,nonzero_eigenvalue_pairs=r,
                                 degree_on_Cayley_open=int(P.degree())))
    return rows


def full_matrix_fiber_checks():
    n=2; K=4**n
    basis=sp_basis(n)
    derivatives=[[exterior_derivative(H,k) for k in range(5)] for H in basis]
    G=matrix(QQ,[[(H*J).trace() for J in basis] for H in basis])
    actual_G=matrix(QQ,[[sum((a*b).trace() for a,b in zip(ds,es))
                        for es in derivatives] for ds in derivatives])
    assert actual_G==(K/4)*G
    A0=matrix(QQ,[[2,1],[0,2]])
    jordan=matrix(QQ,4,4)
    jordan[:2,:2]=A0
    jordan[2:,2:]=-A0.transpose()
    nilpotent=matrix(QQ,4,4)
    nilpotent[0,1]=1; nilpotent[3,2]=-1; nilpotent[1,3]=1
    assert nilpotent**4==0 and nilpotent**3!=0
    singular=matrix.diagonal(QQ,[0,2,0,-2])
    targets=[('invertible_nonsemisimple',jordan),('regular_nilpotent',nilpotent),
             ('singular_semisimple',singular)]
    rows=[]
    for kind,Y in targets:
        R=PolynomialRing(QQ,'s'); s=R.gen()
        delta=(matrix.identity(R,4)-s*Y).det()
        for label,F in [('power_1',d),('power_3',d**3),('removed_trivial_square',(d-1)**2)]:
            q=int(F.degree())
            H=sum(k*F[k]*K**(k-1)*delta**(q-k) for k in range(1,q+1))
            P=F.derivative()(K)*delta**q-2*s*H
            algebra=R.quotient(P,'sbar'); sbar=algebra.gen()
            I=matrix.identity(algebra,4)
            g=(I+sbar*Y)*(I-sbar*Y).inverse()
            exterior=[]
            for k in range(5):
                sets=list(combinations(range(4),k))
                exterior.append(matrix(algebra,[[g.matrix_from_rows_and_columns(list(U),list(V)).det()
                    for V in sets] for U in sets]))
            chi=sum(M.trace() for M in exterior)
            v=vector(algebra,[sum((M*D).trace() for M,D in zip(exterior,ds)) for ds in derivatives])
            coordinates=actual_G.change_ring(algebra).inverse()*v*F.derivative()(chi)/F.derivative()(K)
            image=sum((c*H for c,H in zip(coordinates,basis)),matrix(algebra,4,4))
            assert image==Y.change_ring(algebra)
            if kind=='invertible_nonsemisimple': assert P.degree()==4*q
            rows.append(dict(target=kind,character=label,fiber_algebra_dimension=int(P.degree()),
                             full_exterior_matrix_trace_verified=True))
    return rows


def main():
    scalar=scalar_checks()
    trace=trace_checks()
    original=original_fiber_checks()
    singular=singular_target_checks()
    matrices=full_matrix_fiber_checks()
    paths=['SYMPLECTIC_CAYLEY_DEGREES.md','tools/verify_symplectic_cayley_degree.py',
           'tools/monomial_staircase.py','tools/verify_mixed_regularity.py',
           'SYMPLECTIC_HARMONIC_LOCALIZATION.md']
    report=dict(result='passed',date='2026-10-22',sage_version=SAGE_VERSION,
        theorem='deg(theta_F)=2*n*deg(F) for actual nontrivial characters F(det(I+g)) on Sp_(2n)(C)',
        tensor_power_degree='2*n*m',polynomial_character_extension=True,
        surjectivity_proved=True,finite_flat_over_invertible_targets=True,
        scalar_fiber_cases=len(scalar),independent_weight_fiber_cases=len(trace),
        original_Laurent_Groebner_fibers=len(original),singular_target_degree_cases=len(singular),
        full_exterior_matrix_fiber_cases=len(matrices),
        scalar_fibers=scalar,independent_trace_fibers=trace,original_fibers=original,
        singular_targets=singular,full_matrix_fibers=matrices,novelty_claim=False,externally_reviewed=False,
        runtime_API_changed=False,
        scope='Finite exact corroboration; all-rank degree proof is by a rational function on each generic radial line. No minimum-Cayley-degree or general-representation classification.',
        source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'symplectic-cayley-degrees-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in
        ('source_sha256','scalar_fibers','independent_trace_fibers','original_fibers','singular_targets',
         'full_matrix_fibers')},indent=2))


if __name__=='__main__': main()
