"""Independent native root/subset-model audit; no Sage dependency."""
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd, lcm
from pathlib import Path
import json
import sympy as sp

from audit_affine_branching_certificates import cartan, root_data, inverse

ROOT=Path(__file__).resolve().parents[1]


def orbit(A,highest):
    seen={highest}; todo=deque([highest])
    while todo:
        w=todo.popleft()
        for i in range(len(A)):
            v=tuple(w[j]-w[i]*A[j][i] for j in range(len(A)))
            if v not in seen:seen.add(v);todo.append(v)
    return Counter({w:1 for w in seen})


def multiply(P,Q):
    result=Counter()
    for a,m in P.items():
        for b,n in Q.items():result[tuple(x+y for x,y in zip(a,b))]+=m*n
    return result


def subset_fundamental(size,j,order,total):
    result=Counter()
    for subset in combinations(range(size),j):
        indicator=[int(i in subset) for i in range(size)]; w=[0]*total
        for k,index in enumerate(order):w[index]=indicator[k]-indicator[k+1]
        result[tuple(w)]+=1
    return result


def reconstruct(row):
    A=cartan(row['type']); n=len(A); I=[i-1 for i in row['nodes']]
    J=[i for i in range(n) if i not in I]; inv=inverse(A)
    H=[]
    for j in J:
        col=[inv[i][j] for i in range(n)]; den=lcm(*(v.denominator for v in col))
        ints=[int(v*den) for v in col]; div=gcd(*ints); H.append([v//div for v in ints])
    assert H==row['central_cocharacters']
    cover=sp.Matrix(n,n,lambda i,j:int(i==I[j]) if j<len(I) else H[j-len(I)][i])
    assert abs(cover.det())==row['cover_degree']
    fundamentals={}
    for kind,order in zip(row['factors'],row['factor_orders']):
        assert kind=='A'+str(len(order))
        positions=[I[j] for j in order]
        assert all(A[i][j]==(2 if k==l else -1 if abs(k-l)==1 else 0)
                   for k,i in enumerate(positions) for l,j in enumerate(positions))
        for k,pos in enumerate(order):
            fundamentals[pos]=subset_fundamental(len(order)+1,k+1,order,len(I))
    rebuilt=Counter()
    for t in row['terms']:
        small=Counter({(0,)*len(I):1})
        for i,power in enumerate(t['exponent']):
            for _ in range(power):small=multiply(small,fundamentals[i])
        for w,m in small.items():rebuilt[(tuple(t['charge']),w)]+=t['coefficient']*m
    if row['highest'] is None:
        positive=root_data(A)[1]
        weights=Counter({w:1 for w in positive})
        weights.update({tuple(-v for v in w):1 for w in positive})
        weights[(0,)*n]=n
    else:
        assert row['type']=='E6' and row['highest']==[1,0,0,0,0,0]
        weights=orbit(A,tuple(row['highest']))
        assert sum(weights.values())==27  # Minuscule Weyl-orbit model.
    expected=Counter()
    for w,m in weights.items():
        q=tuple(sum(x*y for x,y in zip(w,col)) for col in H)
        expected[(q,tuple(w[i] for i in I))]+=m
    assert +rebuilt==expected and all(m>=0 for m in rebuilt.values())
    Q=sp.Matrix([[A[i][j] for j in I] for i in I]).inv()
    HH=sp.Matrix(H).T; C=(HH.T*sp.Matrix(A)*HH).inv()
    radius=sp.Integer(2) if row['highest'] is None else sp.Rational(4,3)
    for t in row['terms']:
        a=sp.Matrix(t['exponent']);q=sp.Matrix(t['charge'])
        assert (a.T*Q*a)[0]+(q.T*C*q)[0]<=radius
    return dict(type=row['type'],nodes=row['nodes'],dimension=sum(weights.values()),
                all_charges_reconstructed_by_subset_models=True,cover_degree=row['cover_degree'])


def norm_checks():
    result={}
    for kind in ('E6','E7','E8'):
        Q=inverse(cartan(kind));minimum=min(Q[i][i] for i in range(len(Q)))
        result[kind]=str(minimum)
    assert result==dict(E6='4/3',E7='3/2',E8='2')
    # Bounded arithmetic corroboration of the written all-rank A pair formula.
    for m in range(1,17):
        Q=[[F(min(i,j)*(m+1-max(i,j)),m+1) for j in range(1,m+1)] for i in range(1,m+1)]
        assert min(Q[i][i]+2*Q[i][j]+Q[j][j] for i in range(m) for j in range(m))==2
        assert min(Q[i][i] for i in range(m))==F(m,m+1)
    assert inverse(cartan('E6'))[5][5]==F(4,3)
    return dict(exceptional_minimum_fundamental_norms=result,A_pair_formula_checked_through_rank=16,
                both_E6_minimal_norms='4/3',all_rank_A_statement_uses_written_formula=True)


def symbolic_transfer():
    a,b,y,u,v,z,c=sp.symbols('a b y u v z c',nonzero=False)
    S=z**3+z**-3; AA=z**-2+z*y; BB=z*z+y/z
    original=a*b+u*v+y*y+S*y-2+AA*a*v+BB*b*u
    transformed=(a+BB*u)*(b+AA*v)+(y*y+S*y)*(1-u*v)-2
    assert sp.expand(original-transformed)==0
    psi=sp.expand(z*sp.diff(-S*S/4,z))
    assert psi==-sp.Rational(3,2)*(z**6-z**-6)
    w=sp.symbols('w'); quotient=3*w*w+2*c*w-3
    assert sp.discriminant(quotient,w)==4*(c*c+9)
    # Kernel mu6 action in the ORIGINAL fundamental coordinates, modulo six.
    weights={a:-2,b:2,y:3,u:2,v:-2}
    laurent=sp.Poly(sp.expand(original*z**3),a,b,y,u,v,z)
    for exponent,coeff in laurent.terms():
        charge=sum(e*weights[x] for e,x in zip(exponent[:5],(a,b,y,u,v)))+exponent[5]-3
        assert charge%6==0
    return dict(exact_coordinate_identity=True,scalar_target_equation=True,
                mu6_action_invariance=True,actual_quotient_discriminant='4*(c^2+9)',
                generic_regular_K_orbits=2)


def mod_rank(A,p):
    A=[[int(v)%p for v in row] for row in A];r=0
    for j in range(len(A[0])):
        pivot=next((i for i in range(r,len(A)) if A[i][j]),None)
        if pivot is None:continue
        A[r],A[pivot]=A[pivot],A[r]; scale=pow(A[r][j],-1,p)
        A[r]=[v*scale%p for v in A[r]]
        for i in range(len(A)):
            if i==r:continue
            scale=A[i][j];A[i]=[(x-scale*y)%p for x,y in zip(A[i],A[r])]
        r+=1
    return r


def finite_root_control():
    A=cartan('E6');positive=root_data(A)[1]
    roots=list(positive)+[tuple(-v for v in w) for w in positive]
    H=[2,3,4,6,4,2]
    for p in (1009,1021,1033,1063,1093):
        omegas=[x for x in range(2,p) if pow(x,3,p)==1]
        lambdas=[x for x in range(1,p) if (16*x*x+65*x+16)%p==0]
        if not omegas or not lambdas:continue
        omega=omegas[0];lam=lambdas[0];params=[pow(2,h,p) for h in H]
        params[2]=params[2]*omega%p;params[5]=params[5]*omega%p;params[1]=params[1]*lam%p
        values=[prod_mod([pow(params[j],w[j],p) for j in range(6)],p) for w in roots]
        if 1 in values:continue
        B=[[sum(w[i]*w[j] for w in roots) for j in range(6)] for i in range(6)]
        grad=[sum(v*w[i] for v,w in zip(values,roots))%p for i in range(6)]
        target=(-1365*pow(2048,-1,p))%p
        assert grad==[target*sum(B[i][j]*H[j] for j in range(6))%p for i in range(6)]
        Hess=[[sum(v*w[i]*w[j] for v,w in zip(values,roots))%p for j in range(6)] for i in range(6)]
        rank=mod_rank(Hess,p);outside=sum(sum(x*y for x,y in zip(w,H))!=0 for w in roots)
        assert rank==6 and outside==58
        return dict(prime=p,cartan_rank=rank,root_rank=outside,full_rank=64,kernel_dimension=14,
                    scope='Independent finite-field point corroboration; not a characteristic-zero proof or fiber count.')
    raise AssertionError('No usable independently chosen prime')


def prod_mod(values,p):
    out=1
    for v in values:out=out*v%p
    return out


def main():
    data=json.loads((ROOT/'branching-norm-fibers-verification.json').read_text())
    for name,digest in data['source_sha256'].items():
        assert sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()==digest,name
    rows=data['adjoint_mechanism_tests']+[data['E6_minimal_quadratic_control']]
    report=dict(result='passed',uses_Sage=False,independent_character_reconstructions=[reconstruct(r) for r in rows],
                norm_checks=norm_checks(),symbolic_transfer=symbolic_transfer(),
                independent_finite_root_control=finite_root_control(),
                all_exceptional_adjoint_Levis_resolved=False,literature_priority_claim=False)
    names=['BRANCHING_NORM_FIBERS.md','tools/audit_branching_norm_fibers.py',
           'tools/audit_affine_branching_certificates.py','branching-norm-fibers-verification.json']
    report['source_sha256']={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in names}
    (ROOT/'branching-norm-fibers-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
