"""Native integer audit: no Sage, Borel library, or previous proof checker.

Permutation determinants construct Jacobi--Trudi. Pieri chains compute
products. Stored Borel columns provide an independent embedding check.
"""
from collections import defaultdict, Counter
from functools import lru_cache
from itertools import permutations
from math import factorial, comb
from pathlib import Path
from hashlib import sha256
import json
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form

ROOT=Path(__file__).resolve().parents[1]

@lru_cache(None)
def partitions(m,A,l):
    if m==0:return ((),)
    if l==0:return ()
    return tuple((i,)+p for i in range(min(A,m),0,-1)
                 for p in partitions(m-i,i,l-1))

def parity(p):
    return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))

@lru_cache(None)
def jt(p):
    result=defaultdict(int)
    for perm in permutations(range(len(p))):
        indices=tuple(p[i]-i+perm[i] for i in range(len(p)))
        if min(indices,default=0)<0:continue
        result[tuple(sorted((x for x in indices if x),reverse=True))]+=parity(perm)
    return {p:c for p,c in result.items() if c}

def multiply_dict(a,b):
    result=defaultdict(int)
    for x,c in a.items():
        for y,d in b.items():result[tuple(sorted(x+y,reverse=True))]+=c*d
    return {p:c for p,c in result.items() if c}

@lru_cache(None)
def construct(A,l,m):
    ps=tuple(sorted(partitions(m,A,l),key=lambda p:(sum(x*x for x in p),p)))
    index={p:i for i,p in enumerate(ps)}
    r,c=(A+l)//2,(A+l+1)//2
    exps=[]; hrows=[]
    for p in ps:
        g=p+(0,)*(l-len(p)); beads=[g[i]+l-i-1 for i in range(l)]
        d=sum(x<r for x in beads);t=sum(x>=c+d for x in beads);b=l-t
        valid=[]
        for j in range(max(0,b-r),min(b,c,(b+(A+l)%2)//2)+1):
            bounds=(c-j,)*j+(r-b+j,)*(b-j)
            if all(x<=y for x,y in zip(g[t:],bounds)):valid.append(j)
        assert valid
        j=min(valid)
        factors=[tuple(x for x in z if x) for z in (g[:t],g[t:t+j],g[t+j:])]
        expansion=multiply_dict({factors[0]:1},multiply_dict(jt(factors[1]),jt(factors[2])))
        row=[0]*len(ps)
        for q,v in expansion.items():
            if q and q[0]>A:continue
            assert len(q)<=l
            row[index[q]]+=v
        hrows.append(row);exps.append(t)
    U=sp.Matrix(hrows)
    assert all(U[i,j]==0 for i in range(len(ps)) for j in range(i))
    assert all(U[i,i]==1 for i in range(len(ps)))
    K=sp.Matrix([[jt(p).get(q,0) for q in ps] for p in ps])
    V=U*K.inv();F=V.inv()
    assert all(x.q==1 for x in V) and all(x.q==1 for x in F)
    return ps,tuple(exps),V,F

@lru_cache(None)
def schur_product(A,l,x,y):
    if len(x)>l or len(y)>l or (x and x[0]>A) or (y and y[0]>A):
        return {}
    ans=defaultdict(int)
    for hp,sign in jt(y).items():
        state={x:1}
        for z in hp:
            nxt=defaultdict(int)
            for p,v in state.items():
                pp=p+(0,)*(l-len(p))
                for q in partitions(sum(p)+z,A,l):
                    qq=q+(0,)*(l-len(q))
                    if all(qq[i]>=pp[i] for i in range(l)) and all(pp[i]>=qq[i+1] for i in range(l-1)):
                        nxt[q]+=v
            state=nxt
        for p,v in state.items():ans[p]+=sign*v
    return {p:v for p,v in ans.items() if v}

def products(A,l):
    count=0
    for u in range(A*l+1):
        ps,es,V,F=construct(A,l,u)
        for v in range(u,A*l-u+1):
            qs,ds,W,G=construct(A,l,v)
            rs,ts,X,J=construct(A,l,u+v)
            raw={(i,j):sp.Matrix([schur_product(A,l,p,q).get(r,0) for r in rs])
                 for i,p in enumerate(ps) for j,q in enumerate(qs)}
            for i in range(len(ps)):
                for j in range(len(qs)):
                    z=sp.zeros(len(rs),1)
                    for a in range(len(ps)):
                        for b in range(len(qs)):
                            if F[a,i] and G[b,j]:z+=F[a,i]*G[b,j]*raw[a,b]
                    coeff=X*z
                    for k,value in enumerate(coeff):
                        scaled=value*sp.Rational(2)**(ts[k]-es[i]-ds[j])
                        assert scaled.q==1
                        count+=bool(value)
    return count

def conjugate(p):return tuple(sum(x>=i for x in p) for i in range(1,max(p,default=0)+1))

def main():
    recorded_path=ROOT/'experiments/type_b_multiplicity_search_results.json'
    recorded=json.loads(recorded_path.read_text())
    checked=[]
    for case in recorded['cases']:
        A,l=case['k'],case['a']
        for row in case['degrees']:
            if row['status']!='computed':continue
            m=row['half_polynomial_degree'];ps,es,V,F=construct(A,l,m)
            data=row['B_matrix'];B=sp.zeros(data['rows'],data['columns'])
            for z in data['entries']:B[z['row'],z['column']]=int(z['value'])
            old={conjugate(p):i for i,p in enumerate(row['partitions'])}
            head=(-1)**m*B[:,[old[p] for p in ps]]
            basis=head*F*sp.diag(*[sp.Rational(1,2**e) for e in es])
            assert all(x.q==1 for x in basis)
            basis=sp.Matrix(basis.rows,basis.cols,[int(x) for x in basis])
            assert hermite_normal_form(basis.T)==sp.eye(len(ps))
            assert basis*sp.diag(*[2**e for e in es])*V==head
            checked.append(dict(A=A,ell=l,m=m,rank=len(ps),primitive=True))
    product_checks=[]
    for A,l in [(3,2),(2,3),(4,2),(3,3),(4,3)]:
        product_checks.append(dict(A=A,ell=l,nonzero_constants=products(A,l)))
        print(json.dumps(product_checks[-1]),flush=True)
    corners=[]
    for A,l in [(2,2),(3,2),(4,2),(3,3),(4,3),(5,3),(4,4),(6,3)]:
        c=(A+l+1)//2
        family=[p for m in range(A*l+1) for p in partitions(m,A,l)
                if not p or min(p)>=c]
        for nu in family:
            ps,es,V,F=construct(A,l,sum(nu))
            # Monomial/complete Hall duality and permutation determinant,
            # independently of Sage's monomial symmetric-function conversion.
            col=sp.Matrix([jt(p).get(nu,0) for p in ps])
            assert V*col==sp.Matrix([int(p==nu) for p in ps])
            assert es[ps.index(nu)]==len(nu)
        count=0
        for i,lam in enumerate(family):
            for mu in family[i:]:
                m=sum(lam)+sum(mu)
                if m>A*l:continue
                xs=partitions(sum(lam),A,l);ys=partitions(sum(mu),A,l)
                out=defaultdict(int)
                for x in xs:
                    for y in ys:
                        weight=jt(x).get(lam,0)*jt(y).get(mu,0)
                        if weight:
                            for q,z in schur_product(A,l,x,y).items():out[q]+=weight*z
                nu=tuple(sorted(lam+mu,reverse=True))
                coefficient=1
                for j in range(c,A+1):coefficient*=comb(lam.count(j)+mu.count(j),lam.count(j))
                for p in partitions(m,A,l):
                    expected=coefficient*jt(p).get(nu,0) if len(nu)<=l else 0
                    assert out[p]==expected
                count+=1
        corners.append(dict(A=A,ell=l,classes=len(family),products=count))
    report=dict(status='passed',method='Native permutation Jacobi--Trudi, Pieri chains and integer Hermite; no Sage.',
                embeddings=checked,products=product_checks,
                divided_power_corners=corners,
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                input_sha256=sha256(recorded_path.read_bytes()).hexdigest())
    (ROOT/'type-b-adapted-basis-audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',embedding_cases=len(checked))),flush=True)

if __name__=='__main__':main()
