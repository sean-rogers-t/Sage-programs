"""Independent arithmetic and witness audit; imports neither Sage nor builder."""
from pathlib import Path
from hashlib import sha256
from math import comb,gcd
from functools import reduce
import json
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

ROOT=Path(__file__).resolve().parents[1]
def digest(p):return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def valuation(n,p):
    e=0
    while n%p==0:e+=1;n//=p
    return e
def main():
    reports=json.loads((ROOT/'type-a-k-annihilator-witnesses.json').read_text())
    for name,expected in reports['source_sha256'].items():assert digest(ROOT/name)==expected
    finite=[]
    for row in reports['cases']:
        m=row['m'];d=json.loads((ROOT/('experiments/type-a-k-specialization-%d.json'%m)).read_text())
        variables={f'c{i}':sp.Symbol(f'c{i}') for i in range(1,m)}
        values={variables[f'c{i}']:comb(m,i) for i in range(1,m)}
        ev=lambda x:int(sp.sympify(x.replace('^','**'),locals=variables).subs(values))
        B=sp.Matrix([[ev(x) for x in r] for r in d['presentation']])
        R=sp.Matrix([[ev(x) for x in r] for r in d['syzygies']])
        k=row['kernel_rank']
        K=sp.Matrix(row['K']) if k else sp.zeros(B.cols,0)
        W=sp.Matrix(row['W']) if k else sp.zeros(R.cols,0)
        assert B*K==sp.zeros(B.rows,k)
        assert B.rank()==m and k==B.cols-m
        assert all(x.is_Integer for x in K) and all(x.is_Integer for x in W)
        D=smith_normal_form(K,domain=sp.ZZ)
        assert all(abs(D[i,i])==1 for i in range(k))
        assert R*W==m*K
        finite.append(dict(m=m,primitive_full_kernel_basis=True,explicit_annihilation_identity=True,kernel_rank=k))
    local_cases=0;quadratics=0
    for m in range(2,65):
        contents=[1]
        for j in range(1,m):
            contents.append(reduce(gcd,(comb(m,i)*contents[j-i] for i in range(1,j+1))))
            assert contents[j]==m//gcd(m,j)
        for p,s in sp.factorint(m).items():
            b=p**s;w=[0]+[max(0,s-valuation(j,p)) for j in range(1,m)]
            for j in range(1,m):
                q,r=divmod(j,b)
                coefficient=(comb(m,b)**q if q else 1)*(comb(m,r) if r else 1)
                assert valuation(coefficient,p)==w[j]
                assert comb(m,j)%(p**w[j])==0
            for i in range(1,m):
                for j in range(i,m):
                    if i+j<m:
                        assert w[i]+w[j]-w[i+j]>=0
                        assert p**w[i]*p**w[j]==p**(w[i]+w[j]-w[i+j])*p**w[i+j]
                    quadratics+=1
            local_cases+=1
    a,L=sp.symbols('a L')
    assert sp.expand((a+3)**2*L-(3+2*a)*3*L-a*a*L)==0
    assert 3*a*a*L==a*a*(3*L)
    c=sp.Symbol('c')
    assert sp.rem(2*c,c*c,domain=sp.GF(3))!=0
    n,c1,c2,c3=sp.symbols('n c1 c2 c3')
    f1=n*L-c1
    f2=n*(n-1)*L**2/2-(n-1)*c1*L+c2
    f3=n*(n-1)*(n-2)*L**3/6-(n-1)*(n-2)*c1*L**2/2+(n-2)*c2*L-c3
    C2=(n-1)*c1**2-2*n*c2
    D3=-c1*c2+3*n*c3/(n-2)
    assert sp.cancel(f1**2-2*n*f2/(n-1)-C2/(n-1))==0
    assert sp.cancel(f1*f2-3*n*f3/(n-2)-C2*L-D3)==0
    assert sp.cancel(n*(f1*f2-3*n*f3/(n-2))-C2*f1-C2*c1-n*D3)==0
    aug={c1:n,c2:n*(n-1)/2,c3:n*(n-1)*(n-2)/6}
    assert sp.cancel(C2.subs(aug))==0 and sp.cancel(D3.subs(aug))==0
    paths=['TYPE_A_K_ANNIHILATOR.md','tools/audit_type_a_k_annihilator.py','type-a-k-annihilator-witnesses.json']
    out=dict(result='passed',full_image_cases=finite,ordinary_content_ranks=63,
             local_basis_checks=local_cases,quadratic_table_checks=quadratics,
             cubic_second_order_control=True,two_uniform_low_degree_lifts_symbolically_verified=True,
             arbitrary_rank_annihilator_proved=False,
             source_sha256={name:digest(ROOT/name) for name in paths},
             historical_method_date='2026-10-27',superseded_by='TYPE_A_K_INTEGRAL_RECURRENCE.md',
             scope='Independent finite integer witnesses and arithmetic checks. This earlier method does not prove the uniform annihilator; the subsequent trace recurrence does. The negative flag records this method only.')
    (ROOT/'type-a-k-annihilator-audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('full_image_cases','source_sha256')}))
if __name__=='__main__':main()
