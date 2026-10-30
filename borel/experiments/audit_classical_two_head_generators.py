"""Integer symbolic checks of the uniform two-head straightening proof.

All reductions below use monic integer divisors. Finite parameters
corroborate the recurrence proof, not replace its arbitrary-rank scope.
"""
from hashlib import sha256
from pathlib import Path
import json
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
u,v=sp.symbols('u v')

def sequences(maximum):
    h=[sp.Integer(1),u];p=[sp.Integer(2),u]
    for j in range(2,maximum+1):
        h.append(sp.expand(u*h[-1]-v*h[-2]))
        p.append(sp.expand(u*p[-1]-v*p[-2]))
    return h,p

def monic_reduce(f,basis):
    f=sp.Poly(f,u,v,domain=sp.ZZ);remainder=sp.Poly(0,u,v,domain=sp.ZZ)
    divisors=[sp.Poly(g,u,v,domain=sp.ZZ) for g in basis]
    assert all(g.LC()==1 for g in divisors)
    while not f.is_zero:
        power,coeff=f.terms()[0]
        for g in divisors:
            lead=g.monoms()[0]
            if all(a>=b for a,b in zip(power,lead)):
                f-=sp.Poly(coeff*u**(power[0]-lead[0])*v**(power[1]-lead[1]),u,v)*g
                break
        else:
            term=sp.Poly(coeff*u**power[0]*v**power[1],u,v)
            remainder+=term;f-=term
    return remainder.as_expr()

def main():
    checks=[]
    for n in range(4,19):
        c=(n+1)//2;h,p=sequences(3*n)
        assert sp.expand(h[2*c-1]-p[c]*h[c-1])==0
        assert sp.expand(h[2*c]-p[c]*h[c]+v**c)==0
        assert sp.expand(h[2*c-2]-p[c]*h[c-2]-v**(c-1))==0
        stop=c if n%2==0 else c-1
        basis=[sp.expand(v**j*p[c-j]) for j in range(stop)]+[v**stop]
        pairs=0
        for i,f in enumerate(basis):
            a=sp.Poly(f,u,v).monoms()[0]
            for g in basis[i+1:]:
                b=sp.Poly(g,u,v).monoms()[0];lcm=tuple(max(x,y) for x,y in zip(a,b))
                S=u**(lcm[0]-a[0])*v**(lcm[1]-a[1])*f-u**(lcm[0]-b[0])*v**(lcm[1]-b[1])*g
                assert monic_reduce(S,basis)==0
                pairs+=1
        ambient=[sp.expand(v**j*h[n-1-j]) for j in range(n)]
        products=[p[c]**2-2*v**c,p[c]*p[c+1]-u*v**c,
                  p[c+1]**2-2*v**(c+1),p[c]*v**c,p[c+1]*v**c,v**(2*c)]
        assert all(monic_reduce(f,ambient)==0 for f in products)
        hilbert={}
        for a in range(c):
            for b in range(c-a):
                if n%2 and (a,b)==(0,c-1):continue
                weight=a+2*b;hilbert[weight]=hilbert.get(weight,0)+1
        expected={}
        # U_(n,2): partitions in (c-2)^2, plus q^(c-1)[floor(n/2)].
        for a in range(c-1):
            for b in range(a+1):expected[a+b]=expected.get(a+b,0)+1
        for j in range(n//2):expected[c-1+j]=expected.get(c-1+j,0)+1
        assert hilbert==expected
        checks.append(dict(n=n,c=c,monic_S_pairs=pairs,quotient_rank=sum(hilbert.values()),
                           recurrence_identities=True,all_products=True,hilbert_matches=True))
    result=dict(result='passed',scope='Monic integer reductions and independent Hilbert counts for n=4,...,18; arbitrary rank follows from the written recurrence and adjacent-pair proof.',
        cases=checks,source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    (ROOT/'classical-two-head-generator-audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(result='passed',cases=len(checks),S_pairs=sum(c['monic_S_pairs'] for c in checks))),flush=True)

if __name__=='__main__':main()
