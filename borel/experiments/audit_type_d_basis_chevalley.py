"""Independent actual D embedding via native signed permutations/Chevalley.

No Sage, BGG operators, root-system library or imported basis constructor.
The same source polynomials are evaluated by divisor multiplication on
the complete flag. Compare whole coefficient rows up to row permutation.
"""
from collections import Counter,defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import combinations_with_replacement
from pathlib import Path
import json
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]

class DFlag:
    def __init__(self,n):
        self.n=n
        self.identity=tuple(range(1,n+1))
        self.roots=tuple((i,j,s) for i in range(n) for j in range(i+1,n) for s in (-1,1))

    @lru_cache(None)
    def length(self,w):
        count=0
        for i,j,s in self.roots:
            # First nonzero Cartesian root coefficient determines positivity.
            a,b=w[i],w[j]
            coeff=(1 if a>0 else -1) if abs(a)<abs(b) else s*(1 if b>0 else -1)
            count+=coeff<0
        return count

    @staticmethod
    def reflect(w,i,j,s):
        v=list(w)
        v[i],v[j]=-s*w[j],-s*w[i]
        return tuple(v)

    @lru_cache(None)
    def divisor_edges(self,w,i):
        d=self.length(w);out=[]
        for a,b,s in self.roots:
            if i not in (a,b):continue
            v=self.reflect(w,a,b,s)
            if self.length(v)==d+1:out.append((v,1 if i==a else s))
        return tuple(out)

    @lru_cache(None)
    def monomial(self,powers):
        if not any(powers):return {self.identity:1}
        i=next(i for i,e in enumerate(powers) if e)
        prev=list(powers);prev[i]-=1
        out=defaultdict(int)
        for w,c in self.monomial(tuple(prev)).items():
            for v,a in self.divisor_edges(w,i):out[v]+=a*c
        return {w:c for w,c in out.items() if c}

    def polynomial(self,poly,variables):
        out=defaultdict(int)
        for powers,c in sp.Poly(poly,*variables).terms():
            for w,a in self.monomial(powers).items():out[w]+=int(c)*a
        return {w:c for w,c in out.items() if c}

def schur(lam,ys):
    if not lam:return sp.Integer(1)
    @lru_cache(None)
    def h(d):
        if d<0:return sp.Integer(0)
        if d==0:return sp.Integer(1)
        return sum((sp.prod(ys[i] for i in choice)
                    for choice in combinations_with_replacement(range(len(ys)),d)),sp.Integer(0))
    return sp.expand(sp.Matrix(len(lam),len(lam),lambda i,j:h(lam[i]-i+j)).det())

def main():
    path=ROOT/'type-d-adapted-basis-verification.json'
    data=json.loads(path.read_text());engines={};checks=[]
    for row in data['actual_D']:
        K,H,d=row['K'],row['H'],row['d'];n=K+H
        if n not in engines:engines[n]=DFlag(n)
        engine=engines[n];x=sp.symbols('x0:'+str(n));ys=tuple(t*t for t in x[K:])
        cols=[]
        for sector,lam in row['source']:
            f=schur(tuple(lam),ys)*(sp.prod(x[K:]) if sector else 1)
            cols.append(engine.polynomial(f,x))
        words=set().union(*(set(col) for col in cols))
        native=Counter(tuple(col.get(w,0) for col in cols) for w in words)
        expected=[[0]*row['B_columns'] for _ in range(row['B_rows'])]
        for i,j,a in row['B_entries']:expected[i][j]=a
        expected=Counter(tuple(r) for r in expected if any(r))
        assert native==expected,(K,H,d,native-expected,expected-native)
        checks.append(dict(K=K,H=H,d=d,source_columns=len(cols),coefficient_rows=len(words),exact_match=True))
    result=dict(result='passed',method='Native signed-permutation positive roots and full-flag Chevalley divisor multiplication; exact row multisets agree with BGG.',
        actual_D=checks,independent_actual_embedding=True,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256=sha256(path.read_bytes()).hexdigest())
    (ROOT/'type-d-basis-chevalley-audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(result='passed',actual_D=len(checks))),flush=True)

if __name__=='__main__':main()
