"""Finite exact checks of the type D divided-difference special-class and Euler lemmas.

This verifies complete Schubert expansions, including omitted coefficients.
Run from WSL with sage -python check_d.py.
"""
from pathlib import Path
import sys, json
from sage.all import prod, QQ

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'papers'/'table_audit'))
from recompute import Calculator

def complete(xs,d):
    values=[xs[0].parent()(1)]+[xs[0].parent()(0)]*d
    for x in xs:
        for p in range(1,d+1):
            values[p]+=x*values[p-1]
    return values[d]

def expected(calc,n,r,k,p):
    if p==0:return {calc.W.one():QQ(1)}
    if p<k:
        return {calc.W.from_reduced_word(list(range(r+p-1,r-1,-1))):QQ(1)}
    if p==k:
        tail=list(range(n-2,r-1,-1))
        return {calc.W.from_reduced_word([j]+tail):QQ(1) for j in (n-1,n)}
    word=list(range(n+k-p,n-1))+[n,n-1]+list(range(n-2,r-1,-1))
    return {calc.W.from_reduced_word(word):QQ(2)}

rows=[]
for n in range(4,7):
    calc=Calculator('D'+str(n));w=calc.x
    x=[w[0]]+[w[j]-w[j-1] for j in range(1,n-2)]
    x+=[w[n-2]+w[n-1]-w[n-3],w[n-1]-w[n-2]]
    for k in range(2,n):
        r=n-k
        for p in range(n+k):
            got=calc.expansion(complete(x[:r],p),r)
            want=expected(calc,n,r,k,p)
            assert got==want,(n,k,p,got,want)
            rows.append({'n':n,'k':k,'degree':p,'lemma':'special','passed':True})
        got=calc.expansion(prod(x[r:]),r)
        tail=list(range(n-2,r-1,-1))
        want={calc.W.from_reduced_word([n-1]+tail):QQ((-1)**k),
              calc.W.from_reduced_word([n]+tail):QQ((-1)**(k-1))}
        assert got==want,(n,k,'Euler',got,want)
        rows.append({'n':n,'k':k,'degree':k,'lemma':'Euler','passed':True})
        print('D',n,'k',k,'all special degrees and Euler PASS',flush=True)

Path(__file__).with_name('d_results.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
print('TOTAL',len(rows),'complete expansion checks PASS',flush=True)
