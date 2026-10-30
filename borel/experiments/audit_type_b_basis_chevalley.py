"""Recompute recorded B embeddings by a separate native Chevalley engine.

This targets algorithmic independence, not a new rank census. The D audit
provides generic polynomial/divisor machinery; B adds its short coroots.
"""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import sympy as sp
from audit_type_d_basis_chevalley import DFlag,schur

ROOT=Path(__file__).resolve().parents[1]

class BFlag(DFlag):
    @lru_cache(None)
    def length(self,w):
        return super().length(w)+sum(a<0 for a in w)

    @lru_cache(None)
    def divisor_edges(self,w,i):
        out=list(super().divisor_edges(w,i))
        v=list(w);v[i]=-v[i];v=tuple(v)
        if self.length(v)==self.length(w)+1:out.append((v,2))
        return tuple(out)

def main():
    path=ROOT/'experiments/type_b_multiplicity_search_results.json'
    data=json.loads(path.read_text());checks=[]
    for case in data['cases']:
        n,a,k=case['n'],case['a'],case['k']
        if n>6:continue
        engine=BFlag(n);x=sp.symbols('x0:'+str(n));ys=tuple(z*z for z in x[a:])
        for row in case['degrees']:
            if row['status']!='computed':continue
            cols=[engine.polynomial(schur(tuple(lam),ys),x) for lam in row['partitions']]
            words=set().union(*(set(col) for col in cols))
            native=Counter(tuple(col.get(w,0) for col in cols) for w in words)
            mat=row['B_matrix'];expected=[[0]*mat['columns'] for _ in range(mat['rows'])]
            for ent in mat['entries']:expected[ent['row']][ent['column']]=int(ent['value'])
            expected=Counter(tuple(r) for r in expected if any(r))
            assert native==expected,(n,a,k,row['polynomial_degree'],native-expected,expected-native)
            result=dict(n=n,a=a,k=k,polynomial_degree=row['polynomial_degree'],
                        source_columns=len(cols),coefficient_rows=len(words),exact_match=True)
            checks.append(result);print(json.dumps(result),flush=True)
        engine.monomial.cache_clear();engine.length.cache_clear();engine.divisor_edges.cache_clear()
        DFlag.length.cache_clear();DFlag.divisor_edges.cache_clear()
    result=dict(result='passed',method='Native B signed permutations and Chevalley divisor multiplication, including short coroot coefficient 2; full coefficient rows compared with recorded BGG.',
        actual_B=checks,independent_actual_embedding=True,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256=sha256(path.read_bytes()).hexdigest())
    (ROOT/'type-b-basis-chevalley-audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(result='passed',actual_B=len(checks))),flush=True)

if __name__=='__main__':main()
