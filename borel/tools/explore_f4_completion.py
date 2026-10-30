"""Inventory full F4 adjoint characters with coroot convention explicit."""
from collections import Counter,defaultdict
from pathlib import Path
import json
from sage.all import QQ,ZZ,RootSystem,matrix,vector,lcm,gcd,prod
from verify_affine_branching_classification import levi_data,subsets
from verify_quadratic_branching_fibers import convolve,fundamental_local
ROOT=Path(__file__).resolve().parents[1]

def restriction(I):
    ambient=RootSystem('F4').ambient_space();cr=ambient.simple_coroots();sr=ambient.simple_roots()
    A=matrix(QQ,[[sr[i+1].scalar(cr[j+1]) for j in range(4)] for i in range(4)])
    central=[]
    for j in range(4):
        if j in I:continue
        v=A.inverse().column(j);v*=lcm(x.denominator() for x in v);v/=gcd(list(v));central.append(v)
    groups=defaultdict(Counter)
    weights=Counter({tuple(int(a.scalar(cr[j+1])) for j in range(4)):1 for a in ambient.roots()});weights[(0,)*4]=4
    for w,m in weights.items():
        q=tuple(int(vector(QQ,w).dot_product(v)) for v in central)
        groups[q][tuple(w[i] for i in I)]+=m
    components,height=levi_data('F4',I);terms=[]
    for q,original in sorted(groups.items()):
        residual=original.copy()
        while residual:
            dominant=[w for w in residual if all(v>=0 for v in w)];assert dominant
            nu=max(dominant,key=lambda w:(sum(h*v for h,v in zip(height,w)),w));c=residual[nu]
            small={(0,)*len(I):1}
            for i,e in enumerate(nu):
                for _ in range(e):small=convolve(small,fundamental_local('F4',I,i))
            for w,m in small.items():
                residual[w]-=c*m
                if not residual[w]:del residual[w]
            if c:terms.append(dict(exponent=list(nu),charge=list(q),coefficient=int(c)))
    return dict(type='F4',nodes=[i+1 for i in I],factors=[f for f,o in components],
                factor_orders=[list(o) for f,o in components],
                central_cocharacters=[list(map(int,v)) for v in central],terms=terms,
                degree=max(sum(v['exponent']) for v in terms))

if __name__=='__main__':
    rows=[restriction(I) for I in subsets(4)]
    for r in rows:print(r['nodes'],r['factors'],'degree',r['degree'],flush=True)
    (ROOT/'f4-adjoint-case-exploration.json').write_text(json.dumps(dict(status='exploratory',cases=rows,generic_smoothness_claim=False),indent=2)+'\n')
