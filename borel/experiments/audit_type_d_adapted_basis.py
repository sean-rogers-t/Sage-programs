"""Native full-D replay using Frobenius/Pieri and binary Gaussian elimination.

Unlike the Sage check, no plethysm is computed: squares in characteristic
two determine the exact binary gluing data. Actual exported D columns
are checked by integer Hermite, not their recorded Smith factors.
"""
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import json
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form
from audit_type_b_adapted_basis import construct as sector_basis, partitions, schur_product

ROOT=Path(__file__).resolve().parents[1]

def binary_rref(rows,ncols):
    rows=[[int(x)%2 for x in row] for row in rows]
    k=0;pivots=[]
    for j in range(ncols):
        p=next((i for i in range(k,len(rows)) if rows[i][j]),None)
        if p is None:continue
        rows[k],rows[p]=rows[p],rows[k]
        for i in range(len(rows)):
            if i!=k and rows[i][j]:rows[i]=[a^b for a,b in zip(rows[i],rows[k])]
        pivots.append(j);k+=1
    return rows[:k],pivots

@lru_cache(None)
def construct(K,H,d):
    source=[];blocks=[];exps=[]
    for sec in (0,1):
        rem=d-sec*H
        if rem<0 or rem%2 or (K-sec==0 and rem):continue
        ps,es,V,F=sector_basis(K-sec,H,rem//2)
        source.extend((sec,p) for p in ps);blocks.append(V);exps.extend(es)
    V0=sp.diag(*blocks)
    targets=partitions(d,K,H)
    cols=[]
    for sec,p in source:
        square=schur_product(K-sec,H,p,p)
        if sec:
            square={tuple((q+(0,)*(H-len(q)))[i]+1 for i in range(H)):v
                    for q,v in square.items()}
        cols.append([square.get(t,0)%2 for t in targets])
    M0=sp.Matrix(len(targets),len(source),lambda i,j:cols[j][i])
    coords=M0*V0.inv()
    zero=[i for i,e in enumerate(exps) if not e];positive=[i for i,e in enumerate(exps) if e]
    assert all(int(coords[i,j])%2==0 for i in range(coords.rows) for j in positive)
    reduced,pivots=binary_rref([[coords[i,j] for j in zero] for i in range(coords.rows)],len(zero))
    rows=[];weights=[]
    for r in reduced:
        row=[0]*len(source)
        for j,k in enumerate(zero):row[k]=r[j]
        rows.append(row);weights.append(0)
    for j,k in enumerate(zero):
        if j not in pivots:
            row=[0]*len(source);row[k]=1;rows.append(row);weights.append(1)
    for k in positive:
        row=[0]*len(source);row[k]=1;rows.append(row);weights.append(exps[k])
    V=sp.Matrix(rows)*V0
    assert abs(V.det())==1
    return source,weights,V,len(zero)-len(pivots)

def main():
    path=ROOT/'type-d-adapted-basis-verification.json'
    report=json.loads(path.read_text())
    matches=0
    for row in report['constructed']:
        source,weights,V,g=construct(row['K'],row['H'],row['d'])
        assert len(source)==row['rank'] and g==row['gluing_count']
        assert weights.count(0)==row['zero_prefix_rank']
        matches+=1
    actual=[]
    for row in report['actual_D']:
        source,weights,V,g=construct(row['K'],row['H'],row['d'])
        assert [[s,list(p)] for s,p in source]==row['source']
        B=sp.zeros(row['B_rows'],row['B_columns'])
        for i,j,x in row['B_entries']:B[i,j]=x
        C=B*V.inv()*sp.diag(*[sp.Rational(1,2**e) for e in weights])
        assert all(x.q==1 for x in C)
        C=sp.Matrix(C.rows,C.cols,[int(x) for x in C])
        assert hermite_normal_form(C.T)==sp.eye(len(source))
        assert C*sp.diag(*[2**e for e in weights])*V==B
        actual.append(dict(K=row['K'],H=row['H'],d=row['d'],primitive=True))
    products=[]
    for item in report['multiplication']:
        K,H=item['K'],item['H'];classes=[];count=0
        for d in range(2*K*H+1):
            src,es,V,g=construct(K,H,d)
            if not src:continue
            F=V.inv()
            for j,e in enumerate(es):
                classes.append((d,e,{s:int(F[i,j]) for i,s in enumerate(src) if F[i,j]}))
        for i,(d,e,f) in enumerate(classes):
            for D,E,g in classes[i:]:
                if d+D>2*K*H:continue
                out={}
                for (s,p),a in f.items():
                    for (t,q),b in g.items():
                        sec=(s+t)%2
                        for lam,v in schur_product(K-sec,H,p,q).items():
                            if s+t==2:
                                lam=tuple(x+1 for x in lam+(0,)*(H-len(lam)))
                                if lam[0]>K:continue
                            label=(sec,lam);out[label]=out.get(label,0)+a*b*v
                src,es,V,glue=construct(K,H,d+D)
                if not src:
                    assert not any(out.values()), (K,H,d,D,f,g,out)
                    continue
                col=V*sp.Matrix([out.get(label,0) for label in src])
                for j,x in enumerate(col):
                    assert (x*sp.Rational(2)**(es[j]-e-E)).q==1
                    count+=bool(x)
        assert len(classes)==item['classes'] and count==item['nonzero_products']
        products.append(dict(K=K,H=H,classes=len(classes),nonzero_products=count,integral=True))
    result=dict(status='passed',method='Native Pieri/Frobenius, binary row reduction and integer Hermite; no Sage.',
                constructed_matches=matches,actual_D=actual,multiplication=products,
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                input_sha256=sha256(path.read_bytes()).hexdigest())
    (ROOT/'type-d-adapted-basis-audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',constructed=matches,actual=len(actual))),flush=True)

if __name__=='__main__':main()
