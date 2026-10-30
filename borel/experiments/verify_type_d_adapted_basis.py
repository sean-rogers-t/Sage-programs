"""Construct full D bases using ordinary Grassmannian rows and binary gluing.

No Smith algorithm is used to construct the basis. Fresh Cartesian D
divided differences check the ambient integral classes afterward.
"""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
from sage.all import (GF, QQ, ZZ, Partitions, Partition, matrix, block_diagonal_matrix,
                      diagonal_matrix, PolynomialRing, prod)
from verify_type_b_adapted_basis import basis, SF, S
from verify_type_d_candidate import full_candidate
from verify_type_d_square_comparison import schur_polynomials
from type_d_cartesian_frontier import CartesianDFrontier
from sage_borel import classical_borel
from sage_borel.operations import elementary_polynomial

ROOT=Path(__file__).resolve().parents[1]

@lru_cache(None)
def construct(K,H,d):
    sources=[]; blocks=[]; exps=[]
    for sector in (0,1):
        weight=d-sector*H
        if weight<0 or weight%2:continue
        if K-sector==0 and weight!=0:continue
        ps,es,V,F,fs=basis(K-sector,H,weight//2)
        sources.extend((sector,p) for p in ps);blocks.append(V);exps.extend(es)
    if not sources:return None
    V0=block_diagonal_matrix(blocks)
    targets=tuple(tuple(p) for p in Partitions(d,max_part=K,max_length=H))
    eh=S((1,)*H);mon=SF.monomial()
    polys=[]
    if targets:
        for sector,p in sources:
            doubled=sum((coeff*mon(tuple(2*x for x in nu)) for nu,coeff in mon(S(p))),mon.zero())
            polys.append(S(doubled)*(eh if sector else S.one()))
    M0=matrix(ZZ,len(targets),len(sources),lambda i,j:polys[j][Partition(targets[i])])
    coords=(M0*V0.inverse()).change_ring(ZZ)
    zero=[i for i,e in enumerate(exps) if e==0]
    positive=[i for i,e in enumerate(exps) if e>0]
    assert all(coords[i,j]==0 for i in range(coords.nrows()) for j in positive)
    reduced=coords.matrix_from_columns(zero).change_ring(GF(2)).echelon_form()
    pivots=list(reduced.pivots());rank=len(pivots)
    rows=[];weights=[]
    for i in range(rank):
        row=[ZZ(0)]*len(sources)
        for j,col in enumerate(zero):row[col]=ZZ(reduced[i,j])
        rows.append(row);weights.append(0)
    for j,col in enumerate(zero):
        if j not in pivots:
            row=[ZZ(0)]*len(sources);row[col]=1
            rows.append(row);weights.append(1)
    for col in positive:
        row=[ZZ(0)]*len(sources);row[col]=1
        rows.append(row);weights.append(exps[col])
    U=matrix(ZZ,rows);V=U*V0
    assert abs(V.det())==1
    assert Counter(weights)==full_candidate(K,H,d)
    # Check the full row lattice, not just the predicted Smith multiset.
    old=matrix(ZZ,[[2**max(1,e)*x for x in V0.row(i)] for i,e in enumerate(exps)])
    expected=M0.stack(old)
    new=diagonal_matrix(ZZ,[2**e for e in weights])*V
    assert expected.row_module()==new.row_module()
    return sources,tuple(weights),V,dict(K=K,H=H,d=d,rank=len(sources),
                zero_prefix_rank=rank,gluing_count=len(zero)-rank,
                binary_row_construction=True,full_row_modules_equal=True)

def actual_matrix(K,H,d,cache):
    if (K,H) not in cache:
        b=classical_borel('D',K+H,H)
        ring=PolynomialRing(QQ,K+H,names=['x%d'%i for i in range(K+H)])
        cache[K,H]=CartesianDFrontier(b),ring
    frontier,ring=cache[K,H];x=ring.gens();cols=[];labels=[]
    for sector in (0,1):
        remaining=d-sector*H
        if remaining<0 or remaining%2:continue
        if K-sector==0 and remaining:continue
        formal,parts,polys=schur_polynomials(K-sector,H,remaining//2)
        hom=formal.hom([elementary_polynomial(ring,[z*z for z in x[K:]],i)
                       for i in range(1,H+1)],ring)
        factor=prod(x[K:]) if sector else ring.one()
        for p,f in zip(parts,polys):
            labels.append((sector,tuple(p)));cols.append(frontier.expand(factor*hom(f)))
    words=sorted(set(w for col in cols for w in col),key=str)
    return labels,matrix(ZZ,len(words),len(cols),lambda i,j:cols[j].get(words[i],0))

def product_checks(K,H):
    classes=[]
    for d in range(2*K*H+1):
        data=construct(K,H,d)
        if data is None:continue
        source,weights,V,meta=data;F=V.inverse().change_ring(ZZ)
        for j,e in enumerate(weights):
            components=tuple(sum((F[i,j]*S(p) for i,(s,p) in enumerate(source) if s==sec),S.zero())
                             for sec in (0,1))
            classes.append((d,e,components))
    count=0;eh=S((1,)*H)
    for i,(d,e,(f,g)) in enumerate(classes):
        for D,E,(F,G) in classes[i:]:
            if d+D>2*K*H:continue
            target=construct(K,H,d+D)
            plus=f*F+eh*g*G;minus=f*G+g*F
            if target is None:
                assert not any(coeff for p,coeff in plus if len(p)<=H and (not p or p[0]<=K))
                assert not any(coeff for p,coeff in minus if len(p)<=H and (not p or p[0]<=K-1))
                continue
            sources,es,V,meta=target
            col=matrix(ZZ,len(sources),1,[(plus if s==0 else minus)[Partition(p)] for s,p in sources])
            coeff=V*col
            for j,x in enumerate(coeff.column(0)):
                assert (QQ(x)*QQ(2)**(es[j]-e-E)).denominator()==1
                count+=bool(x)
    return dict(K=K,H=H,classes=len(classes),nonzero_products=count,integral=True)

def main():
    constructed=[]; actual=[];cache={}
    for K,H in [(1,2),(1,3),(1,4),(2,2),(2,3),(2,4),(3,2),(3,3),(3,4),(4,2),(4,3),(4,4)]:
        for d in range(2*K*H+1):
            value=construct(K,H,d)
            if value is not None:constructed.append(value[3])
    for K,H in [(1,3),(2,2),(3,2),(2,3)]:
        for d in range(2*K*H+1):
            value=construct(K,H,d)
            if value is None:continue
            labels,B=actual_matrix(K,H,d,cache)
            source,weights,V,meta=value
            order={label:i for i,label in enumerate(labels)}
            B=B.matrix_from_columns([order[label] for label in source])
            primitive=(B*V.inverse()*diagonal_matrix(QQ,[QQ(1)/2**e for e in weights])).change_ring(ZZ)
            assert all(abs(x)==1 for x in primitive.elementary_divisors() if x)
            assert primitive*diagonal_matrix(ZZ,[2**e for e in weights])*V==B
            actual.append(dict(K=K,H=H,d=d,integral=True,primitive=True,
                source=[[int(s),list(map(int,p))] for s,p in source],
                B_rows=B.nrows(),B_columns=B.ncols(),
                B_entries=[[i,j,int(B[i,j])] for i in range(B.nrows()) for j in range(B.ncols()) if B[i,j]]))
        print(json.dumps(dict(K=K,H=H,actual_checked=True)),flush=True)
    products=[product_checks(K,H) for K,H in [(1,2),(1,3),(2,2),(3,2),(2,3),(2,4),(3,3)]]
    report=dict(status='passed',constructed=constructed,actual_D=actual,multiplication=products,
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    (ROOT/'type-d-adapted-basis-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status='passed',constructed=len(constructed),actual=len(actual))),flush=True)

if __name__=='__main__':main()
