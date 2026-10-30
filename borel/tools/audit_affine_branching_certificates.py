"""Sage-free audit of every exceptional affinity certificate."""
from fractions import Fraction as F
from functools import lru_cache
from collections import deque
from pathlib import Path
from hashlib import sha256
import json

ROOT=Path(__file__).resolve().parents[1]


def transpose(A): return tuple(zip(*A))


def multiply(A,B):
    return tuple(tuple(sum(x*y for x,y in zip(row,col)) for col in transpose(B)) for row in A)


def inverse(A):
    n=len(A); rows=[[F(x) for x in row]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for j in range(n):
        pivot=next(i for i in range(j,n) if rows[i][j])
        rows[j],rows[pivot]=rows[pivot],rows[j]
        scale=rows[j][j]; rows[j]=[x/scale for x in rows[j]]
        for i in range(n):
            if i==j: continue
            scale=rows[i][j]; rows[i]=[x-scale*y for x,y in zip(rows[i],rows[j])]
    return tuple(tuple(row[n:]) for row in rows)


def cartan(kind):
    n=int(kind[1:]); A=[[2*int(i==j) for j in range(n)] for i in range(n)]
    if kind=='G2': A[0][1]=-3; A[1][0]=-1
    elif kind=='F4':
        for i in range(3): A[i][i+1]=A[i+1][i]=-1
        A[2][1]=-2
    else:
        assert kind in ('E6','E7','E8')
        for i,j in [(1,3),(3,4),(2,4)]+[(i,i+1) for i in range(4,n)]:
            A[i-1][j-1]=A[j-1][i-1]=-1
    return tuple(tuple(row) for row in A)


def root_gram(A):
    n=len(A); d=[None]*n
    for start in range(n):
        if d[start] is not None: continue
        d[start]=F(1); todo=[start]
        while todo:
            i=todo.pop()
            for j in range(n):
                if i==j or not A[i][j]: continue
                value=d[i]*F(A[i][j],A[j][i])
                if d[j] is None: d[j]=value; todo.append(j)
                else: assert d[j]==value
    B=tuple(tuple(d[i]*A[i][j] for j in range(n)) for i in range(n))
    assert B==transpose(B)
    return B


def submatrix(A,I): return tuple(tuple(A[i][j] for j in I) for i in I)


@lru_cache(None)
def root_data(A):
    n=len(A); ai=inverse(A); Q=multiply(transpose(ai),multiply(root_gram(A),ai))
    columns=list(transpose(A)); roots=set(columns); queue=deque(columns)
    while queue:
        w=queue.popleft()
        for i,col in enumerate(columns):
            v=tuple(w[j]-w[i]*col[j] for j in range(n))
            if v not in roots: roots.add(v); queue.append(v)
        assert len(roots)<=240
    positives=[]
    for root in roots:
        coeff=tuple(sum(x*y for x,y in zip(row,root)) for row in ai)
        assert all(x>=0 for x in coeff) or all(x<=0 for x in coeff)
        if all(x>=0 for x in coeff): positives.append(root)
    return Q,tuple(positives)


@lru_cache(None)
def dimension(A,highest):
    if not A: return 1
    Q,positive=root_data(A); total=F(1)
    for root in positive:
        qr=tuple(sum(x*y for x,y in zip(row,root)) for row in Q)
        total*=sum((x+1)*y for x,y in zip(highest,qr))/sum(qr)
    assert total.denominator==1
    return int(total)


def main():
    data=json.loads((ROOT/'affine-branching-classification-verification.json').read_text())
    assert data['result']=='passed'; negative=positive=words=0; seen_pairs=set()
    for record in data['exceptional_certificates']:
        kind=record['type']; k=record['node']; A=cartan(kind); n=len(A)
        pair=(kind,k); assert pair not in seen_pairs; seen_pairs.add(pair)
        highest=tuple(int(i==k-1) for i in range(n)); Q,_=root_data(A)
        assert record['dimension']==dimension(A,highest)
        bound=Q[k-1][k-1]; covered=set()
        for cert in record['negative_certificates']:
            w=highest
            for node in cert['weyl_word']:
                i=node-1; w=tuple(w[j]-w[i]*A[j][i] for j in range(n))
            assert list(w)==cert['ambient_extremal_weight']
            assert sum(w[i]*Q[i][j]*w[j] for i in range(n) for j in range(n))==bound
            words+=1
            for nodes in cert['covered_levi_node_sets']:
                I=tuple(j-1 for j in nodes)
                assert I not in covered; covered.add(I)
                assert all(w[i]>=0 for i in I) and sum(w[i] for i in I)>1
                negative+=1
        B=root_gram(A)
        for cert in record['positive_restrictions']:
            I=tuple(j-1 for j in cert['levi_nodes'])
            assert I not in covered; covered.add(I); positive+=1
            if I:
                M=submatrix(A,I); mi=inverse(M)
                gram=multiply(transpose(mi),multiply(submatrix(B,I),mi))
                assert all(x>=0 for row in gram for x in row)
                gap=min(gram[i][i]+2*gram[i][j]+gram[j][j] for i in range(len(I)) for j in range(i,len(I)))
                assert bound<gap
                assert F(cert['norm_gap']['ambient_highest_norm_squared'])==bound
                assert F(cert['norm_gap']['minimum_nonfundamental_norm_squared'])==gap
                assert [[F(x) for x in row] for row in cert['norm_gap']['levi_fundamental_Gram']]==[list(row) for row in gram]
                dim=sum(c['multiplicity']*dimension(M,tuple(c['highest'])) for c in cert['constituents'])
                assert dim==record['dimension']
                assert all(sum(c['highest'])<=1 for c in cert['constituents'])
            else:
                assert cert['constituents']==[dict(highest=[],multiplicity=record['dimension'])]
        expected={tuple(i for i in range(n) if mask>>i&1) for mask in range(2**n-1)}
        assert covered==expected
    assert seen_pairs=={(kind,k) for kind in ('G2','F4','E6','E7','E8') for k in range(1,int(kind[1:])+1)}
    assert negative==3240 and positive==133
    report=dict(result='passed',date='2026-10-23',uses_Sage=False,
                exceptional_pairs=negative+positive,negative_pairs=negative,affine_pairs=positive,
                reflection_word_certificates=words,cartan_matrices_independently_specified=True,
                norms_and_Weyl_dimensions_recomputed_with_fractions=True,
                source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for name in ['AFFINE_BRANCHING_CLASSIFICATION.md',
                                            'affine-branching-classification-verification.json',
                                            'tools/audit_affine_branching_certificates.py']})
    (ROOT/'affine-branching-classification-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__': main()
