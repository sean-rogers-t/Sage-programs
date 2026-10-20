"""Independent E7 matrices from the minuscule weight graph.

No WeylCharacterRing, SpringerCharacterMap, or geometry interpolation code is
used here. Chevalley/Serre relations certify the 56-dimensional module.
"""
from sage.all import QQ, ZZ, matrix, diagonal_matrix, vector

C = matrix(ZZ, [[2,0,-1,0,0,0,0], [0,2,0,-1,0,0,0],
                [-1,0,2,-1,0,0,0], [0,-1,-1,2,-1,0,0],
                [0,0,0,-1,2,-1,0], [0,0,0,0,-1,2,-1],
                [0,0,0,0,0,-1,2]])


def matrices():
    found={(0,0,0,0,0,0,1)}
    pending=list(found)
    while pending:
        weight=pending.pop()
        for i in range(7):
            other=tuple(weight[j]-weight[i]*C[j,i] for j in range(7))
            if other not in found:
                found.add(other)
                pending.append(other)
    weights=tuple(sorted(found))
    assert len(weights)==56 and all(a in (-1,0,1) for w in weights for a in w)
    positions={w:i for i,w in enumerate(weights)}
    roots=tuple(tuple(C.column(i)) for i in range(7))
    H=[diagonal_matrix(QQ,[w[i] for w in weights],sparse=True) for i in range(7)]
    E=[]
    F=[]
    for i,a in enumerate(roots):
        e=matrix(QQ,56,56,sparse=True)
        f=matrix(QQ,56,56,sparse=True)
        for j,w in enumerate(weights):
            upper=tuple(w[k]+a[k] for k in range(7))
            lower=tuple(w[k]-a[k] for k in range(7))
            if upper in positions:
                e[positions[upper],j]=1
            if lower in positions:
                f[positions[lower],j]=1
        E.append(e)
        F.append(f)
    bracket=lambda x,y:x*y-y*x
    relation_checks=0
    for i in range(7):
        for j in range(7):
            assert bracket(H[i],E[j])==C[i,j]*E[j]
            assert bracket(H[i],F[j])==-C[i,j]*F[j]
            assert bracket(E[i],F[j])==(H[i] if i==j else matrix(QQ,56,56,sparse=True))
            relation_checks+=3
            if i!=j:
                for xs in (E,F):
                    value=xs[j]
                    for _ in range(1-int(C[i,j])):
                        value=bracket(xs[i],value)
                    assert not value
                    relation_checks+=1
    def root_vectors(generators,sign):
        result={tuple(sign*a for a in root):x for root,x in zip(roots,generators)}
        queue=list(result)
        while queue:
            beta=queue.pop()
            for i,x in enumerate(generators):
                gamma=tuple(beta[j]+sign*roots[i][j] for j in range(7))
                if gamma not in result:
                    value=bracket(x,result[beta])
                    if value:
                        result[gamma]=value
                        queue.append(gamma)
        assert len(result)==63
        return result
    positive=root_vectors(E,1)
    negative=root_vectors(F,-1)
    labels=tuple(sorted(positive))+tuple(sorted(negative))
    basis=H+[positive.get(a,negative.get(a)) for a in labels]
    labels=((0,)*7,)*7+labels
    for i in range(7):
        for x,a in zip(basis[7:],labels[7:]):
            assert bracket(H[i],x)==a[i]*x
    # Distinct root weights and nonzero root vectors, together with rank(H)=7,
    # prove linear independence of this 133-element Lie basis.
    assert matrix(QQ,weights).rank()==7
    return dict(weights=weights,basis=basis,labels=labels,relations=relation_checks)


def adjoint_matrices(oracle):
    basis=oracle['basis']
    labels=oracle['labels']
    lookup={a:i for i,a in enumerate(labels) if i>=7}
    pivot={i:next(iter(basis[i].dict())) for i in range(7,133)}
    W=matrix(QQ,oracle['weights'])
    rows=W.transpose().pivots()
    inverse=matrix(QQ,[W.row(i) for i in rows]).inverse()
    result=[]
    brackets=0
    for i,x in enumerate(basis):
        ad=matrix(QQ,133,133,sparse=True)
        for j,y in enumerate(basis):
            value=x*y-y*x
            label=tuple(a+b for a,b in zip(labels[i],labels[j]))
            if any(label):
                if label not in lookup:
                    assert not value
                else:
                    k=lookup[label]
                    coefficient=value[pivot[k]]/basis[k][pivot[k]]
                    assert value==coefficient*basis[k]
                    ad[k,j]=coefficient
            else:
                coefficients=inverse*vector(QQ,[value[k,k] for k in rows])
                assert value==sum((a*h for a,h in zip(coefficients,basis[:7])),
                                  matrix(QQ,56,56,sparse=True))
                for k,a in enumerate(coefficients):
                    ad[k,j]=a
            brackets+=1
        result.append(ad)
    return result,brackets
