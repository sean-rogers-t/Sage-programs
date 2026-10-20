"""Dominance bounds, orbits, and character evaluations for E7 geometry."""
from sage.all import QQ, ZZ, matrix, vector, CartanType, WeylCharacterRing, prod

C=matrix(QQ,CartanType('E7').cartan_matrix())
CINV=C.inverse()
MARKS=(2,2,3,4,3,2,1)
DEGREES=(2,6,8,10,12,14,18)


def orbit(initial):
    seen,queue={tuple(initial)},[tuple(initial)]
    while queue:
        w=queue.pop()
        for i in range(7):
            v=tuple(w[j]-int(C[j,i])*w[i] for j in range(7))
            if v not in seen:
                seen.add(v)
                queue.append(v)
    return tuple(sorted(seen))


WEIGHTS=orbit((0,0,0,0,0,0,1))
ROOTS=orbit(tuple(map(int,C.column(0))))
assert len(WEIGHTS)==56 and len(ROOTS)==126


def weights_912():
    R=WeylCharacterRing('E7',style='coroots')
    coroots=R.space().simple_coroots()
    return tuple(sorted((tuple(int(w.inner_product(coroots[i+1])) for i in range(7)),int(m))
                         for w,m in R((0,1,0,0,0,0,0)).weight_multiplicities().items()))


def candidates(highest,degree):
    target=CINV*vector(highest)*degree
    columns=CINV.columns()
    answer=[]
    def walk(i,left,prefix):
        if i==7:
            if all(a.denominator()==1 for a in left):
                answer.append(tuple(prefix))
            return
        bound=min(a//b for a,b in zip(left,columns[i]))
        for k in range(int(bound)+1):
            walk(i+1,left-k*columns[i],prefix+[k])
    walk(0,target,[])
    return tuple(sorted(answer,key=lambda e:(sum(e),e)))


def parent_plan(exponents):
    positions={e:i for i,e in enumerate(exponents)}
    parents=[]
    factors=[]
    costs=[]
    for e in exponents:
        costs.append(sum(a*b for a,b in zip(e,MARKS)))
        if not any(e):
            parents.append(-1)
            factors.append((-1,-1))
            continue
        if sum(e)==1 and tuple(0 for _ in e) not in positions:
            parents.append(-1)
            factors.append((e.index(1),-1))
            continue
        even=[i for i in (0,2,3,5) if e[i]]
        if even:
            chosen=(even[0],)
        else:
            odd=[i for i in (1,4,6) for _ in range(e[i])]
            chosen=tuple(odd[:2])
        p=list(e)
        for i in chosen:
            p[i]-=1
        parents.append(positions[tuple(p)])
        factors.append((chosen[0],chosen[1] if len(chosen)>1 else -1))
    return parents,factors,costs


def character_values(point,field,extra):
    ts=tuple(map(field,point))
    def value(w):
        return prod(t**a for t,a in zip(ts,w))
    values=[value(w) for w in WEIGHTS]
    p=[sum(v**k for v in values) for k in (1,2,3,4)]
    a1=7+sum(value(w) for w in ROOTS)
    a2=7+sum(value(tuple(2*x for x in w)) for w in ROOTS)
    e2=(p[0]**2-p[1])/2
    e3=(p[0]**3-3*p[0]*p[1]+2*p[2])/6
    e4=(p[0]**4-6*p[0]**2*p[1]+3*p[1]**2+8*p[0]*p[2]-6*p[3])/24
    return (a1,sum(m*value(w) for w,m in extra),(a1*a1-a2)/2-a1,
            e4-e2,e3-p[0],e2-1,p[0])


def integer_row(point,exponents,extra,plan=None):
    chars=character_values(point,QQ,extra)
    B=ZZ(prod(abs(t) for t in point))
    numerators=[ZZ(c*B**a) for c,a in zip(chars,MARKS)]
    parents,factors,costs=parent_plan(exponents) if plan is None else plan
    K=max(costs)
    base_powers=[B**k for k in range(K+1)]
    monomials=[]
    for parent,(i,j) in zip(parents,factors):
        if parent<0:
            monomials.append(ZZ.one() if i<0 else numerators[i])
        else:
            v=monomials[parent]*numerators[i]
            if j>=0:
                v*=numerators[j]
            monomials.append(v)
    return [m*base_powers[K-c] for m,c in zip(monomials,costs)],B**K


def target_values(point,name,degrees,mode='coordinates'):
    support=WEIGHTS if name=='minimal' else ROOTS
    index=12 if name=='minimal' else 36
    values=[prod(QQ(t)**a for t,a in zip(point,w)) for w in support]
    if mode=='determinant':
        return [matrix(QQ,7,7,lambda i,j:sum(w[i]*w[j]*v for w,v in zip(support,values))).det()]
    moment=vector(QQ,[sum(w[i]*v for w,v in zip(support,values)) for i in range(7)])
    h=CINV*moment/index
    scale=24 if name=='minimal' else 36
    eigenvalues=[sum(a*b for a,b in zip(w,h))*scale for w in WEIGHTS]
    return [sum(v**d for v in eigenvalues) for d in degrees]
