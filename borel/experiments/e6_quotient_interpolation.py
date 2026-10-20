"""Exact, finite-dimensional interpolation of higher E6 quotient coordinates.

The dominance bound proves completeness; a nonsingular evaluation matrix
over QQ determines the coefficients exactly. No probabilistic polynomial
identity claim or rational reconstruction is used.
"""
from pathlib import Path
from time import perf_counter
import json
import random
import sys
from sage.all import (QQ, ZZ, GF, CartanType, matrix, vector, prod, lcm)
from e6_geometry_support import candidates

C = matrix(QQ, CartanType('E6').cartan_matrix())
CINV = C.inverse()


def orbit(initial):
    found, queue = {tuple(initial)}, [tuple(initial)]
    while queue:
        w = queue.pop()
        for i in range(6):
            a = tuple(w[j]-w[i]*int(C[j,i]) for j in range(6))
            if a not in found:
                found.add(a)
                queue.append(a)
    return tuple(sorted(found))


WEIGHTS = orbit((1,0,0,0,0,0))
ROOTS = orbit(tuple(map(int,C.column(0))))
assert len(WEIGHTS)==27 and len(ROOTS)==72


def character_values(point, field):
    ts = tuple(map(field,point))
    powers=[{k:t**k for k in range(-2,3)} for t in ts]
    monomial=lambda w: prod(powers[i][int(a)] for i,a in enumerate(w))
    values=[monomial(w) for w in WEIGHTS]
    p1=sum(values)
    p2=sum(v*v for v in values)
    p3=sum(v**3 for v in values)
    dual=[1/v for v in values]
    d1=sum(dual)
    d2=sum(v*v for v in dual)
    return (p1, 6+sum(monomial(a) for a in ROOTS), (p1*p1-p2)/2,
            (p1**3-3*p1*p2+2*p3)/6, (d1*d1-d2)/2, d1)


def row_values(chars, exponents):
    bounds=[max(e[i] for e in exponents) for i in range(6)]
    powers=[[x**k for k in range(b+1)] for x,b in zip(chars,bounds)]
    return [prod(powers[i][k] for i,k in enumerate(e)) for e in exponents]


def target_values(point, name, degrees, mode='coordinates'):
    support=WEIGHTS if name=='minimal' else ROOTS
    index=6 if name=='minimal' else 24
    scale=18 if name=='minimal' else 24
    values=[prod(QQ(t)**a for t,a in zip(point,w)) for w in support]
    if mode=='determinant':
        return [matrix(QQ,6,6,lambda i,j:sum(w[i]*w[j]*v for w,v in zip(support,values))).det()]
    moment=vector(QQ,[sum(w[i]*v for w,v in zip(support,values)) for i in range(6)])
    h=CINV*moment/index
    eigenvalues=[sum(a*b for a,b in zip(w,h))*scale for w in WEIGHTS]
    return [sum(v**d for v in eigenvalues) for d in degrees]


def interpolate(name, degrees, *, seed=20261020, mode='coordinates'):
    start=perf_counter()
    highest=(1,0,0,0,0,0) if name=='minimal' else (0,1,0,0,0,0)
    exponents=candidates(highest,max(degrees))
    n=len(exponents)
    field=GF(1000003)
    rng=random.Random(seed+max(degrees)+(name=='adjoint'))
    points=[]
    rows=[]
    print(name,degrees,'interpolation size',n,flush=True)
    for _ in range(4):
        while len(points)<2*n:
            point=tuple(rng.choice((-5,-4,-3,-2,-1,1,2,3,4,5)) for _ in range(6))
            points.append(point)
            rows.append(row_values(character_values(point,field),exponents))
        A=matrix(field,rows)
        pivots=A.transpose().pivots()
        if len(pivots)==n:
            points=[points[i] for i in pivots]
            break
        print('rank',len(pivots),flush=True)
        nrows=len(points)+n
        while len(points)<nrows:
            point=tuple(rng.choice((-5,-4,-3,-2,-1,1,2,3,4,5)) for _ in range(6))
            points.append(point)
            rows.append(row_values(character_values(point,field),exponents))
    else:
        raise ArithmeticError('could not find a unisolvent point set')
    print('selected points',len(points),'seconds',perf_counter()-start,flush=True)
    exact_rows=[]
    targets=[]
    denominators=[]
    for point in points:
        row=row_values(character_values(point,QQ),exponents)
        rhs=target_values(point,name,degrees,mode)
        denominator=lcm(v.denominator() for v in row+rhs)
        denominators.append(int(denominator))
        exact_rows.append([ZZ(denominator*v) for v in row])
        targets.append([ZZ(denominator*v) for v in rhs])
    A=matrix(ZZ,exact_rows)
    B=matrix(ZZ,targets)
    print('exact matrix ready; max bits',max(v.nbits() for v in A.list()),'seconds',perf_counter()-start,flush=True)
    solution=A.solve_right(B)
    assert A*solution==B and solution.denominator()==1
    print('exact solve passed','seconds',perf_counter()-start,flush=True)
    results=[]
    for column,d in enumerate(degrees):
        terms=[(list(e),str(solution[i,column])) for i,e in enumerate(exponents) if solution[i,column]]
        allowed=set(candidates(highest,d))
        assert all(tuple(e) in allowed for e,c in terms)
        scale=(18 if name=='minimal' else 24)**d if mode=='coordinates' else 3*(6 if name=='minimal' else 24)**6
        results.append(dict(degree=d,terms=terms,scale=scale))
    return dict(representation=name,degrees=degrees,candidate_count=len(exponents),
                modulus=int(field.characteristic()),points=points,
                cleared_denominators=list(map(str,denominators)),coordinates=results,
                exact_equations_checked=True,modular_evaluation_rank=len(exponents),
                seconds=perf_counter()-start)


def main():
    name=sys.argv[1] if len(sys.argv)>1 else 'minimal'
    mode=sys.argv[2] if len(sys.argv)>2 else 'coordinates'
    blocks=(((2,5,8),(6,9,12)) if name=='minimal' else ((2,5,6,8,9,12),)) if mode=='coordinates' else ((6,),)
    cases=[]
    for block in blocks:
        case=interpolate(name,block,mode=mode)
        cases.append(case)
        suffix='quotient' if mode=='coordinates' else 'determinant'
        Path(__file__).with_name('e6-%s-%s-interpolation.json'%(name,suffix)).write_text(
            json.dumps(dict(result='passed',cases=cases),indent=2)+'\n')


if __name__=='__main__':
    main()
