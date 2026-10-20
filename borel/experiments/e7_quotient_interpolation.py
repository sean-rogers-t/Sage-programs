"""CRT discovers coefficients; exact unisolvent equations certify them.

The dominance bound in E7_GEOMETRY.md proves that the finite candidate
space contains the true answer. Modular agreement alone never passes.
"""
from pathlib import Path
from time import perf_counter
import json
import random
import re
import sys
import numpy as np
from sage.all import ZZ, GF, matrix, vector, previous_prime
from sage.misc.cython import cython_import
from e7_geometry_support import (CINV, WEIGHTS, ROOTS, candidates, parent_plan,
    weights_912, integer_row, target_values)


def certificate_text(data):
    """Keep each point and coefficient on one readable line."""
    source=json.dumps(data,indent=2)
    source=re.sub(r'\[\n\s*((?:-?\d+,\n\s*)*-?\d+)\n\s*\]',
                  lambda m:'['+re.sub(r'\s+','',m.group(1))+']',source)
    source=re.sub(r'\[\n\s*(\[[^\n]*\]),\n\s*("[^"\n]*")\n\s*\]',r'[\1,\2]',source)
    assert json.loads(source)==data
    return source+'\n'


def interpolate(name, degrees, mode='coordinates'):
    start=perf_counter()
    highest=(0,0,0,0,0,0,1) if name=='minimal' else (1,0,0,0,0,0,0)
    exponents=candidates(highest,max(degrees))
    plan=parent_plan(exponents)
    n=len(exponents)
    extra=weights_912()
    assert max(abs(a) for w,m in extra for a in w)<=4
    kernel=cython_import(str(Path(__file__).with_name('e7_interpolation_kernels.pyx')),use_cache=True)
    support=WEIGHTS if name=='minimal' else ROOTS
    pairings=[[int((2 if name=='minimal' else 1)*(vector(mu)*CINV*vector(w)))
               for w in support] for mu in WEIGHTS]
    rng=random.Random(20261020+max(degrees)+(name=='adjoint'))
    points=[]
    seen=set()
    # Independent prime valuations distinguish the torus coordinates. Signs
    # in the first six slots give ample samples without enormous denominators.
    primes=(2,3,5,7,11,13,17)
    def add_points(count):
        while len(points)<count:
            p=tuple(q**rng.randrange(3)*(rng.choice((-1,1)) if i<6 else 1)
                    for i,q in enumerate(primes))
            if p not in seen:
                seen.add(p)
                points.append(p)
    modulus=ZZ.one()
    coefficients=None
    prime=8388593
    rank_certificate=None
    print(name,degrees,mode,'candidates',n,flush=True)
    for iteration in range(30):
        field=GF(prime)
        if coefficients is None:
            add_points(n+64)
        tick=perf_counter()
        chars,targets=kernel.base_rows(points,WEIGHTS,ROOTS,extra,pairings,degrees,prime)
        if mode=='determinant':
            targets=np.asarray([[int(field(target_values(p,name,degrees,mode)[0]))]
                                for p in points],dtype=np.int64)
        array=kernel.evaluation_rows(chars,plan[0],plan[1],prime)
        A=matrix(field,array)
        del array,chars
        if coefficients is None:
            pivots=A.transpose().pivots()
            if len(pivots)!=n:
                print('insufficient rank',len(pivots),'of',n,flush=True)
                add_points(len(points)+n//4+64)
                continue
            points=[points[i] for i in pivots]
            A=A.matrix_from_rows(pivots)
            targets=targets[list(pivots),:]
            rank_certificate=int(prime)
        B=matrix(field,targets)
        del targets
        try:
            solution=A.solve_right(B)
        except ValueError:
            prime=int(previous_prime(prime-1))
            continue
        assert A*solution==B
        del A,B
        if coefficients is None:
            coefficients=[[ZZ(solution[i,j]) for j in range(solution.ncols())] for i in range(n)]
            modulus=ZZ(prime)
        else:
            inverse=pow(int(modulus%prime),-1,int(prime))
            for i,row in enumerate(coefficients):
                for j,c in enumerate(row):
                    row[j]=c+modulus*(((int(solution[i,j])-int(c%prime))*inverse)%prime)
            modulus*=prime
        del solution
        signed=[[c if c<=modulus//2 else c-modulus for c in row] for row in coefficients]
        print('CRT',iteration+1,'bits',modulus.nbits(),'iteration seconds',perf_counter()-tick,flush=True)
        # Cheap exact rejection before the complete deterministic certificate.
        row,denominator=integer_row((2,-3,2,-3,2,-3,2),exponents,extra,plan)
        expected=target_values((2,-3,2,-3,2,-3,2),name,degrees,mode)
        if all(sum(a*c[j] for a,c in zip(row,signed))==denominator*b
               for j,b in enumerate(expected)):
            print('candidate passes probe; certifying all',n,'exact equations',flush=True)
            for k,p in enumerate(points):
                row,denominator=integer_row(p,exponents,extra,plan)
                rhs=target_values(p,name,degrees,mode)
                assert all(sum(a*c[j] for a,c in zip(row,signed))==denominator*b
                           for j,b in enumerate(rhs))
                if (k+1)%1000==0:
                    print('exact equations',k+1,'seconds',perf_counter()-start,flush=True)
            results=[]
            for j,d in enumerate(degrees):
                terms=[(list(e),str(c[j])) for e,c in zip(exponents,signed) if c[j]]
                allowed=set(candidates(highest,d))
                assert all(tuple(e) in allowed for e,c in terms)
                scale=(24 if name=='minimal' else 36)**d if mode=='coordinates' else 2*(12 if name=='minimal' else 36)**7
                results.append(dict(degree=d,terms=terms,scale=scale))
            return dict(representation=name,degrees=degrees,mode=mode,
                candidate_count=n,rank_prime=rank_certificate,modular_evaluation_rank=n,
                discovery_modulus=str(modulus),points=points,coordinates=results,
                exact_equations_checked=n,seconds=perf_counter()-start)
        prime=int(previous_prime(prime-1))
    raise ArithmeticError('CRT budget exhausted; no certified formula')


if __name__=='__main__':
    name=sys.argv[1]
    block=sys.argv[2] if len(sys.argv)>2 else 'low'
    mode='determinant' if block=='determinant' else 'coordinates'
    degrees=(7,) if mode=='determinant' else ((18,) if block=='high' else (2,6,8,10,12,14))
    result=interpolate(name,degrees,mode)
    path=Path(__file__).with_name('e7-%s-%s-interpolation.json'%(name,block))
    path.write_text(certificate_text(dict(result='passed',cases=[result])))
    print('saved',path.name,'seconds',result['seconds'],flush=True)
