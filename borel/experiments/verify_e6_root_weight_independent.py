"""Independent SymPy/stdlib audit. No Sage, runtime, or prior audit imports.

Rebuilds the D5 coordinates, generator BGG coefficients by sparse binomial
divided differences, all 17 relation/image lattices by determinantal
divisors, and all 32 D5 parabolic Poincare polynomials by signed permutations.
It checks the general transfer theorem's finite inputs, not its geometric
Leray--Hirsch justification, which is supplied in the proof note.
"""
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
from math import comb
from pathlib import Path
from time import perf_counter
import json

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
C = sp.Matrix([[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],
               [0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]])
Y = sp.Matrix([[3,0,0,0,0,0],[-1,-2,2,0,0,0],[-1,2,2,-2,0,0],
               [-1,0,0,2,-2,0],[-1,0,0,0,2,-2],[-1,0,0,0,0,2]])
ZERO = (0,)*6


def add(*terms):
    result = {}
    for coefficient, f in terms:
        for e,c in f.items():
            result[e] = result.get(e,0) + coefficient*c
    return {e:c for e,c in result.items() if c}


def mul(f,g):
    result = {}
    for a,c in f.items():
        for b,d in g.items():
            e=tuple(x+y for x,y in zip(a,b))
            result[e]=result.get(e,0)+c*d
    return {e:c for e,c in result.items() if c}


@lru_cache(None)
def root_power(i,k):
    if k == 0:
        return {ZERO:1}
    alpha={tuple(int(j==r) for j in range(6)):int(C[i,r]) for r in range(6) if C[i,r]}
    return mul(root_power(i,k-1),alpha)


@lru_cache(None)
def dd_monomial(e,i):
    out={}
    for k in range(1,e[i]+1):
        lower=list(e); lower[i]-=k
        for a,c in root_power(i,k-1).items():
            b=tuple(x+y for x,y in zip(lower,a))
            out[b]=out.get(b,0)+(-1)**(k+1)*comb(e[i],k)*c
    return {a:c for a,c in out.items() if c}


def coefficient(f,word):
    for i in reversed(word):
        out={}
        for e,c in f.items():
            for a,b in dd_monomial(e,i-1).items():
                out[a]=out.get(a,0)+c*b
        f={e:c for e,c in out.items() if c}
    assert all(e==ZERO for e in f)
    return f.get(ZERO,0)


def v3(x):
    x=sp.Rational(x)
    if x == 0:
        return 10**9
    a,b=abs(int(x.p)),int(x.q)
    result=0
    while a%3==0: a//=3;result+=1
    while b%3==0: b//=3;result-=1
    return result


def islocal(A):
    return all(sp.denom(x)%3 for x in A)


def monomials(weights,total):
    if total < 0:
        return []
    result=[()]
    used=[0]
    for w in weights:
        following=[]; totals=[]
        for es,u in zip(result,used):
            for a in range((total-u)//w+1):
                following.append(es+(a,)); totals.append(u+a*w)
        result,used=following,totals
    return [es for es,u in zip(result,used) if u==total]


def d5_polynomials():
    roots=[]
    for i in range(5):
        for j in range(i+1,5):
            for sign in (-1,1):
                v=[0]*5;v[i]=1;v[j]=sign;roots.append(tuple(v))
    simple=[tuple(int(j==i)-int(j==i+1) for j in range(5)) for i in range(4)]
    simple.append((0,0,0,1,1))
    counts=[[0]*21 for _ in range(32)]
    def negative(v):
        return next(c for c in v if c) < 0
    for perm in permutations(range(5)):
        for signs in product((-1,1),repeat=5):
            if signs.count(-1)%2:
                continue
            def act(v):
                out=[0]*5
                for i,c in enumerate(v):out[perm[i]]=signs[i]*c
                return out
            length=sum(negative(act(a)) for a in roots)
            desc=sum(1<<i for i,a in enumerate(simple) if negative(act(a)))
            for mask in range(32):
                if not mask&desc:counts[mask][length]+=1
    for c in counts:
        while c[-1]==0:c.pop()
    return counts


def main():
    start=perf_counter()
    path=ROOT/'experiments/e6-root-weight-comparison.json'
    data=json.loads(path.read_text())
    assert data['result']=='passed'
    assert C.det()==3 and Y.det()==-192
    q=Y*C.inv()
    assert all(x.q==1 for x in q)
    # The y1, y2,...,y6 basis of Q has index 64. Replacing y1=3h
    # by h gives an index-64 basis of X. All indices are 3-local units.
    assert abs(q.det())==64
    assert abs(sp.Matrix([[1,0,0,0,0,0]]+Y.tolist()[1:]).det())==64
    def reflect_weight(v,i):
        return tuple(v[j]-v[i]*int(C[i,j]) for j in range(6))
    assert all(reflect_weight(tuple(Y.row(0)),i)==tuple(Y.row(0)) for i in range(1,6))
    for i in range(1,6):
        for row in Y.tolist()[1:]:
            value=reflect_weight(row,i)
            assert value in [tuple(r) for r in Y.tolist()[1:]]+[
                tuple(-c for c in r) for r in Y.tolist()[1:]]
    distances={(1,0,0,0,0,0):0}
    frontier=list(distances)
    for v in frontier:
        for i in range(6):
            u=reflect_weight(v,i)
            if u not in distances:
                distances[u]=distances[v]+1;frontier.append(u)
    assert len(distances)==27
    seen=set()
    for d,words in data['schubert_labels'].items():
        for word in words:
            value=(1,0,0,0,0,0)
            for i in reversed(word):value=reflect_weight(value,i-1)
            assert distances[value]==len(word)==int(d) and value not in seen
            seen.add(value)
    assert seen==set(distances)
    linear=[{tuple(int(j==i) for j in range(6)):int(c) for i,c in enumerate(row) if c}
            for row in Y.tolist()]
    xs=[mul(f,f) for f in linear[1:]]
    def elementary(k):
        result={}
        for subset in combinations(xs,k):
            f={ZERO:1}
            for g in subset:f=mul(f,g)
            result=add((1,result),(1,f))
        return result
    euler={ZERO:1}
    for f in linear[1:]:euler=mul(euler,f)
    inputs=[linear[0],elementary(1),elementary(2),euler,elementary(3),elementary(4)]
    names=('y1','e1','e2','euler','e3','e4')
    bgg=0
    for name,f in zip(names,inputs):
        d=sum(next(iter(f)))
        got={tuple(term['word']):sp.Rational(term['coefficient']) for term in data['generator_schubert'][name]}
        for word in data['schubert_labels'][str(d)]:
            assert coefficient(f,word)==got.get(tuple(word),0),(name,word)
            bgg+=1
    print('Independent sparse BGG coefficients passed',bgg,flush=True)
    h,s=sp.symbols('h s')
    gs=[3*h,-3*h*h,42*h**4-96*s,-21*h**5+48*h*s,
        282*h**6-672*h*h*s,-1371*h**8+3552*h**4*s-768*s*s]
    r9=2*h**9-6*h**5*s+3*h*s*s
    r12=h**12-6*h**4*s*s+s**3
    matrix_checks=0
    for row in data['degrees']:
        d=row['degree']
        exps=monomials((1,4),d)
        chosen=[e for e in exps if e[0]<9 and e[1]<3]
        # Follow the report's basis order only after verifying its content.
        reported=[tuple(e) for e in row['monomial_basis']]
        assert set(chosen)==set(reported)
        chosen=reported
        basis=[h**i*s**j for i,j in chosen]
        relations=[sp.expand(r*h**i*s**j) for r,k in ((r9,9),(r12,12))
                   for i,j in monomials((1,4),d-k)]
        def vec(f):
            p=sp.Poly(f,h,s)
            return sp.Matrix([p.coeff_monomial(h**i*s**j) for i,j in exps])
        A=sp.Matrix.hstack(*(vec(f) for f in basis+relations))
        assert A.rows==A.cols and v3(A.det())==0
        coordinate=A.inv()[:len(basis),:]
        scaling=sp.diag(*[1 if j==0 else 3 for i,j in chosen])
        def image(f):return scaling.inv()*coordinate*vec(f)
        source=monomials((1,2,4,5,6,8),d)
        assert set(source)==set(map(tuple,row['root_monomials']))
        source=[tuple(e) for e in row['root_monomials']]
        cols=[image(sp.expand(sp.prod(g**a for g,a in zip(gs,es)))) for es in source]
        T=sp.Matrix.hstack(*cols)
        assert T==sp.Matrix(row['image_columns_in_weight_basis']) and islocal(T)
        H=sp.Matrix(row['root_basis_in_weight_basis'])
        assert H==T[:,row['selected_root_monomials']]
        assert islocal(H.inv()*T)
        assert T*sp.Matrix(row['root_basis_witness'])==H
        # Recompute determinantal valuations from all minors of H, without
        # invoking Sage's Smith implementation or accepting its exponents.
        previous=0; factors=[]
        for k in range(1,H.rows+1):
            value=min(v3(H.extract(rs,cs).det()) for rs in combinations(range(H.rows),k)
                      for cs in combinations(range(H.cols),k))
            factors.append(value-previous);previous=value
        assert factors==row['smith_exponents']
        matrix_checks+=1
    counts=d5_polynomials()
    for index,row in enumerate(data['d5_refinements']):
        mask=sum(1<<(i-1) for i in row['d5_retained'])
        assert counts[mask]==row['fiber_betti']
        for a in (1,2,3):
            pol=[0]*(17+len(counts[mask])-1)
            for base in data['degrees']:
                for k,c in enumerate(counts[mask]):
                    pol[k+base['degree']]+=c*base['smith_exponents'].count(a)
            assert pol==row['multiplicities'][str(a)]
    assert sum(counts[0])==1920
    # Independently enumerate each Levi on the faithful 27-weight orbit.
    weights=sorted(distances); positions={v:i for i,v in enumerate(weights)}
    reflections=[tuple(positions[reflect_weight(v,i)] for v in weights) for i in range(6)]
    def modrank(A):
        A=[[int(c)%3 for c in row] for row in A];r=0
        for j in range(len(A[0]) if A else 0):
            k=next((k for k in range(r,len(A)) if A[k][j]),None)
            if k is None:continue
            A[r],A[k]=A[k],A[r];a=pow(A[r][j],-1,3)
            A[r]=[a*c%3 for c in A[r]]
            for k in range(r+1,len(A)):
                a=A[k][j];A[k]=[(b-a*c)%3 for b,c in zip(A[k],A[r])]
            r+=1
        return r
    q=sp.symbols('q')
    D5=sum(c*q**i for i,c in enumerate(counts[0]))
    full=[sp.expand(D5*sum(q**r['degree']*r['smith_exponents'].count(a) for r in data['degrees'])) for a in (1,2,3)]
    good_count=0; exceptions=[]
    for row in data['proper_parabolics']:
        I=[i-1 for i in row['e6_retained']]
        identity=tuple(range(27));elts={identity:0};queue=[identity]
        for w in queue:
            for i in I:
                t=tuple(reflections[i][j] for j in w)
                if t not in elts:elts[t]=elts[w]+1;queue.append(t)
        cs=[0]*(max(elts.values())+1)
        for length in elts.values():cs[length]+=1
        assert cs==row['levi_length_polynomial']
        good=modrank([[C[j,i] for i in I] for j in range(6)])==len(I)
        assert good==row['common_fiber_basis']
        divisor=sum(c*q**i for i,c in enumerate(cs))
        for a,f in zip((1,2,3),full):
            quo,rem=sp.div(f,divisor,q)
            assert (rem==0)==good
            if good:
                coeff=[int(sp.expand(quo).coeff(q,i)) for i in range(sp.degree(quo,q)+1)]
                assert all(c>=0 for c in coeff) and coeff==row['multiplicities'][str(a)]
        if good:good_count+=1
        else:exceptions.append(tuple(row['e6_retained']))
    assert good_count==60
    assert set(exceptions)=={(1,3,5,6),(1,2,3,5,6),(1,3,4,5,6)}
    result=dict(result='passed', sage_imported=False, runtime_imported=False,
                generator_bgg_coefficients=bgg, independently_checked_minuscule_labels=27,
                independently_rebuilt_image_lattices=matrix_checks,
                independently_enumerated_d5_refinements=32,
                independently_enumerated_e6_levi_polynomials=63,
                sharp_deconvolution_cases=good_count, exceptional_retained_sets=exceptions,
                scope='Independent arithmetic, source-coordinate and BGG checks; geometric source/transfer theorems are separately proved, not inferred from the finite checks.',
                seconds=round(perf_counter()-start,3),
                source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for p in (path,Path(__file__).resolve())})
    (ROOT/'experiments/e6-root-weight-independent.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASSED',result,flush=True)


if __name__=='__main__':main()
