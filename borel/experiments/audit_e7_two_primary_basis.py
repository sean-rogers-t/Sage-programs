"""Independent all-degree E7/P7 basis replay.

No runtime or main-experiment imports. Girard partitions replace Newton
recurrences; SymPy grlex and a sparse reducer replace Singular; Python
bit elimination replaces the primary mod-2 rank algorithm. Sage exact
polynomial/matrix arithmetic remains shared for the large products.
"""
from functools import lru_cache
from hashlib import sha256
from itertools import product
from math import factorial,comb
from pathlib import Path
from time import perf_counter
import json
import sympy as sp
from gmpy2 import mpq
from sage.all import QQ,ZZ,PolynomialRing,matrix,vector,prod,CartanMatrix

ROOT=Path(__file__).resolve().parents[1]
WS=(1,2,5,6,8,9,12)


def partitions(d,j=1):
    if j==7:
        if d==0:yield ()
        return
    for n in range(d//j+1):
        for e in partitions(d-j*n,j+1):yield (n,)+e


def weighted_monomials(ws,d):
    if not ws:
        if d==0:yield ()
        return
    for n in range(d//ws[-1]+1):
        for e in weighted_monomials(ws[:-1],d-n*ws[-1]):yield e+(n,)


def inputs():
    T=PolynomialRing(QQ,7,names=('E1','E2','E3','E4','E5','E6','u'))
    es=list(T.gens()[:6]);E=es[0];u=T.gen(6);ps=[T(6)]
    for k in range(1,13):
        f=T.zero()
        for e in partitions(k):
            length=sum(e);c=QQ((-1)**(k+length)*k*factorial(length-1))/prod(factorial(a) for a in e)
            f+=c*prod(g**a for g,a in zip(es,e))
        ps.append(f)
    result=[(3*u-E)/2]
    for d in WS[1:]:
        f=T.zero()
        for k in range(d+1):
            f+=comb(d,k)*6**k*ps[k]*((2*E-6*u)**(d-k)+(-4*E+6*u)**(d-k))
            pair=(sum(comb(k,j)*ps[j]*ps[k-j] for j in range(k+1))-2**k*ps[k])/2
            f+=comb(d,k)*(-6)**k*(2*E)**(d-k)*pair
        result.append(T(f))
    return T,result,[T(f(*(es+[E+2*u]))) for f in result]


def rational(x):return QQ((int(x.numerator),int(x.denominator)))
def bitrank(rows):
    pivots={}
    for row in rows:
        bits=0
        for j,c in enumerate(row):
            if int(c)%2:bits|=1<<j
        while bits:
            p=bits.bit_length()-1
            if p in pivots:bits^=pivots[p]
            else:pivots[p]=bits;break
    return len(pivots)


def coefficient_certificate(gs,es,B):
    fs=[prod(g**a for g,a in zip(gs,e)) for e in es];ds=[f.dict() for f in fs]
    support=set().union(*(set(f) for f in ds))
    C=matrix(QQ,sorted(set(tuple(f.get(e,0) for f in ds) for e in support)))
    den=C.denominator()*B.denominator();p=ZZ(2)**den.valuation(2)
    K=matrix(ZZ,C.denominator()*C)*matrix(ZZ,B.denominator()*B)
    assert all(c%p==0 for c in K.list())
    assert bitrank([[c//p for c in row] for row in K])==B.ncols()
    return C.nrows()


def target():
    h,s,t=sp.symbols('h s t')
    gb=sp.groebner([s*s-2*h*t,2*s*t-9*h**4*s*s+6*h**9*s-h**14,
                   t*t+10*h**3*s**3-9*h**8*s*s+2*h**13*s],h,s,t,order='grlex',domain=sp.QQ)
    leads=[f.terms(order='grlex')[0][0] for f in gb.polys]
    bounds=[min(e[i] for e in leads if e[i] and sum(e)==e[i]) for i in range(3)]
    standard=[e for e in product(*(range(n) for n in bounds))
              if not any(all(a>=b for a,b in zip(e,lm)) for lm in leads)]
    assert len(standard)==56 and max(a+5*b+9*c for a,b,c in standard)==27
    rules=[]
    for f in gb.polys:
        lm,lc=f.terms(order='grlex')[0]
        rules.append((lm,{e:mpq(-int((c/lc).p),int((c/lc).q)) for e,c in f.terms() if e!=lm}))
    @lru_cache(None)
    def normal(e):
        for lm,tail in rules:
            if all(a>=b for a,b in zip(e,lm)):
                rest=tuple(a-b for a,b in zip(e,lm));out={}
                for u,c in tail.items():
                    for v,k in normal(tuple(a+b for a,b in zip(rest,u))):out[v]=out.get(v,mpq(0))+c*k
                return tuple((v,c) for v,c in out.items() if c)
        return ((e,mpq(1)),)
    expressions=[h,-144*h*h,403200*h**5-933120*s,4720896*h**6-11197440*h*s,
                 -622420992*h**8+1463132160*h**3*s,
                 -1289723904*h**9-188116992*h**4*s+15237476352*t,
                 3066946633728*h**12-23932781199360*h**7*s+78538306682880*h**3*t]
    polynomials=[{e:mpq(int(c.p),int(c.q)) for e,c in sp.Poly(f,h,s,t).terms()} for f in expressions]
    @lru_cache(None)
    def image(es):
        if not any(es):return (((0,0,0),mpq(1)),)
        j=next(i for i,a in enumerate(es) if a);lower=list(es);lower[j]-=1;out={}
        for u,c in image(tuple(lower)):
            for v,k in polynomials[j].items():
                for e,a in normal(tuple(x+y for x,y in zip(u,v))):out[e]=out.get(e,mpq(0))+c*k*a
        return tuple((e,c) for e,c in out.items() if c)
    return normal,image,standard


def main():
    start=perf_counter();T,root,weight=inputs();normal,image,standard=target()
    primary_path=ROOT/'experiments/e7-two-primary-basis.json';primary=json.loads(primary_path.read_text());assert primary['result']=='passed'
    A=2*root[0];B=-root[1]/144;F5=root[2]/5760;F9=root[5]/387072
    assert all(c.denominator()==1 for f in (A,B,F5,F9) for c in f.coefficients())
    C=(35*A*B**2-F5)/81;D=(F9+1666*A*B**4+243*B**2*C)/19683
    for name,f in (('2h',A),('h2',B),('2s',C),('2t',D)):
        assert T(primary['source_lifts'][name])==f and all(c.denominator()%2 for c in f.coefficients())
    rpath=ROOT/'experiments/e7-character-form-results.json';wpath=ROOT/'experiments/e7-adjoint-comparison-results.json'
    rr=json.loads(rpath.read_text());ww=json.loads(wpath.read_text());assert rr['result']==ww['result']=='passed'
    # Independent actual-lattice decomposition: E6 simple roots plus h or 2h.
    C=matrix(QQ,CartanMatrix('E7'));h=vector(QQ,[0]*6+[1]);V=list(C.rows()[:6])
    XB=matrix(QQ,V+[h]);QB=matrix(QQ,V+[2*h])*C.inverse()
    assert XB.det()==3 and QB.det()==3
    for i in range(6):assert h[i]==0
    out=dict(result='running',degrees=[],root_complement_index=3,weight_complement_index=3,
             explicit_source_lifts_checked=4,refinement_exponent_bound=2,refinement_count=64,
             no_runtime_import=True,no_main_experiment_import=True)
    for d in range(28):
        target_es=sorted(e for e in standard if e[0]+5*e[1]+9*e[2]==d);pos={e:i for i,e in enumerate(target_es)}
        def col(terms):
            v=vector(QQ,len(pos))
            for e,c in terms:assert e in pos;v[pos[e]]=rational(c)
            return v
        labels=[(i,0) for i in range(28) if i==d]+[(i,b) for b in (5,9) for i in range(14) if i+b==d]
        columns=[]
        for i,b in labels:
            e=(i,int(b==5),int(b==9));columns.append((1 if b==0 else 2)*col(normal(e)))
        L=matrix(QQ,columns).transpose();assert L.is_square() and L.det()!=0
        es=[tuple(e) for e in rr['degrees'][d]['monomial_exponents']]
        assert set(es)==set(weighted_monomials(WS,d))
        assert es==[tuple(e) for e in ww['degrees'][d]['monomial_exponents']]
        A=L.inverse()*matrix(QQ,[col(image(e)) for e in es]).transpose()
        for kind,gs,row,key,factor in (('M',root,rr['degrees'][d],'source_basis_matrix',12),
                                       ('N',weight,ww['degrees'][d],'source_basis',36)):
            B=matrix(QQ,row[key])/(factor**d)
            coefficient_certificate(gs,es,B)
            Z=A*B
            scales=[2 if i%2 and kind=='M' else 1 for i,b in labels]
            W=matrix(QQ,[[Z[i,j]/scales[i] for j in range(Z.ncols())] for i in range(Z.nrows())])
            assert W.denominator()%2
            assert bitrank([[int(c.numerator())*pow(int(c.denominator()),-1,2)%2 for c in row] for row in W])==len(labels)
        multiplicity=sum(i%2 for i,b in labels)
        assert multiplicity==primary['degrees'][d]['multiplicity_2']
        out['degrees'].append(dict(degree=d,rank=len(labels),multiplicity_2=multiplicity,full_sources_checked=2))
        print('Independent E7',d,'Z/2',multiplicity,round(perf_counter()-start,2),flush=True)
    out.update(result='passed',seconds=round(perf_counter()-start,3),
               method='Closed Girard sums; SymPy grlex and separate sparse normal forms; Python bit elimination; Sage exact polynomial/matrix products shared.',
               source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                              for p in (Path(__file__).resolve(),primary_path,rpath,wpath)})
    (ROOT/'experiments/e7-two-primary-basis-audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASSED',out['seconds'],flush=True)


if __name__=='__main__':main()
