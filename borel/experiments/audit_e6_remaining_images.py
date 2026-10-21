"""Independent proof-certificate replay for the remaining E6 comparisons.

No sage_borel or e6_remaining_images imports. SymPy recomputes Groebner
bases from literal 27 branching weights in a different monomial order.
An independent sparse reducer evaluates every Borel column. Newton sums
reconstruct the root-integrality conditions, and primitive-column checks
certify each proposed full source and image lattice. Smith multiplicities
come from a separate bounded local diagonalization, not Sage Smith reduction.
Sage is used only for polynomial and matrix arithmetic in this audit.
"""
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from time import perf_counter
import argparse,json,gzip
import sympy as sp
from gmpy2 import mpq
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector,prod,binomial

ROOT=Path(__file__).resolve().parents[1]


def powers(ws,d):
    if not ws:
        if d==0:yield ()
        return
    for a in range(d//ws[-1]+1):
        for e in powers(ws[:-1],d-a*ws[-1]):yield e+(a,)


def coefficient_inputs(kind):
    T=PolynomialRing(QQ,5,'u');u=T.gens()
    def moments(es,count,bound):
        out=[T(count)]
        for k in range(1,bound+1):
            out.append(sum((-1)**(j-1)*(k*es[j] if j==k else es[j]*out[k-j])
                           for j in range(1,min(k,len(es)-1)+1)))
        return out
    def shifted_elementary(es,count,shift,bound):
        raw=moments(es,count,bound)
        ps=[sum(binomial(k,j)*(-shift)**(k-j)*raw[j] for j in range(k+1)) for k in range(bound+1)]
        out=[T.one()]
        for k in range(1,bound+1):out.append(sum((-1)**(j-1)*out[k-j]*ps[j] for j in range(1,k+1))/k)
        return out
    if kind=='A5':
        E=[T.one()]+list(u)+[T.zero()]
        cs=shifted_elementary(E,6,u[0]/6,6)
        return T,cs[2:],[2,3,4,5,6]
    shift=(u[0]+u[3])/6
    a=shifted_elementary([T.one(),u[0],u[1],u[2]],3,shift,3)
    b=shifted_elementary([T.one(),u[3],u[4],T.zero()],3,shift,3)
    assert a[1]+b[1]==0
    return T,[a[1],a[2],a[3],b[2],b[3]],[1,2,3,2,3]


def independent_ring(kind):
    z=sp.symbols('z');cs=sp.symbols('c2 c3 c4 c5 c6')
    independent=json.loads((ROOT/'experiments/e6-independent-rational-base.json').read_text())
    assert independent['result']=='passed' and not independent['sage_imported']
    relations=list(map(sp.sympify,independent['relations']))
    BG=sp.groebner(relations,*cs,order='grlex',domain=sp.QQ)
    baseleads=[f.terms(order='grlex')[0][0] for f in BG.polys]
    base=[]
    for d in range(21):
        for e in powers((2,3,4,5,6),d):
            if not any(all(a>=b for a,b in zip(e,lm)) for lm in baseleads):base.append((d,e))
    assert len(base)==36
    reducers=[]
    if kind=='A5':
        variables=(z,)+cs;weights=(1,2,3,4,5,6)
        generators=list(variables);sourceweights=list(weights)
        reducers.append(((2,0,0,0,0,0),{(0,1,0,0,0,0):mpq(4)}))
        for f in BG.polys:
            lm,lc=f.terms(order='grlex')[0]
            tail={((0,)+e):mpq(-int((c/lc).p),int((c/lc).q)) for e,c in f.terms() if e!=lm}
            reducers.append(((0,)+lm,tail))
        basis_exponents=[(k,)+e for k in (0,1) for d,e in base]
    else:
        a1,a2,a3=sp.symbols('a1 a2 a3');c2,c3,c4,c5,c6=cs
        b2=c2-a2+a1*a1;b3=c3-a3-a1*c2+2*a1*a2-a1**3
        relative=[a1*b3+a2*b2-a1*a3-c4,a2*b3+a3*b2-c5,a3*b3-c6]
        RG=sp.groebner(relative,a1,a2,a3,domain=sp.QQ.frac_field(*cs),order='grlex')
        variables=(a1,a2,a3)+cs;weights=(1,2,3,2,3,4,5,6)
        generators=[a1,a2,a3,b2,b3,4*c2];sourceweights=[1,2,3,2,3,2]
        for f in RG.polys:
            lm,lc=f.terms(order='grlex')[0];lead=lm+(0,)*5
            p=sp.Poly(f.as_expr()/lc,*variables,domain=sp.QQ)
            tail={e:mpq(-int(c.p),int(c.q)) for e,c in p.terms() if e!=lead}
            reducers.append((lead,tail))
        assert {lm[:3] for lm,tail in reducers}==set(powers((1,1,1),4))
        for f in BG.polys:
            lm,lc=f.terms(order='grlex')[0]
            tail={((0,0,0)+e):mpq(-int((c/lc).p),int((c/lc).q)) for e,c in f.terms() if e!=lm}
            reducers.append(((0,0,0)+lm,tail))
        basis_exponents=[e+f for total in range(4) for e in powers((1,1,1),total) for d,f in base]
    zero=(0,)*len(variables)
    @lru_cache(None)
    def normal(e):
        for lm,tail in reducers:
            if all(a>=b for a,b in zip(e,lm)):
                rest=tuple(a-b for a,b in zip(e,lm));out={}
                for v,c in tail.items():
                    for w,k in normal(tuple(a+b for a,b in zip(rest,v))):out[w]=out.get(w,mpq(0))+c*k
                return tuple((w,c) for w,c in out.items() if c)
        return ((e,mpq(1)),)
    gp=[]
    for f in generators:
        gp.append({e:mpq(int(c.p),int(c.q)) for e,c in sp.Poly(f,*variables,domain=sp.QQ).terms()})
    @lru_cache(maxsize=20000)
    def image(e):
        if not any(e):return ((zero,mpq(1)),)
        i=next(i for i,a in enumerate(e) if a);lower=list(e);lower[i]-=1;out={}
        for a,c in image(tuple(lower)):
            for b,k in gp[i].items():
                for v,t in normal(tuple(x+y for x,y in zip(a,b))):out[v]=out.get(v,mpq(0))+c*k*t
        return tuple((v,c) for v,c in out.items() if c)
    return dict(normal=normal,image=image,weights=weights,sourceweights=sourceweights,
                base_groebner_count=len(BG.polys),variables=len(variables),basis_exponents=basis_exponents)


def local_smith(V):
    """Independent Python-list elimination modulo 81, with global pivots."""
    modulus=81;n=V.nrows()
    A=[[int(c.numerator())*pow(int(c.denominator()),-1,modulus)%modulus for c in row] for row in V]
    def valuation(c):
        if not c:return 4
        value=0
        while c%3==0:c//=3;value+=1
        return value
    answer=[]
    for r in range(n):
        v,i,j=min((valuation(A[i][j]),i,j) for i in range(r,n) for j in range(r,n))
        assert v<=3
        A[r],A[i]=A[i],A[r]
        for row in A:row[r],row[j]=row[j],row[r]
        p=3**v;u=pow(A[r][r]//p,-1,modulus)
        A[r]=[u*c%modulus for c in A[r]]
        for i in range(r+1,n):
            q=A[i][r]//p;A[i]=[(a-q*b)%modulus for a,b in zip(A[i],A[r])]
        for j in range(r+1,n):
            q=A[r][j]//p
            for i in range(n):A[i][j]=(A[i][j]-q*A[i][r])%modulus
        answer.append(v)
    return answer


def local(A):return A.denominator()%3!=0


def primitive_columns(A):
    assert local(A)
    return matrix(GF(3),A).rank()==A.nrows()


def audit_source(kind,row):
    T,gs,ws=coefficient_inputs(kind)
    es=[tuple(e) for e in row['exponents']]
    assert set(es)==set(powers(ws,row['degree']))
    fs=[prod(g**a for g,a in zip(gs,e)) for e in es]
    dictionaries=[f.dict() for f in fs]
    support=set().union(*(set(D) for D in dictionaries))
    coefficients=sorted(set(tuple(D.get(e,0) for D in dictionaries) for e in support))
    C=matrix(QQ,coefficients)
    n=len(es)
    B=(matrix(QQ,row['basis']) if 'basis' in row else
       matrix(QQ,n,n,{(i,j):QQ(c) for i,j,c in row['basis_sparse']}))
    assert C.nrows()==row['coefficient_rows'] and local(B)
    # Exact integral arithmetic checks every root coefficient. Column
    # primitivity modulo 3 proves this is the full local source lattice.
    den=C.denominator()*B.denominator()
    K=matrix(ZZ,C.denominator()*C)*matrix(ZZ,B.denominator()*B)
    p=ZZ(3)**den.valuation(3)
    assert all(c%p==0 for c in K.list())
    mod=matrix(GF(3),K.nrows(),K.ncols(),[GF(3)(c//p) for c in K.list()])
    assert mod.rank()==n
    return es,B


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--case',choices=('A5','A2A2A1'),default='A5')
    parser.add_argument('--max-degree',type=int);parser.add_argument('--resume',action='store_true');args=parser.parse_args()
    started=perf_counter();path=ROOT/('experiments/e6-remaining-'+args.case.lower()+'.json')
    data=json.loads(path.read_text());archive=None
    if 'certificate_archive' in data:
        archive=ROOT/data['certificate_archive']
        assert sha256(archive.read_bytes()).hexdigest()==data['certificate_archive_sha256']
        data=json.loads(gzip.decompress(archive.read_bytes()))
    print('AUDIT building independent SymPy ring',args.case,flush=True)
    ring=independent_ring(args.case)
    print('AUDIT rational basis ready',round(perf_counter()-started,2),flush=True)
    top=data['top_degree'];bound=min(top,args.max_degree if args.max_degree is not None else top)
    sources={};codehash=sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    result=dict(result='running',case=args.case,source_rows=[],degree_rows=[],
                runtime_imported=False,main_experiment_imported=False,auditor_sha256=codehash)
    target={d:[] for d in range(top+1)}
    for e in ring['basis_exponents']:
        d=sum(a*w for a,w in zip(e,ring['weights']))
        assert d<=top
        target[d].append(e)
    for d in target:target[d].sort()
    assert sum(map(len,target.values()))==data['rational_dimension']
    out=ROOT/('experiments/e6-remaining-'+args.case.lower()+'-audit.json')
    def rowhash(row):return sha256(json.dumps(row,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    first=0
    if args.resume and out.exists():
        result=json.loads(out.read_text());assert result['auditor_sha256']==codehash
        first=len(result['degree_rows'])
        for d in range(first):
            assert result['source_rows'][d]['input_sha256']==rowhash(data['sources'][d])
            assert result['degree_rows'][d]['input_sha256']==rowhash(data['degrees'][d])
            row=data['sources'][d];n=len(row['exponents'])
            B=matrix(QQ,row['basis']) if 'basis' in row else matrix(QQ,n,n,{(i,j):QQ(c) for i,j,c in row['basis_sparse']})
            sources[d]=([tuple(e) for e in row['exponents']],B)
    for d in range(first,bound+1):
        row=data['sources'][d];es,B=audit_source(args.case,row);sources[d]=(es,B)
        result['source_rows'].append(dict(degree=d,rank=len(es),complete=True,input_sha256=rowhash(row)))
        basis=target[d];positions={e:i for i,e in enumerate(basis)}
        def col(terms):
            v=vector(QQ,len(basis))
            for e,c in terms:
                assert e in positions
                v[positions[e]]=QQ((int(c.numerator),int(c.denominator)))
            return v
        ns=[];ms=[];cw=1 if args.case=='A5' else 2
        for k in range(d//cw+1):
            exps,BB=sources[d-cw*k]
            formal=[((k,)+e if args.case=='A5' else e+(k,)) for e in exps]
            A=(matrix(QQ,[col(ring['image'](tuple(e))) for e in formal]).transpose()
               if formal else matrix(QQ,len(basis),0))
            ns.append(A);ms.append(A*BB)
        N=ns[0];M=ms[0]
        for A in ns[1:]:N=N.augment(A)
        for A in ms[1:]:M=M.augment(A)
        prior=data['degrees'][d]
        change=matrix(QQ,[col(ring['normal'](tuple(e))) for e in prior['normal_exponents']]).transpose()
        assert change.is_square() and change.det()!=0
        NH=change*matrix(QQ,prior['weight_image_basis'])
        if 'root_image_basis' in prior:MH=change*matrix(QQ,prior['root_image_basis'])
        else:
            oldM=matrix(QQ,prior['root_image_columns']);den=oldM.denominator()
            MH=change*(matrix(ZZ,den*oldM).transpose().hermite_form(include_zero_rows=False).transpose()/den)
        assert primitive_columns(NH.inverse()*N)
        assert primitive_columns(MH.inverse()*M)
        V=NH.inverse()*MH
        assert local(V) and local(27*V.inverse())
        diagonal=local_smith(V)
        multiplicities=[diagonal.count(a) for a in (1,2,3)]
        assert all(a>=0 for a in multiplicities)
        expected=[prior['smith_exponents'].count(a) for a in (1,2,3)]
        assert multiplicities==expected
        result['degree_rows'].append(dict(degree=d,rank=len(basis),multiplicities=multiplicities,
                                          input_sha256=rowhash(data['degrees'][d])))
        result['seconds']=round(perf_counter()-started,3)
        out.write_text(json.dumps(result,indent=2)+'\n')
        print('AUDIT',args.case,d,multiplicities,result['seconds'],flush=True)
    result['result']='passed' if bound==top and data['result']=='passed' else 'partial'
    result['source_sha256']={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                             for p in (path,Path(__file__).resolve(),ROOT/'experiments/e6-independent-rational-base.json')}
    if archive is not None:result['raw_sha256']={archive.relative_to(ROOT).as_posix():sha256(archive.read_bytes()).hexdigest()}
    result['method']='SymPy grlex Groebner bases, independent sparse normal forms, Newton coefficient lattices, primitive lattice certificates, and independent Python local diagonalization modulo 81; Sage arithmetic is shared.'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['result'].upper(),args.case,result['seconds'],flush=True)


if __name__=='__main__':main()
