"""Exact root/weight image comparison for the three bad E6 Levis.

Coordinates use the full-rank A5+A1 subsystem. Rational Borel ideals
come from six independently checked E6 trace invariants. Root-integral
invariants are certified by compressed coefficient lattices, retaining
the complete rational kernel before forming each image.
"""
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from time import perf_counter
import argparse
import json
import gzip

from sage.all import (QQ, ZZ, GF, PolynomialRing, matrix, vector, prod,
                      binomial, CartanMatrix, diagonal_matrix, TermOrder, inverse_mod,
                      pari)

ROOT=Path(__file__).resolve().parents[1]
C=matrix(ZZ,CartanMatrix('E6'))


def powers(weights,total):
    if total<0:return
    if not weights:
        if total==0:yield ()
        return
    for e in range(total//weights[0]+1):
        for tail in powers(weights[1:],total-e*weights[0]):
            yield (e,)+tail


def rows(A):return [[str(c) for c in r] for r in A]
def local(A):return all(c.denominator()%3 for c in A.list())


def a5_coordinates(blocks=False):
    alpha=[vector(QQ,C.row(i)) for i in range(6)]
    if not blocks:
        chain=[alpha[i-1] for i in (1,3,4,5,6)]
        z=vector(QQ,[0,1,0,0,0,0])
    else:
        gamma=vector(QQ,[0,0,-1,1,-1,0])
        chain=[alpha[0],alpha[2],gamma,alpha[4],alpha[5]]
        z=alpha[1]
    x6=-sum(((i+1)*a for i,a in enumerate(chain)),vector(QQ,[0]*6))/6
    xs=[x6+sum(chain[i:],vector(QQ,[0]*6)) for i in range(5)]+[x6]
    assert not sum(xs,vector(QQ,[0]*6))
    rootbasis=matrix(QQ,[z]+[x-xs[-1] for x in xs[:-1]])
    weightbasis=matrix(QQ,[z]+xs[:-1])
    assert local(rootbasis*C.inverse()) and local(weightbasis)
    assert (rootbasis*C.inverse()).det().valuation(3)==0
    assert weightbasis.det().valuation(3)==0
    return xs,z


def trace_invariants(P,elementary,zsquare):
    p=[P(6),P.zero()]
    for k in range(2,13):
        answer=P.zero()
        for j in range(1,min(k,6)+1):
            if j<k:answer+=(-1)**(j-1)*elementary[j]*p[k-j]
            elif j==k:answer+=(-1)**(j-1)*k*elementary[j]
        p.append(answer)
    result=[]
    for k in (2,5,6,8,9,12):
        pair=(sum(binomial(k,j)*p[j]*p[k-j] for j in range(k+1))-2**k*p[k])/2
        doublet=2*sum(binomial(k,j)*(-1)**(k-j)*p[k-j]*zsquare**(j//2)/2**j
                      for j in range(0,k+1,2))
        result.append(P(pair+doublet))
    return result


def context(kind):
    if kind=='A5':
        P=PolynomialRing(QQ,6,names=('z','c2','c3','c4','c5','c6'),order='degrevlex')
        z,*cs=P.gens();weights=(1,2,3,4,5,6)
        elementary=[P.one(),P.zero()]+cs
        relations=trace_invariants(P,elementary,z*z)
        T=PolynomialRing(QQ,5,names=('E1','E2','E3','E4','E5'))
        E=[T.one()]+list(T.gens())+[T.zero()];mu=E[1]/6
        inputs=[T(sum(binomial(6-j,k-j)*(-mu)**(k-j)*E[j] for j in range(k+1)))
                for k in range(2,7)]
        baseweights=(2,3,4,5,6);compressionweights=(1,2,3,4,5)
        I=(1,3,4,5,6);top=21;dimension=72
    else:
        P=PolynomialRing(QQ,6,names=('h','a2','a3','b2','b3','t'),order='degrevlex')
        h,a2,a3,b2,b3,t=P.gens();weights=(1,2,3,2,3,2)
        A=[P.one(),h,a2,a3];B=[P.one(),-h,b2,b3]
        elementary=[sum(A[i]*B[k-i] for i in range(4) if 0<=k-i<4) for k in range(7)]
        relations=trace_invariants(P,elementary,t)
        T=PolynomialRing(QQ,5,names=('U1','U2','U3','V1','V2'))
        u1,u2,u3,v1,v2=T.gens();mu=(u1+v1)/6
        inputs=[(u1-v1)/2,u2-2*mu*u1+3*mu**2,
                u3-mu*u2+mu**2*u1-mu**3,
                v2-2*mu*v1+3*mu**2,-mu*v2+mu**2*v1-mu**3]
        inputs=list(map(T,inputs));baseweights=(1,2,3,2,3);compressionweights=(1,2,3,1,2)
        I=(1,2,3,5,6);top=29;dimension=720
    source_inputs=P.gens();normal_weights=weights
    print(kind,'building rational Borel basis',flush=True)
    if kind=='A5':
        G=P.ideal(relations).groebner_basis()
        assert P.ideal(G).vector_space_dimension()==dimension
        def poly(es):return prod(x**a for x,a in zip(P.gens(),es))
        normal={d:[e for e in powers(normal_weights,d) if poly(e).reduce(G)==poly(e)] for d in range(top+1)}
        assert not any(poly(e).reduce(G) for e in powers(normal_weights,top+1))
    else:
        # Full-rank A5+A1 invariant quotient has dimension 36. Adjoin the
        # universal Gr(3,6) Chern presentation with its 20-element basis.
        S=PolynomialRing(QQ,5,names=('c2','c3','c4','c5','c6'))
        cs=S.gens();br=trace_invariants(S,[S.one(),S.zero()]+list(cs),4*cs[0])
        assert br[0]==0
        BG=S.ideal(br[1:]).groebner_basis()
        assert S.ideal(BG).vector_space_dimension()==36
        base=[]
        for d in range(21):
            for e in powers((2,3,4,5,6),d):
                f=prod(x**a for x,a in zip(cs,e))
                if f.reduce(BG)==f:base.append((d,e))
        assert len(base)==36
        P8=PolynomialRing(QQ,8,names=('a1','a2','a3','c2','c3','c4','c5','c6'),
                          order=TermOrder('degrevlex',3)+TermOrder('degrevlex',5))
        a1,a2,a3,c2,c3,c4,c5,c6=P8.gens()
        b2=c2-a2+a1*a1;b3=c3-a3-a1*c2+2*a1*a2-a1**3
        relative=[a1*b3+a2*b2-a1*a3-c4,a2*b3+a3*b2-c5,a3*b3-c6]
        RG=P8.ideal(relative).groebner_basis()
        assert {tuple(f.lm().exponents()[0]) for f in RG}=={
            e+(0,)*5 for e in powers((1,1,1),4)}
        embed=S.hom(P8.gens()[3:],P8)
        G=list(RG)+[embed(f) for f in BG]
        source_inputs=(a1,a2,a3,b2,b3,4*c2)
        original=P.hom(source_inputs,P8)
        assert all(original(f).reduce(G)==0 for f in relations)
        P=P8;normal_weights=(1,2,3,2,3,4,5,6)
        normal={d:[] for d in range(top+1)}
        for total in range(4):
            for e in powers((1,1,1),total):
                degree=sum(a*w for a,w in zip(e,(1,2,3)))
                for birth,f in base:normal[degree+birth].append(e+f)
        for d in normal:normal[d].sort()
    assert sum(map(len,normal.values()))==dimension
    return dict(P=P,T=T,weights=weights,baseweights=baseweights,
                compressionweights=compressionweights,inputs=inputs,relations=relations,
                groebner=G,normal=normal,retained=I,top=top,dimension=dimension,
                image_inputs=source_inputs,normal_weights=normal_weights)


def root_source(ctx,d):
    started=perf_counter()
    es=list(powers(ctx['baseweights'],d));T=ctx['T'];inputs=ctx['inputs']
    fs=[T(prod(g**a for g,a in zip(inputs,e))) for e in es]
    dictionaries=[f.dict() for f in fs]
    support=sorted(set().union(*(set(D) for D in dictionaries)))
    # Unique coefficient rows preserve precisely the same integrality conditions.
    coefficients=sorted(set(tuple(D.get(e,0) for D in dictionaries) for e in support))
    M=matrix(QQ,coefficients)
    # Only the 3-local lattice is needed. Reduce the dual row lattice
    # modulo 3^k, adding 3^k I to impose weight-integral coefficients.
    den=ZZ(3)**max(0,int(M.denominator().valuation(3)))
    if den==1:
        B=diagonal_matrix(QQ,[1]*len(es))
    else:
        reduced=matrix(ZZ,M.nrows(),M.ncols(),
            [ZZ(c.numerator()*inverse_mod(c.denominator(),den))%den for c in (den*M).list()])
        print('source',d,'local pivots',M.dimensions(),round(perf_counter()-started,2),flush=True)
        # PARI's modular HNF computes exactly the lattice generated by
        # these columns and den*I, without unbounded intermediate entries.
        H=matrix(ZZ,pari(reduced.transpose()).mathnfmodid(den)).transpose()
        B=den*H.change_ring(QQ).inverse()
    assert local(B)
    print('source',d,'basis ready',round(perf_counter()-started,2),flush=True)
    denominator=M.denominator()*B.denominator()
    integer=matrix(ZZ,M.denominator()*M)*matrix(ZZ,B.denominator()*B)
    ppower=ZZ(3)**denominator.valuation(3)
    assert all(c%ppower==0 for c in integer.list())
    unit=GF(3)(denominator//ppower)**(-1)
    mod=matrix(GF(3),integer.nrows(),integer.ncols(),
               [unit*GF(3)(c//ppower) for c in integer.list()])
    # Completeness over Z_(3): reduction of these coefficient columns has
    # full rank. This is a proof certificate for the whole local source.
    pivots=list(mod.transpose().pivots())
    assert len(pivots)==len(es)
    minor=mod.matrix_from_rows(pivots)
    assert minor.det()!=0
    print('source',d,'primitive certificate',round(perf_counter()-started,2),flush=True)
    return dict(exponents=es,basis=B,coefficient_matrix=M,
                pivot_rows=pivots,minor_mod3=(rows(minor) if len(es)<=64 else None))


def verify_trace_branching(kind):
    xs,z=a5_coordinates(kind!='A5')
    highest=vector(QQ,[1,0,0,0,0,0]);found={tuple(highest)};queue=[highest]
    for w in queue:
        for i in range(6):
            v=w-w[i]*vector(QQ,C.row(i))
            if tuple(v) not in found:found.add(tuple(v));queue.append(v)
    proposed={tuple(xs[i]+xs[j]) for i,j in combinations(range(6),2)}
    proposed|={tuple(-x+sign*z/2) for x in xs for sign in (-1,1)}
    if proposed!=found:
        assert {tuple(-c for c in w) for w in proposed}==found
    return dict(x_coordinates=[list(map(str,x)) for x in xs],z_coordinates=list(map(str,z)))


def bound_for_resume(ctx,args):
    return min(ctx['top'],args.max_degree if args.max_degree is not None else ctx['top'])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--case',choices=('A5','A2A2A1'),default='A5')
    parser.add_argument('--max-degree',type=int);parser.add_argument('--resume',action='store_true');args=parser.parse_args()
    started=perf_counter();ctx=context(args.case)
    if args.case!='A5':pari.allocatemem(512*1024**2,silent=True)
    P=ctx['P'];G=ctx['groebner'];weights=ctx['weights'];xs=P.gens()
    output=ROOT/('experiments/e6-remaining-'+args.case.lower()+'.json')
    report=dict(result='running',case=args.case,retained=list(ctx['retained']),top_degree=ctx['top'],
                rational_dimension=ctx['dimension'],grading=list(weights),normal_grading=list(ctx['normal_weights']),
                relations=[str(f) for f in ctx['relations']],
                branching=verify_trace_branching(args.case),sources=[],degrees=[])
    source={}
    first=0
    if args.resume and output.exists():
        report=json.loads(output.read_text())
        if 'certificate_archive' in report:
            archive=ROOT/report['certificate_archive']
            assert sha256(archive.read_bytes()).hexdigest()==report['certificate_archive_sha256']
            report=json.loads(gzip.decompress(archive.read_bytes()))
        first=len(report['degrees'])
        if first>bound_for_resume(ctx,args):
            print('Already complete through',first-1,flush=True);return
        for r in report['sources']:
            n=len(r['exponents'])
            basis=matrix(QQ,r['basis']) if 'basis' in r else matrix(QQ,n,n,{(i,j):QQ(c) for i,j,c in r['basis_sparse']})
            source[r['degree']]=dict(exponents=[tuple(e) for e in r['exponents']],basis=basis)
    bound=bound_for_resume(ctx,args)
    for d in range(first,bound+1):
        sd=root_source(ctx,d);source[d]=sd
        n=sd['basis'].nrows()
        encoded=({'basis':rows(sd['basis'])} if n<=128 else
                 {'basis_sparse':[[int(i),int(j),str(c)] for (i,j),c in sd['basis'].dict().items()]})
        report['sources'].append(dict(degree=d,exponents=[list(e) for e in sd['exponents']],
           **encoded,coefficient_rows=sd['coefficient_matrix'].nrows(),
           pivot_rows=sd['pivot_rows'],minor_mod3=sd['minor_mod3']))
        target=ctx['normal'][d]
        @lru_cache(None)
        def normal(es):return prod(x**a for x,a in zip(ctx['image_inputs'],es)).reduce(G)
        columns=[]; rootcols=[]; labels=[]
        # The A1 coordinate (z or z^2) is independent of the A5 lattice.
        central_weight=1 if args.case=='A5' else 2
        for k in range(d//central_weight+1):
            j=d-central_weight*k;entry=source[j]
            formal=[((k,)+e if args.case=='A5' else e+(k,)) for e in entry['exponents']]
            fs=[normal(tuple(e)) for e in formal]
            D=[{tuple(e):c for e,c in f.dict().items()} for f in fs]
            A=matrix(QQ,[[f.get(tuple(e),0) for f in D] for e in target])
            columns.append(A);rootcols.append(A*entry['basis']);labels+=formal
        N=columns[0];M=rootcols[0]
        for a in columns[1:]:N=N.augment(a)
        for a in rootcols[1:]:M=M.augment(a)
        den=N.denominator();NH=matrix(ZZ,den*N).transpose().hermite_form(include_zero_rows=False).transpose()/den
        assert NH.is_square() and NH.nrows()==len(target),(d,N.dimensions(),NH.dimensions(),target,columns)
        V=NH.inverse()*M
        assert local(V)
        diagonal=matrix(ZZ,V.denominator()*V).elementary_divisors()
        vals=[int(x.valuation(3)) for x in diagonal]
        assert max(vals,default=0)<=3
        md=M.denominator()
        MH=matrix(ZZ,md*M).transpose().hermite_form(include_zero_rows=False).transpose()/md
        report['degrees'].append(dict(degree=d,rank=len(target),normal_exponents=[list(e) for e in target],
            smith_exponents=vals,nonunit_factors=[3**a for a in vals if a],
            weight_image_basis=rows(NH),root_image_basis=rows(MH),source_labels=[list(e) for e in labels]))
        print(args.case,d,'source',len(sd['exponents']),'rows',sd['coefficient_matrix'].nrows(),
              'factors',report['degrees'][-1]['nonunit_factors'],'seconds',round(perf_counter()-started,2),flush=True)
        report['seconds']=round(perf_counter()-started,3)
        output.write_text(json.dumps(report,indent=2)+'\n')
    report['result']='passed' if bound==ctx['top'] else 'partial'
    report['source_sha256']={Path(__file__).relative_to(ROOT).as_posix():sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()}
    output.write_text(json.dumps(report,indent=2)+'\n')
    print(report['result'].upper(),args.case,report['seconds'],flush=True)


if __name__=='__main__':main()
