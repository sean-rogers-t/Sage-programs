"""Explicit simultaneous 2-local image bases on E7/P7.

Replays full invariant sources before applying the rational Borel map.
The historical source bases are proposals: fresh coefficient membership
and primitivity checks certify every degree, without trusting their Smith lists.
"""
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from time import perf_counter
import argparse,json
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector,prod,binomial

ROOT=Path(__file__).resolve().parents[1]
WEIGHTS=(1,2,5,6,8,9,12)


def digest(p):return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def local(A):return all(c.denominator()%2 for c in A.list())
def rows(A):return [[str(c) for c in r] for r in A]


def powers(ws,d):
    if not ws:
        if d==0:yield ()
        return
    for a in range(d//ws[0]+1):
        for e in powers(ws[1:],d-a*ws[0]):yield (a,)+e


def coefficient_inputs():
    T=PolynomialRing(QQ,7,names=('E1','E2','E3','E4','E5','E6','u'))
    es=list(T.gens()[:6]);E=es[0];u=T.gen(6)
    ps=[T(6)]
    for k in range(1,13):
        ps.append(T(sum((-1)**(j-1)*es[j-1]*ps[k-j] for j in range(1,min(k,7)))
                    +((-1)**(k-1)*k*es[k-1] if k<=6 else 0)))
    X=[sum(binomial(k,j)*6**j*(-E)**(k-j)*ps[j] for j in range(k+1)) for k in range(13)]
    Z=3*E-6*u;inputs=[(3*u-E)/2]
    for d in WEIGHTS[1:]:
        inputs.append(T(2*sum(binomial(d,j)*X[d-j]*Z**j for j in range(0,d+1,2))
                 +QQ((-1)**d)/2*(sum(binomial(d,j)*X[j]*X[d-j] for j in range(d+1))-2**d*X[d])))
    # Same symbols, but last coordinate now denotes the integral weight z.
    weight=[T(f(*(es+[2*u+E]))) for f in inputs]
    return T,inputs,weight


def ring():
    P=PolynomialRing(QQ,3,names=('h','s','t'));h,s,t=P.gens()
    relations=[s*s-2*h*t,2*s*t-9*h**4*s*s+6*h**9*s-h**14,
               t*t+10*h**3*s**3-9*h**8*s*s+2*h**13*s]
    ideal=P.ideal(relations);G=ideal.groebner_basis();assert ideal.vector_space_dimension()==56
    S=2*s;T=2*t
    basis=[(i,0,h**i) for i in range(28)]+[(i,5,h**i*S) for i in range(14)]+[(i,9,h**i*T) for i in range(14)]
    bydegree={d:[(i,b,f) for i,b,f in basis if i+b==d] for d in range(28)}
    structures={}
    for d,bs in bydegree.items():
        fs=[f.reduce(G) for i,b,f in bs];support=sorted(set().union(*(set(f.dict()) for f in fs)))
        L=matrix(QQ,[[f.dict().get(e,0) for f in fs] for e in support]);assert L.is_square() and L.det()!=0
        structures[d]=(support,L,L.inverse())
    @lru_cache(None)
    def coords(f,d):
        reduced=P(f).reduce(G)
        if d>27:assert not reduced;return vector(QQ,0)
        support,L,inv=structures[d]
        assert set(reduced.dict())<=set(support)
        return inv*vector(QQ,[reduced.dict().get(e,0) for e in support])
    # Published catalogue Borel images in its unchanged normalization.
    images=[h,-144*h*h,403200*h**5-933120*s,
            4720896*h**6-11197440*h*s,-622420992*h**8+1463132160*h**3*s,
            -1289723904*h**9-188116992*h**4*s+15237476352*t,
            3066946633728*h**12-23932781199360*h**7*s+78538306682880*h**3*t]
    return P,G,bydegree,coords,images


def certify_source(T,inputs,exps,B):
    fs=[T(prod(g**a for g,a in zip(inputs,e))) for e in exps]
    ds=[f.dict() for f in fs];support=set().union(*(set(f) for f in ds))
    cs=sorted(set(tuple(f.get(e,0) for f in ds) for e in support))
    C=matrix(QQ,cs);den=C.denominator()*B.denominator()
    K=matrix(ZZ,C.denominator()*C)*matrix(ZZ,B.denominator()*B)
    p=ZZ(2)**den.valuation(2);assert all(c%p==0 for c in K.list())
    mod=matrix(GF(2),K.nrows(),K.ncols(),[GF(2)(c//p) for c in K.list()])
    assert mod.rank()==len(exps)
    return dict(rank=len(exps),coefficient_rows=C.nrows(),integral=True,primitive_mod2=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--max-degree',type=int,default=27);args=ap.parse_args()
    start=perf_counter();T,root_inputs,weight_inputs=coefficient_inputs()
    P,G,basis,coords,images=ring();h,s,t=P.gens()
    w,p2,p5,p6,p8,p9,p12=root_inputs
    A=2*w;B=-p2/144;F5=p5/5760;F9=p9/387072
    C=(35*A*B**2-F5)/81
    D=(F9+1666*A*B**4+243*B**2*C)/19683
    assert all(local(matrix(QQ,[f.coefficients()])) for f in (A,B,F5,F9,C,D))
    a,b=2*h,h*h;c,dd=2*s,2*t
    assert ((35*a*b**2-images[2]/5760)/81-c).reduce(G)==0
    assert ((images[5]/387072+1666*a*b**4+243*b*b*c)/19683-dd).reduce(G)==0
    # Prove that the proposed simultaneous modules are algebras, independently
    # of the full-source calculation. Check multiplication by algebra generators.
    closure=[]
    for kind,gens in (('N',[(h,1),(2*s,5),(2*t,9)]),('M',[(2*h,1),(h*h,2),(2*s,5),(2*t,9)])):
        checks=0
        for d,bs in basis.items():
            for i,birth,f in bs:
                factor=2 if kind=='M' and i%2 else 1
                for g,k in gens:
                    v=coords(factor*f*g,d+k)
                    if kind=='M' and d+k<=27:
                        v=vector(QQ,[x/(2 if ii%2 else 1) for x,(ii,bb,ff) in zip(v,basis[d+k])])
                    assert all(x.denominator()%2 for x in v)
                    checks+=1
        closure.append(dict(module=kind,checks=checks))
    root_path=ROOT/'experiments/e7-character-form-results.json'
    weight_path=ROOT/'experiments/e7-adjoint-comparison-results.json'
    root=json.loads(root_path.read_text());weight=json.loads(weight_path.read_text())
    assert root['result']==weight['result']=='passed'
    out=dict(result='running',coefficient_ring='Z_(2)',degrees=[],closure_checks=closure,
             source_lifts={'2h':str(A),'h2':str(B),'2s':str(C),'2t':str(D)},
             primitive_inputs={'psi5_over_5760':str(F5),'psi9_over_387072':str(F9)})
    bound=min(27,args.max_degree)
    for d in range(bound+1):
        r=root['degrees'][d];wrow=weight['degrees'][d];exps=[tuple(e) for e in r['monomial_exponents']]
        assert set(exps)==set(powers(WEIGHTS,d)) and r['monomial_exponents']==wrow['monomial_exponents']
        nr=len(basis[d]);image=matrix(QQ,[coords(prod(g**a for g,a in zip(images,e)),d) for e in exps]).transpose()
        certificates={};matrices={}
        for kind,entry,key,index,inputs in (('M',r,'source_basis_matrix',12,root_inputs),('N',wrow,'source_basis',36,weight_inputs)):
            BB=matrix(QQ,entry[key])/(index**d)
            certificates[kind]=certify_source(T,inputs,exps,BB)
            W=image*BB
            scale=[2 if i%2 else 1 for i,b,f in basis[d]] if kind=='M' else [1]*nr
            localcoords=matrix(QQ,[[W[i,j]/scale[i] for j in range(W.ncols())] for i in range(nr)])
            assert local(localcoords) and matrix(GF(2),localcoords).rank()==nr
            matrices[kind]=dict(integral=True,primitive_mod2=True)
        multiplicity=sum(i%2 for i,b,f in basis[d])
        out['degrees'].append(dict(degree=d,basis=[dict(h_exponent=i,birth=b,polynomial=str(f),root_multiplier=2 if i%2 else 1) for i,b,f in basis[d]],
                                  source_certificates=certificates,image_certificates=matrices,multiplicity_2=multiplicity,
                                  multiplicity_4=0,multiplicity_8=0))
        print('E7 basis',d,'rank',nr,'Z/2',multiplicity,'seconds',round(perf_counter()-start,2),flush=True)
        (ROOT/'experiments/e7-two-primary-basis.json').write_text(json.dumps(out,indent=2)+'\n')
    out['result']='passed' if bound==27 else 'partial';out['seconds']=round(perf_counter()-start,3)
    out['source_sha256']={p.relative_to(ROOT).as_posix():digest(p) for p in (Path(__file__).resolve(),root_path,weight_path)}
    out['method']='Fresh Newton coefficient expansions; full-source membership and mod-2 primitivity; explicit simultaneous image bases and algebra closure in the published rational ring.'
    (ROOT/'experiments/e7-two-primary-basis.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['result'].upper(),out['seconds'],flush=True)


if __name__=='__main__':main()
