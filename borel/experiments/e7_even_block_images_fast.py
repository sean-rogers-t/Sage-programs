"""Complete A5 even-block images via finite source modules over A5+A2.

Use the audited 672-dimensional intermediary. The full source module is
obtained before applying Borel; no invariant-after-quotient shortcut is used.
"""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
from functools import lru_cache
import argparse,gzip,json
from sage.all import QQ,ZZ,PolynomialRing,TermOrder,matrix,prod
ROOT=Path(__file__).resolve().parents[1]
BW=(2,3,3,4,5,6)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--blocks',choices=['42','222'],required=True);args=ap.parse_args();kind=args.blocks;start=perf_counter()
 bp=ROOT/'experiments/e7-a5a2-images.json';br=json.loads(bp.read_text());ba=ROOT/br['certificate_archive'];assert sha256(ba.read_bytes()).hexdigest()==br['certificate_archive_sha256'];bc=json.loads(gzip.decompress(ba.read_bytes()))
 sp=ROOT/'experiments/a5-even-block-sources.json';sr=json.loads(sp.read_text());sa=ROOT/sr['certificate_archive'];assert sha256(sa.read_bytes()).hexdigest()==sr['certificate_archive_sha256'];sc=next(c for c in json.loads(gzip.decompress(sa.read_bytes()))['cases'] if c['blocks']==kind)
 P0=PolynomialRing(QQ,6,names=('q2','q3','e3','e4','e5','e6'),order=TermOrder('wdegrevlex',BW));baseG=[P0(f) for f in bc['groebner_basis']];norm={int(d):[tuple(e) for e in ee] for d,ee in bc['normal_exponents'].items()};ni={d:{e:i for i,e in enumerate(ee)} for d,ee in norm.items()}
 NB={0:matrix(QQ,1,1,[1])};MB={0:matrix(QQ,1,1,[1])}
 for row in bc['degrees']:
  d=row['degree'];n=len(norm[d]);NB[d]=matrix(QQ,n,n,list(map(QQ,sum(row['image_bases']['N'],[]))));MB[d]=matrix(QQ,n,n,list(map(QQ,sum(row['image_bases']['M'],[]))))
 NI={d:A.inverse() for d,A in NB.items()}
 aux=('s','t') if kind=='42' else ('s2','t2','s','t');aw=(1,2) if kind=='42' else (1,2,1,2);order=TermOrder('wdegrevlex',(1,2))
 if kind=='222':order+=TermOrder('wdegrevlex',(1,2))
 order+=TermOrder('wdegrevlex',BW);P=PolynomialRing(QQ,len(aux)+6,names=aux+P0.variable_names(),order=order);naux=len(aux);base=P.gens()[naux:];q2,q3,e3,e4,e5,e6=base;s,t=P.gens()[naux-2:naux];es=[P.one(),P.zero(),-q2,e3,e4,e5,e6];aa=[P.one()]
 for j in range(1,7):aa.append(es[j]-s*aa[j-1]-(t*aa[j-2] if j>1 else 0))
 FG=list(P.ideal(aa[5:7]).groebner_basis());outer=[(i,j) for j in range(5) for i in range(5-j)];relations=list(aa[5:7])
 if kind=='42':fiber=[e for e in outer];s1=-s;s2=s;t1=aa[2];t2=t
 else:
  s2,t2=P.gens()[:2];bb=[P.one()]
  for j in range(1,5):bb.append(aa[j]-s2*bb[j-1]-(t2*bb[j-2] if j>1 else 0))
  IG=list(P.ideal(bb[3:5]).groebner_basis());FG=IG+FG;relations+=bb[3:5];fiber=[(i,j,k,l) for k,l in outer for j in range(3) for i in range(3-j)];s1=-s-s2;t1=aa[2]-s2*s1-t2
 bmap=P0.hom(base,P);G=FG+[bmap(f) for f in baseG];assert all(g.lc()==1 for g in FG);leads=[tuple(f.lm().exponents()[0][:naux]) for f in FG];assert all(not any(g.lm().exponents()[0][naux:]) for g in FG)
 assert all(not any(all(a<=b for a,b in zip(l,e)) for l in leads) for e in fiber);assert len(fiber)==(15 if kind=='42' else 90)
 fd=[sum(a*b for a,b in zip(e,aw)) for e in fiber];fidx={e:i for i,e in enumerate(fiber)};top=45+max(fd)
 # Quadratic auxiliary w=x6, with w^2=s*w-t. This conversion is in
 # the universal weight source before imposing any Borel relation.
 U=PolynomialRing(P,'w');w=U.gen();rel=w*w-s*w+t
 if kind=='42':
  a=s-2*w;b=aa[2]+3*s*w+6*w*w;c3=aa[3]-2*w*aa[2]-3*s*w*w-4*w**3;c4=aa[4]-w*aa[3]+w*w*aa[2]+s*w**3+w**4;u=-9*s;v=c3-a*b+a**3;c=c4-a*a*b/3
 else:u=2*s1+s2;v=s1+2*s2;a=s-2*w;b=t1-t-u*w;c=t2-t-v*w
 vals=[u,v,a*a,b*b-u*a*b,c*c-v*a*c,2*b-u*a,2*c-v*a,v*b-u*c,2*b*c-a*(v*b+u*c)];inv=[]
 for f in vals:
  f=U(f).mod(rel);assert f.degree()<=0;inv.append(P(f[0]))
 S=PolynomialRing(QQ,9,names=sc['generator_names']);phi=S.hom(inv,P)
 assert all(not phi(S(f)) for f in sc['relations'])
 for row,(f,d) in zip(sc['full_A5_root_generators'],bc['generators']['M'][2:]):assert not (phi(S(row['polynomial']))-bmap(P0(f))).reduce(G)
 modules=[]
 for j,e in enumerate(sc['module_generator_exponents']):
  f=P(prod(inv[i]**k for i,k in enumerate(e) if k)).reduce(G);components={}
  for ex,c0 in f.dict().items():
   fi=fidx[tuple(ex[:naux])];components[fi]=components.get(fi,P0.zero())+P0({tuple(ex[naux:]):c0})
  modules.append((sc['module_generator_degrees'][j],components))
 print('complete module/fiber',kind,len(modules),len(fiber),round(perf_counter()-start,2),flush=True)
 # A cached rational monomial reducer avoids Singular conversions in
 # every small multiplication map.
 reducers=[]
 for f in baseG:
  e=tuple(f.lm().exponents()[0]);reducers.append((e,{tuple(u):-c/f.lc() for u,c in f.dict().items() if tuple(u)!=e}))
 @lru_cache(None)
 def normal_mon(e):
  if sum(a*b for a,b in zip(e,BW))>45:return {}
  for lead,tail in reducers:
   if all(a>=b for a,b in zip(e,lead)):
    shift=tuple(a-b for a,b in zip(e,lead));out={}
    for u,c in tail.items():
     for v,b in normal_mon(tuple(a+b for a,b in zip(u,shift))).items():out[v]=out.get(v,QQ.zero())+c*b
    return {e:c for e,c in out.items() if c}
  return {e:QQ.one()}
 # Cache small base multiplication maps. All rational calculations take
 # place in rank <=30 homogeneous pieces, never in the large total ring.
 @lru_cache(None)
 def block(ftext,bd):
  f=P0(ftext);k=sum(a*b for a,b in zip(f.exponents()[0],BW));outd=bd+k;n=len(norm[outd]);A=matrix(QQ,n,len(norm[bd]))
  for j,e in enumerate(norm[bd]):
   for ex,c0 in f.dict().items():
    for v0,b0 in normal_mon(tuple(a+b for a,b in zip(ex,e))).items():A[ni[outd][tuple(v0)],j]+=c0*b0
  T=NI[outd]*A*MB[bd];assert T.denominator()%2;return [[int(c.numerator()*pow(int(c.denominator()),-1,4))%4 for c in col] for col in T.columns()]
 def ranks4(columns,n):
  # Mod-four column operations; high plane records the even remainder.
  piv={};even=[];mask=(1<<n)-1
  for lo,hi in columns:
   while lo:
    k=lo.bit_length()-1
    if k not in piv:
     if hi>>k&1:hi^=lo
     piv[k]=(lo,hi);break
    pl,ph=piv[k]
    if hi>>k&1:hi^=ph^(lo&pl)
    else:hi^=ph^((~lo)&pl&mask)
    lo^=pl
   else:even.append(hi)
  ep={}
  for hi in even:
   for k,(pl,_) in sorted(piv.items(),reverse=True):
    if hi>>k&1:hi^=pl
   while hi:
    k=hi.bit_length()-1
    if k in ep:hi^=ep[k]
    else:ep[k]=hi;break
  return len(piv),len(ep)
 records=[];columnhash=[]
 for d in range(top+1):
  offsets={};n=0
  for fi,k in enumerate(fd):
   if 0<=d-k<=45:offsets[fi]=n;n+=len(norm[d-k])
  columns=[]
  for md,comp in modules:
   bd=d-md
   if not 0<=bd<=45:continue
   localcols=[[0,0] for _ in norm[bd]]
   for fi,f in comp.items():
    if fi not in offsets:continue
    outd=d-fd[fi]
    if not norm[outd]:continue
    for j,cc in enumerate(block(str(f),bd)):
     shift=offsets[fi]
     localcols[j][0]|=sum((c&1)<<(i+shift) for i,c in enumerate(cc));localcols[j][1]|=sum((c>>1)<<(i+shift) for i,c in enumerate(cc))
   columns.extend(localcols)
  r0,r1=ranks4(columns,n);assert r0+r1==n,(kind,d,n,r0,r1,'higher factors require a new calculation');m=n-r0;record=dict(degree=d,rank=n,multiplicity_2=m,mod_four_unit_rank=r0,mod_four_even_rank=r1,columns=len(columns));records.append(record);columnhash.append(sha256(json.dumps(columns,separators=(',',':')).encode()).hexdigest());print('full even-block image',kind,record,round(perf_counter()-start,2),flush=True)
 cert=dict(blocks=kind,fiber_exponents=fiber,fiber_degrees=fd,fiber_relations=list(map(str,relations)),fiber_groebner_basis=list(map(str,FG)),partial_root_generators=list(map(str,inv)),module_components=[dict(degree=d,components={str(i):str(f) for i,f in cc.items()}) for d,cc in modules],mod_four_column_sha256=columnhash,degrees=records)
 stem='e7-even-block-fast-'+kind;archive=ROOT/('experiments/'+stem+'-certificates.json.gz');archive.write_bytes(gzip.compress(json.dumps(cert,separators=(',',':')).encode(),mtime=0));out=dict(result='passed',complete=True,blocks=kind,total_rank=sum(r['rank'] for r in records),top_polynomial_degree=top,degrees=records,seconds=round(perf_counter()-start,3),scope='Complete primary integral source images for the indicated A5 block Levi times the complementary A2; common A2 fibers transfer the results to standard E7 parabolics. Separate source, lattice and image audits are required.',certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),bp,sp)},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (ba,sa)});(ROOT/('experiments/'+stem+'.json')).write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
