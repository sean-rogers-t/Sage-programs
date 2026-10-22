"""Test a finite free-module certificate for the full root source.

Seven invariant parameters have degrees 1,2,4,6,8,10,32. If 32 invariant
products are independent over their fraction field, with degree numerator
(1+q^5)[16]_(q^2), their free span has the whole rational Hilbert series.
This is a proof only when the recorded degree/Hilbert equality succeeds.
"""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
import gzip,json,heapq
from sage.all import GF,PolynomialRing,TermOrder,matrix,prod
from e7_d5a1_source_algebra import parse
ROOT=Path(__file__).resolve().parents[1]

def standard(G,n):
 leads=[tuple(g.lm().exponents()[0]) for g in G];found={(0,)*n};todo=list(found)
 while todo:
  e=todo.pop()
  for j in range(n):
   u=tuple(a+int(i==j) for i,a in enumerate(e))
   if u not in found and not any(all(a<=b for a,b in zip(l,u)) for l in leads):found.add(u);todo.append(u)
 return sorted(found)

def main():
 start=perf_counter();path=ROOT/'experiments/e7-d5a1-sources.json';data=json.loads(path.read_text());np=ROOT/'experiments/e7-d5a1-norm.json';nd=json.loads(np.read_text());na=ROOT/nd['certificate_archive'];assert sha256(na.read_bytes()).hexdigest()==nd['certificate_archive_sha256'];nc=json.loads(gzip.decompress(na.read_bytes()));P=PolynomialRing(GF(2),names=('h','z','c2','c3','c4','c5','r'));gs=[(parse(f,P),d) for f,d in data['generators']]+[(parse(nc['norm_mod2'],P),32)];parameters=[0,2,3,7,10,12,21];assert [gs[j][1] for j in parameters]==[1,2,4,6,8,10,32]
 # Zero fiber proves finite Noether normalization. Its length is the
 # product-of-degrees ratio for this homogeneous regular sequence.
 zeroG=P.ideal([gs[j][0] for j in parameters]).groebner_basis();zero_standard=standard(zeroG,7);assert len(zero_standard)==512;print('zero fiber',len(zero_standard),round(perf_counter()-start,2),flush=True)
 K=GF(2**8,'a');a=K.gen();S=PolynomialRing(K,6,names=('z','c2','c3','c4','c5','r'),order=TermOrder('wdegrevlex',(1,2,3,4,5,2)));phi=P.hom([S.one()]+list(S.gens()),S);values=[K.one(),a,a**3,a**5,a**7,a**11,a**13];relations=[phi(gs[j][0])-v for j,v in zip(parameters[1:],values[1:])];G=S.ideal(relations).groebner_basis();mons=standard(G,6);assert len(mons)==512;ix={e:i for i,e in enumerate(mons)};gg=[phi(f).reduce(G) for f,d in gs];print('generic finite-field fiber',len(G),round(perf_counter()-start,2),flush=True)
 queue=[(0,(0,)*len(gs),S.one())];seen={(0,)*len(gs)};piv={};chosen=[]
 while queue:
  d,e,f=heapq.heappop(queue);v={ix[tuple(u)]:c for u,c in f.dict().items()}
  for j,(old,_) in sorted(piv.items()):
   if v.get(j):
    c=v[j]
    for k,b in old.items():
     v[k]=v.get(k,K.zero())-c*b
     if not v[k]:v.pop(k,None)
  if not v:continue
  j=min(v);v={k:c/v[j] for k,c in v.items()};piv[j]=(v,len(chosen));chosen.append(dict(degree=d,generator_exponents=list(e),fiber_polynomial=str(f),pivot=j));print('free candidate',len(chosen),'degree',d,round(perf_counter()-start,2),flush=True)
  for i,(g,k) in enumerate(gs):
   if i in parameters or d+k>35:continue
   u=tuple(b+int(i==j) for j,b in enumerate(e))
   if u not in seen:seen.add(u);heapq.heappush(queue,(d+k,u,(f*gg[i]).reduce(G)))
 expected=sorted(list(range(0,32,2))+list(range(5,36,2)));complete=sorted(row['degree'] for row in chosen)==expected
 cert=dict(parameter_indices=parameters,zero_fiber_groebner=list(map(str,zeroG)),zero_fiber_standard=[list(e) for e in zero_standard],field_modulus=str(K.modulus()),field_generator='a',parameter_values=list(map(str,values)),fiber_relations=list(map(str,relations)),fiber_groebner=list(map(str,G)),fiber_standard=[list(e) for e in mons],chosen_products=chosen);archive=ROOT/'experiments/e7-d5a1-free-source-certificates.json.gz';archive.write_bytes(gzip.compress(json.dumps(cert,separators=(',',':')).encode(),mtime=0));report=dict(result='passed' if complete else 'incomplete',full_source_complete_all_degrees=complete,parameter_degrees=[gs[j][1] for j in parameters],ambient_free_rank=512,source_free_rank=len(chosen),basis_degrees=[row['degree'] for row in chosen],expected_basis_degrees=expected,seconds=round(perf_counter()-start,3),scope='A positive result proves full root-source generation by the given invariant products over the parameter algebra; a negative result is not a nonfreeness theorem.',certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),ROOT/'experiments/e7_d5a1_source_algebra.py',path,np)},raw_sha256={na.relative_to(ROOT).as_posix():sha256(na.read_bytes()).hexdigest()});(ROOT/'experiments/e7-d5a1-free-source.json').write_text(json.dumps(report,indent=2)+'\n');print('COMPLETE',complete,'rank',len(chosen),'degrees',report['basis_degrees'],flush=True)

if __name__=='__main__':main()
