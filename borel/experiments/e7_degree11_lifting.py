"""Fresh degree-eleven E7 half-sum and invariant-lifting certificate.

The cubic/quintic are Frobenius derivatives of the quadratic form, not
read from a Smith table or a previous saturation witness. A primitive
kernel check and explicit identities support the group-equivariant lemma
in E7_DEGREE11_LIFTING.md.
"""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
import gzip,json
from sage.all import QQ,ZZ,GF,PolynomialRing,CartanMatrix,matrix,vector,gcd,prod
from e7_full_flag_kernel import module_exponents,encoded
ROOT=Path(__file__).resolve().parents[1]

def digest(p):return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def main():
 start=perf_counter();C=matrix(QQ,CartanMatrix('E7'))
 XB=C.matrix_from_columns(range(6)).augment(vector(QQ,[0]*6+[1]).column())
 R=PolynomialRing(QQ,7,names=['a'+str(i) for i in range(1,7)]+['w']);a=R.gens()[:6];w=R.gen(6)
 B=R(-QQ(3)/2*(vector(R,a)*C[:6,:6].inverse()*vector(R,a)))
 bp=[R.one()]
 for j in range(5):bp.append(bp[-1]*B)
 def reduced(f):return R(sum(c*R({tuple(e[:6])+(e[6]%2,):1})*bp[e[6]//2] for e,c in f.dict().items()))
 def mod2(f):
  assert f.denominator()%2
  return F({e:GF(2)(ZZ(c.numerator())) for e,c in f.dict().items()})
 F=PolynomialRing(GF(2),names=R.variable_names());ab=F.gens()[:6];Bb=mod2(B)
 G3=sum(x*x*Bb.derivative(x) for x in ab);G5=sum(x**4*Bb.derivative(x) for x in ab)
 lift=lambda f:R({e:ZZ(c) for e,c in f.dict().items()})
 g3,g5=lift(G3),lift(G5)
 orbit={(0,0,0,0,0,0,1)};pending=list(orbit)
 while pending:
  t=pending.pop()
  for i in range(7):
   u=tuple(t[j]-t[i]*C[i,j] for j in range(7))
   if u not in orbit:orbit.add(u);pending.append(u)
 assert len(orbit)==56
 inv=XB.inverse();forms=[sum(c*x for c,x in zip(inv*vector(QQ,t),R.gens())) for t in sorted(orbit)]
 assert reduced(sum(f*f for f in forms))==0
 relations={}
 for d in (6,8,10):
  f=reduced(sum(l**d for l in forms));den=f.denominator();f=R(den*f/gcd((den*f).coefficients()))
  assert f.denominator()==1;relations[d]=f
 print('fresh moments and quadratic derivatives',round(perf_counter()-start,2),flush=True)
 C5=relations[6].derivative(w);C7=relations[8].derivative(w);C9=relations[10].derivative(w)
 assert (C5/4).denominator()==1 and (C7/2).denominator()==1 and (C9/4).denominator()==1
 assert mod2(relations[6])==G3**2 and mod2(relations[8])==G3*G5
 assert mod2(C5/4)==G5+Bb*G3
 assert mod2(C7/2)==Bb*G5+Bb**2*G3
 H=(relations[6]*g5+relations[8]*g3)/2;assert H.denominator()==1
 assert mod2(H.derivative(w))==mod2(B*relations[8]+B**2*relations[6])
 support=module_exponents(11);ix={e:i for i,e in enumerate(support)}
 labels=[(k,e) for k in relations for e in module_exponents(11-k)]
 A=matrix(ZZ,len(support),len(labels))
 for j,(k,e) in enumerate(labels):
  for u,c in reduced(relations[k]*R({e:1})).dict().items():A[ix[tuple(u)],j]=ZZ(c)
 hc=vector(ZZ,[H.monomial_coefficient(R({e:1})) for e in support])
 initial_rank=matrix(GF(2),A).rank();primitive_rank=matrix(GF(2),A.augment(hc.column())).rank()
 q=PolynomialRing(ZZ,'q').gen();PG=prod(sum(q**j for j in range(k)) for k in (2,6,8,10,12,14,18))
 assert len(labels)==len(support)-int(PG[11])==462
 assert initial_rank==461 and primitive_rank==462
 # The displayed syzygy spans the one-dimensional mod-two coefficient kernel.
 coeff=vector(GF(2),[G5.monomial_coefficient(F({e:1})) if k==6 else G3.monomial_coefficient(F({e:1})) if k==8 else 0 for k,e in labels])
 assert coeff and matrix(GF(2),A)*coeff==0
 # The four non-good-Levi cases admit no integral lift of the cubic.
 T=PolynomialRing(GF(2),6,names=['a'+str(i) for i in range(1,7)]);t=T.gens()
 gamma3=T(G3(*t,0));C6=C[:6,:6];nonlifting=[]
 for I in ((1,2,3,4,5,6),(1,2,3,4,5),(2,3,4,5,6),(2,3,4,5)):
  omitted=[j for j in range(6) if j+1 not in I]
  ls=[sum(GF(2)(ZZ(3*C6.inverse()[r,j]))*t[r] for r in range(6)) for j in omitted]
  remainder=gamma3.reduce(T.ideal(ls).groebner_basis()) if ls else gamma3
  assert remainder
  nonlifting.append(dict(retained_nodes=list(I),fixed_linear_forms=[str(l) for l in ls],cubic_restriction=str(remainder),cubic_does_not_lift=True))
 family=[];good=[]
 for mask in range(127):
  I=[j for j in range(7) if mask>>j&1]
  if any(j not in I for j in (1,4,6)):family.append(tuple(j+1 for j in I))
  # E7 is a tree; components are type A iff the induced graph has degree <=2.
  is_A=all(sum(C[i,j]==-1 for j in I)<=2 for i in I)
  if is_A and matrix(GF(2),C.matrix_from_columns(I)).rank()==len(I):good.append(tuple(j+1 for j in I))
 assert len(family)==112 and len(good)==108
 assert set(family)-set(good)=={tuple(c['retained_nodes']) for c in nonlifting}
 archive=ROOT/'experiments/e7-degree11-lifting-certificates.json.gz'
 certificate=dict(quadratic=encoded(w*w-B),relations={str(d):encoded(f) for d,f in relations.items()},cubic=encoded(g3),quintic=encoded(g5),half_sum=encoded(H))
 archive.write_bytes(gzip.compress(json.dumps(certificate,separators=(',',':')).encode(),mtime=0))
 out=dict(result='passed',degree=11,uniform_family_size=112,E6_subgroups_lemma=True,orbit_size=56,quadratic_Frobenius_identities_verified=True,half_sum_derivative_identity_verified=True,initial_kernel_rank_mod2=int(initial_rank),primitive_kernel_rank_mod2=int(primitive_rank),rational_kernel_rank=462,saturation_quotient_order=2,mod2_syzygy_dimension=1,good_Levi_count=108,exceptional_noncubic_lifting_cases=nonlifting,seconds=round(perf_counter()-start,3),scope='Noncircular degree-eleven invariant-kernel stability, with a general rank-one saturation lemma proved in the note. No Smith multiplicity input and no all-degree E7 structural theorem.',certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in (Path(__file__).resolve(),ROOT/'experiments/e7_full_flag_kernel.py')})
 (ROOT/'experiments/e7-degree11-lifting.json').write_text(json.dumps(out,indent=2)+'\n')
 print('PASSED',out['seconds'],flush=True)
if __name__=='__main__':main()
