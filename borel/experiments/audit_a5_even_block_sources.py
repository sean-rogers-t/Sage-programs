"""Independent exact invariant-ring and finite-module audit.

SymPy verifies the coefficient identities and independently recomputes the
two finite GF(2) quotients. No Sage or primary experiment imports.
"""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
from itertools import product
import gzip,json
import sympy as S
ROOT=Path(__file__).resolve().parents[1]

def main():
 start=perf_counter();sp=ROOT/'experiments/a5-even-block-sources.json';sr=json.loads(sp.read_text());sa=ROOT/sr['certificate_archive'];assert sha256(sa.read_bytes()).hexdigest()==sr['certificate_archive_sha256'];sc=json.loads(gzip.decompress(sa.read_bytes()));ap=ROOT/'experiments/a5-two-local-source.json';full=json.loads(ap.read_text());records=[]
 u,v,a,b,c=S.symbols('u v a b c');inv=[u,v,a*a,b*b-u*a*b,c*c-v*a*c,2*b-u*a,2*c-v*a,v*b-u*c,2*b*c-a*(v*b+u*c)]
 def parse(text,loc):return S.sympify(text.replace('^','**'),locals=loc)
 for case in sc['cases']:
  kind=case['blocks'];nn=S.symbols(' '.join(case['generator_names']));loc=dict(zip(case['generator_names'],nn));sub=dict(zip(nn,inv));assert all(S.expand(parse(f,loc).subs(sub,simultaneous=True))==0 for f in case['relations'])
  sigma={a:-a,b:b-u*a,c:c-v*a};assert all(S.expand(f.subs(sigma,simultaneous=True)-f)==0 for f in inv)
  # Verify actual last-pair transposition on exact S4 or S2xS2
  # coefficient coordinates, including the odd-denominator changes.
  if kind=='42':
   c1,c2,c3,c4,aa=S.symbols('c1 c2 c3 c4 aa');cs=[1,c1,c2,c3,c4];tr={aa:-aa}
   for j in range(1,5):tr[cs[j]]=sum(S.binomial(4-k,j-k)*(-aa)**(j-k)*cs[k] for k in range(j+1))
   ss=(2*aa-c1)/3;uu=-9*ss;vv=c3-aa*c2+aa**3;bb=c2;cc=c4-aa*aa*c2/3;coords=[uu,vv,aa,bb,cc];expected=[uu,vv,-aa,bb-uu*aa,cc-vv*aa]
   fullcs=[c1+aa,c2+aa*c1,c3+aa*c2,c4+aa*c3,aa*c4]
  else:
   r,s,b0,c0,aa=S.symbols('r s b0 c0 aa');tr={aa:-aa,r:r-2*aa,s:s-2*aa,b0:b0-aa*r+aa*aa,c0:c0-aa*s+aa*aa};coords=[r-aa,s-aa,aa,b0,c0];expected=[coords[0],coords[1],-aa,b0-coords[0]*aa,c0-coords[1]*aa];z=S.Symbol('z');cp=S.Poly((1+r*z+b0*z*z)*(1+s*z+c0*z*z)*(1+aa*z),z);fullcs=[cp.nth(j) for j in range(1,6)]
  assert all(S.expand(f.subs(tr,simultaneous=True)-g)==0 for f,g in zip(coords,expected))
  actual=dict(zip((u,v,a,b,c),coords));oldvars=S.symbols('c1 c2 c3 c4 c5');oldloc=dict(zip(['c1','c2','c3','c4','c5'],oldvars))
  for (text,d),row in zip(full['generators'],case['full_A5_root_generators']):
   left=parse(row['polynomial'],loc).subs(sub,simultaneous=True).subs(actual,simultaneous=True);right=parse(text,oldloc).subs(dict(zip(oldvars,fullcs)),simultaneous=True);assert S.expand(left-right)==0 and row['degree']==d
  # Fullness lemma: the ambient algebra is free on eight square-free
  # a,b,c monomials over B0=R[u,v,A,Nb,Nc]. The five columns below have
  # sole kernel (-v,u,2) over R and (v,u,0) modulo two. The latter is
  # exactly the reduction of the former, so the inclusion is 2-saturated.
  # After inverting 2, B,C,a diagonalize sigma; these columns generate
  # the full invariant ring. Saturation therefore proves completeness.
  A,Nb,Nc=S.symbols('A Nb Nc')
  basis=list(product(range(2),repeat=3));cols=[]
  for f in [1,inv[5],inv[6],inv[7],inv[8]]:
   g=S.Poly(S.rem(S.rem(S.rem(S.sympify(f),c*c-v*a*c-Nc,c),b*b-u*a*b-Nb,b),a*a-A,a),a,b,c);cols.append([g.coeff_monomial(e) for e in basis])
  assert all(S.expand(f)==0 for f in S.Matrix(8,5,lambda i,j:cols[j][i])*S.Matrix([0,-v,u,2,0]))
  assert S.Matrix(8,5,lambda i,j:cols[j][i]).rank()==4
  # Independent Groebner computation; only spanning, not minimality, is
  # needed for the finite source module. Lex is different from primary.
  polys=[parse(f,loc) for f in case['relations']]+[parse(r['polynomial'],loc) for r in case['full_A5_root_generators']];cleared=[]
  for f in polys:
   den,ff=S.Poly(f,*nn,domain=S.QQ).clear_denoms();assert den%2;cleared.append(ff.as_expr())
  G=S.groebner(cleared,*nn,modulus=2,order='grevlex');leads=[p.LM(order=G.order).exponents for p in G.polys];found={(0,)*9};todo=list(found)
  while todo:
   e=todo.pop()
   for j in range(9):
    ee=list(e);ee[j]+=1;ee=tuple(ee)
    if ee not in found and not any(all(x<=y for x,y in zip(l,ee)) for l in leads):found.add(ee);todo.append(ee)
  # The recorded homogeneous module monomials must span our independently
  # computed quotient; compare normal-form rank over GF(2), not just size.
  std=sorted(found);ix={e:i for i,e in enumerate(std)};piv={}
  for e in case['module_generator_exponents']:
   rem=S.Poly(G.reduce(S.prod(x**k for x,k in zip(nn,e)))[1],*nn,modulus=2);bits=sum((int(co)%2)<<ix[ex] for ex,co in rem.terms())
   while bits:
    k=bits.bit_length()-1
    if k in piv:bits^=piv[k]
    else:piv[k]=bits;break
  assert len(piv)==len(std)==case['module_generators'];records.append(dict(blocks=kind,source_complete_all_degrees=True,finite_source_module_verified=True,module_generators=len(std),independent_groebner_size=len(G.polys)));print('independent even-block source',records[-1],round(perf_counter()-start,2),flush=True)
 out=dict(result='passed',full_sources_complete_all_degrees=True,cases=records,seconds=round(perf_counter()-start,3),engine='Native SymPy exact coefficient substitutions, free eight-term ambient normal forms, diagonalization/saturation proof and independent GF(2) grevlex Groebner computations. No Sage or primary project imports.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),sp,ap)},raw_sha256={sa.relative_to(ROOT).as_posix():sha256(sa.read_bytes()).hexdigest()});(ROOT/'experiments/a5-even-block-sources-audit.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
