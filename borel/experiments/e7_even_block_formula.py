"""One exact even-block formula for the completed A5-based E7 cases.

This verifies a finite-family polynomial identity. A geometric folding
explanation or an analogous theorem for other ranks is not asserted.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
import gzip,json
ROOT=Path(__file__).resolve().parents[1]

def trim(a):
 a=list(a)
 while len(a)>1 and not a[-1]:a.pop()
 return a
def add(a,b,shift=0,scale=1):
 out=list(a)+[0]*max(0,len(b)+shift-len(a))
 for i,c in enumerate(b):out[i+shift]+=scale*c
 return trim(out)
def mul(a,b):
 out=[0]*(len(a)+len(b)-1)
 for i,c in enumerate(a):
  for j,d in enumerate(b):out[i+j]+=c*d
 return trim(out)
def bracket(n,step=1):return [int(i%step==0) for i in range((n-1)*step+1)]
@lru_cache(None)
def gaussian(n,k):
 if k==0 or k==n:return (1,)
 return tuple(add(gaussian(n-1,k),gaussian(n-1,k-1),n-k))
def universal(P):
 numerator=[0]+list(P);out=[]
 for i,c in enumerate(numerator):out.append(c-(out[i-1] if i else 0))
 assert not out[-1];return trim(out[:-1])

def main():
 ap=ROOT/'experiments/e7-a5a2-images-audit.json';audit=json.loads(ap.read_text());assert audit['result']=='passed' and audit['complete'];bp=ROOT/'experiments/e7-a5a2-images.json';br=json.loads(bp.read_text());ba=ROOT/br['certificate_archive'];bc=json.loads(gzip.decompress(ba.read_bytes()));P0=[len(bc['normal_exponents'][str(d)]) for d in range(46)];T0=[0]+[row['smith_multiplicities'].get('1',0) for row in bc['degrees']];U0=universal(P0);C0=add(U0,T0,scale=-1);cases=[];paths=[ap,bp]
 for name,composition,fiber,folded in [('42',[2,1],gaussian(6,2),bracket(3,4)),('222',[1,1,1],mul(gaussian(6,2),gaussian(4,2)),mul(bracket(2,4),bracket(3,4)))]:
  pp=ROOT/('experiments/e7-even-block-fast-'+name+'.json');primary=json.loads(pp.read_text());cp=ROOT/('experiments/e7-even-block-'+name+'-comparison-audit.json');checked=json.loads(cp.read_text());assert checked['result']=='passed' and checked['complete'] and checked['primary_and_independent_all_columns_agree'];P=[row['rank'] for row in primary['degrees']];T=checked['multiplicity_coefficients'];assert trim(P)==mul(P0,fiber);assert trim(T)==add(mul(U0,fiber),mul(C0,folded),scale=-1);paths +=[pp,cp];cases.append(dict(blocks=name,halved_composition=composition,ordinary_fiber_polynomial=list(fiber),halved_block_polynomial=folded,multiplicity_coefficients=T,formula_verified=True))
 # The same factors are exact Weyl length-polynomial ratios of C3
 # and its maximal-rank C-block subgroups, evaluated at q^2.
 WC3=mul(mul(bracket(2,2),bracket(4,2)),bracket(6,2));WC2C1=mul(mul(bracket(2,2),bracket(4,2)),bracket(2,2));WC1cubed=mul(mul(bracket(2,2),bracket(2,2)),bracket(2,2));assert mul(WC2C1,cases[0]['halved_block_polynomial'])==WC3;assert mul(WC1cubed,cases[1]['halved_block_polynomial'])==WC3
 out=dict(result='passed',uniform_even_block_formula_verified=True,ambient_type='E7',intermediary='A5+A2',halved_rank=3,base_rank_polynomial=P0,base_multiplicity_polynomial=T0,base_universal_polynomial=U0,base_correction_polynomial=C0,cases=cases,formula='T_a(q) = q/(1+q) P0(q) [6;2a1,...,2ak]_q - C0(q) [3;a1,...,ak]_(q^4), C0=qP0/(1+q)-T0',A2_fiber_multipliers=dict(partial=[1,1,1],full=[1,2,2,1]),folded_Weyl_polynomial_identity_verified=True,scope='An exact common formula for all compositions of three, hence every even-block A5-based case in this E7 family, together with the common A2 fiber transfers. The signed correction is not a Hilbert series. The Weyl-polynomial identity alone gives no geometric folding map or theorem for other ranks; no priority claim.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in [Path(__file__).resolve()]+paths},raw_sha256={ba.relative_to(ROOT).as_posix():sha256(ba.read_bytes()).hexdigest()});(ROOT/'experiments/e7-even-block-formula.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED uniform E7 even-block formula',flush=True)

if __name__=='__main__':main()
