"""Replay full-source lift identities with Python integer polynomial arithmetic.

No Sage, Singular, SymPy, matrix normal form or certificate-builder import.
Forward constructions and final witnesses suffice; lattice completeness
is not assumed by this verifier.
"""
from pathlib import Path
from hashlib import sha256
from math import comb
from fractions import Fraction
import gzip,json

ROOT=Path(__file__).resolve().parents[1]
def digest(p):
    b=p.read_bytes()
    return sha256(b if p.suffix=='.gz' else b.replace(b'\r\n',b'\n')).hexdigest()
def add_scaled(out,poly,n):
    for e,c in poly.items():
        out[e]=out.get(e,0)+n*c
        if out[e]==0:del out[e]
def mul(a,b):
    out={}
    for e,x in a.items():
        for f,y in b.items():add_scaled(out,{tuple(i+j for i,j in zip(e,f)):x*y},1)
    return out
def sparse(row,size):
    assert all(type(j)==int and type(c)==int and 0<=j<size and c for j,c in row)
    assert len({j for j,c in row})==len(row)
    return dict(row)
def linear_combination(coefficients,polys):
    out={}
    for j,c in coefficients.items():add_scaled(out,polys[j],c)
    return out

def main():
    cases=[];bindings={};identity_count=0
    for m in (3,4,7):
        name='experiments/type-a-k-graded-lifts-%d.json.gz'%m
        data=json.loads(gzip.decompress((ROOT/name).read_bytes()))
        assert data['result']=='passed' and data['m']==m
        for source,expected in data['source_sha256'].items():assert digest(ROOT/source)==expected
        bindings[name]=digest(ROOT/name)
        zero=(0,)*m;cs=[{zero:1}]
        for k in range(1,m):
            e=list(zero);e[k]=1;cs.append({tuple(e):1})
        cm={}
        for k in range(m):
            e=list(zero);e[0]=m-k
            add_scaled(cm,mul(cs[k],{tuple(e):1}),(-1)**(m+1+k))
        cs.append(cm);fs=[{zero:1}]
        for j in range(1,m):
            poly={}
            for k in range(j+1):
                e=list(zero);e[0]=j-k
                add_scaled(poly,mul(cs[k],{tuple(e):1}),(-1)**k*comb(m-k,j-k))
            fs.append(poly)
        gens=cs[1:]+fs[1:];weights=list(range(1,m+1))+list(range(1,m))
        J={0:[{zero:1}]};seen=[];counts={'auxiliary_identities':0,'overflow':0,'lower':0}
        assert [r['degree'] for r in data['records']]==list(range(1,2*m-1))
        for record in data['records']:
            d=record['degree'];basis=[tuple(e) for e in record['basis']]
            assert len(set(basis))==len(basis)
            assert all(len(e)==m and all(type(x)==int and x>=0 for x in e) and sum(x*w for x,w in zip(e,[1]+list(range(1,m))))==d for e in basis)
            candidates=[];coefficient_candidates=[]
            for gi,(g,w) in enumerate(zip(gens,weights)):
                if w>d:continue
                block=[mul(g,h) for h in J[d-w]]
                candidates.extend(block)
                if gi<m:coefficient_candidates.extend(block)
            assert len(candidates)==record['candidate_count']
            assert len(coefficient_candidates)==record['coefficient_candidate_count']
            assert len(record['J'])==len(record['construction'])
            J[d]=[]
            for vector,construction in zip(record['J'],record['construction']):
                h={basis[k]:v for k,v in sparse(vector,len(basis)).items()}
                assert linear_combination(sparse(construction,len(candidates)),candidates)==h
                J[d].append(h);counts['auxiliary_identities']+=1
            expected=[(i,d-i) for i in range(1,m) if i<=d-i<m]
            assert [(x['i'],x['j']) for x in record['pairs']]==expected
            assert len(record['lifts'])==len(expected)
            for pair,lift in zip(record['pairs'],record['lifts']):
                i,j=pair['i'],pair['j'];seen.append((i,j))
                ratio=Fraction(comb(m,i)*comb(m,j),comb(m,d)) if d<m else Fraction(0)
                assert pair['denominator']==ratio.denominator and pair['numerator']==ratio.numerator
                assert ratio.denominator%data['prime']
                target={};add_scaled(target,mul(fs[i],fs[j]),m*ratio.denominator)
                if d<m:add_scaled(target,fs[d],-m*ratio.numerator)
                assert linear_combination(sparse(lift,len(coefficient_candidates)),coefficient_candidates)==target
                counts['overflow' if d>=m else 'lower']+=1
        assert sorted(seen)==[(i,j) for i in range(1,m) for j in range(i,m)]
        identity_count+=sum(counts.values());cases.append(dict(m=m,**counts,all_pairs_covered=True))
    translation_checks=0
    for m in range(2,25):
        for i in range(1,m+1):
            for r in range(i+1):
                for ell in range(i-r+1):
                    value=sum((-1)**k*comb(m-k,i-k)*comb(m-r,k-r)*comb(i-k,ell)
                              for k in range(r,i-ell+1))
                    expected=(-1)**r*comb(m-r,i-r) if ell==i-r else 0
                    assert value==expected
                    translation_checks+=1
        for k in range(m+1):
            for r in range(k+1):
                assert sum(comb(m-j,k-j)*(-1)**(j-r)*comb(m-r,j-r)
                           for j in range(r,k+1))==int(k==r)
    for name in ['tools/audit_type_a_k_graded_lifts.py','experiments/type_a_k_graded_lifts.py','TYPE_A_K_GRADED_LIFTS.md']:
        bindings[name]=digest(ROOT/name)
    out=dict(result='passed',cases=cases,integer_polynomial_identities=identity_count,
             translation_coefficient_checks=translation_checks,
             full_defect_annihilator_certified_ranks=[3,4,7],uniform_overflow_annihilation_proved=False,
             source_sha256=bindings,historical_method_date='2026-10-27',superseded_by='TYPE_A_K_INTEGRAL_RECURRENCE.md',
             scope='Independent finite replay of integer polynomial identities in the full universal source. This earlier method does not prove arbitrary-rank witnesses; the subsequent trace recurrence does. The negative flag records this method only.')
    (ROOT/'type-a-k-graded-lifts-audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='source_sha256'}))
if __name__=='__main__':main()
