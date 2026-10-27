"""Reproduce a homogeneous dimension certificate, never sample a target fiber.

Run separately at 32003 and 32009. A truncated standard basis suffices:
every returned polynomial belongs to the original homogeneous ideal.
Only 41 selected members and their leading supports are needed for the bound.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
import gzip,json,sys,time
from sage.all import QQ,GF,PolynomialRing,matrix,prod
from sage.interfaces.singular import singular
ROOT=Path(__file__).resolve().parents[1]

def main():
    prime=int(sys.argv[1]) if len(sys.argv)>1 else 32003
    assert prime in (32003,32009)
    name=f'e8-a1-four-{prime}-certificate.json.gz'
    cert=json.loads(gzip.decompress((ROOT/name).read_bytes()))
    row=next(r for r in json.loads((ROOT/'e8-final-character-models.json').read_text()) if r['seed']==[1,2,5,7])
    Q=PolynomialRing(QQ,names=row['residual_variables']+['h'],order='degrevlex');v=Q.gens();h=v[-1]
    f=sum(QQ(q['coefficient'])*prod(w**e for w,e in zip(v,q['exponent'])) for q in row['residual_terms'])
    g=[f.derivative(w) for w in v[:4]];J=matrix(Q,[[q.derivative(w) for w in v[:-1]] for q in g])
    original=[q.homogenize(h) for q in g+J.minors(4) if q]
    P=Q.change_ring(GF(prime));decode=lambda q:P({tuple(e):c for e,c in q})
    assert set(map(P,original))==set(decode(q) for q in cert['original_homogeneous'])
    I=singular(P.ideal([P(q) for q in original]));singular.eval('degBound=12;')
    start=time.time();print('Replaying original homogeneous ideal at prime',prime,flush=True)
    G=list(I.std().sage().gens());singular.eval('degBound=0;')
    chosen=[decode(q['polynomial']) for q in cert['selected_members']]
    assert all(not q.reduce(G) for q in chosen)
    supports=[]
    for q,record in zip(chosen,cert['selected_members']):
        lm=tuple(q.lm().exponents()[0]);assert lm==tuple(record['leading_exponent'])
        supports.append(set(i for i,e in enumerate(lm) if e))
    # A coordinate subspace of a monomial variety is allowed iff it
    # contains no generator support. Enumerate these independently of CAS dimension.
    allowed=[]
    for k in range(10):
        allowed += [tuple(c) for c in combinations(range(9),k) if all(not s<=set(c) for s in supports)]
    assert max(map(len,allowed))==4
    assert P.ideal([q.lm() for q in chosen]).dimension()==4
    result=dict(result='passed',prime=prime,selected_members=41,partial_basis_size=len(G),
                degree_bound=12,original_ideal_memberships_replayed=True,
                homogeneous_dimension_upper_bound=4,characteristic_zero_affine_dimension_upper_bound=3,
                coordinate_supports_enumerated=512,seconds=round(time.time()-start,3),
                source_sha256={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for n in ['tools/verify_e8_a1_four_certificate.py','e8-final-character-models.json']},
                raw_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest()})
    (ROOT/f'e8-a1-four-{prime}-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','raw_sha256')},indent=2))
if __name__=='__main__':main()
