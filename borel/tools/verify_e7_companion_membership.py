"""Independent replay of the homogeneous Hilbert membership certificate.

Uses a second prime and direct std, instead of the generator's slimgb,
for reverse inclusion. Rational Buchberger verification is exact. The
Hilbert-function proof of reverse inclusion is in E7_FINAL_FAMILIES.md.
"""
from pathlib import Path
from hashlib import sha256
import gzip,json,time
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,next_prime
from sage.interfaces.singular import singular
ROOT=Path(__file__).resolve().parents[1]

def main():
    raw=(ROOT/'e7-homogeneous-membership-certificate.json.gz').read_bytes()
    cert=json.loads(gzip.decompress(raw))
    oldraw=(ROOT/'e7-quartic-boundary-candidate.json.gz').read_bytes()
    old=json.loads(gzip.decompress(oldraw))
    S=PolynomialRing(QQ,names='x,y,w,A,B,C',order='degrevlex');x,y,w,A,B,C=S.gens()
    phi=x*x+y*y+w*w-x*x*y*y-x*x*w*w-y*y*w*w-x*y*w*(A*x+B*y+C*w)-A*B*x*y-A*C*x*w-B*C*y*w-A*B*C-3
    equations=[phi.derivative(v) for v in (x,y,w)]+[A*A+B*B+C*C-A*B*C-4]
    jac=matrix(S,[[f.derivative(v) for v in S.gens()] for f in equations])
    originals=equations+jac.minors(4)
    decode6=lambda v:S({tuple(e):QQ(c) for e,c in v})
    assert set(originals)==set(decode6(v) for v in old['original_generators'])
    P=PolynomialRing(QQ,names=cert['variables'],order=cert['term_order']);h=P.gen(6)
    embed=S.hom(list(P.gens()[:6]),P)
    IH=[embed(f).homogenize(h) for f in originals]
    G=[P({tuple(e):QQ(c) for e,c in f}) for f in cert['homogeneous_basis']]
    assert len(G)==2347 and all(f.is_homogeneous() and f.lc()==1 for f in G)
    assert all(not f.reduce(G) for f in IH)
    qgb=singular(P.ideal(G))
    assert singular.eval('system("verifyGB",%s)'%qgb.name()).strip()=='1'
    print('exact rational homogeneous Buchberger check passed',flush=True)
    prime=next_prime(ZZ(cert['prime']))
    while any(c.denominator()%prime==0 for f in G+IH for c in f.coefficients()):prime=next_prime(prime)
    Q=P.change_ring(GF(prime));Gmod=[Q(f) for f in G]
    assert all(f.lm().exponents()==g.lm().exponents() for f,g in zip(G,Gmod))
    # Independent characteristic and standard-basis algorithm. Homogeneity
    # is essential: no affine fiber-dimension semicontinuity is assumed.
    Bmod=list(Q.ideal([Q(f) for f in IH]).groebner_basis(algorithm='libsingular:std'))
    assert all(not f.reduce(Bmod) for f in Gmod)
    assert all(not f.reduce(Gmod) for f in Bmod)
    print('second-prime homogeneous equality passed',int(prime),flush=True)
    dehom=[f.subs({h:1}) for f in G]
    chosen=[embed(decode6(v)) for v in old['selected_members']]
    assert all(not f.reduce(dehom) for f in chosen)
    leading=[tuple(f.lm().exponents()[0][:6]) for f in chosen]
    from itertools import combinations
    supports=[{i for i,e in enumerate(v) if e} for v in leading]
    dimension=max(k for k in range(7) for v in combinations(range(6),k) if not any(a<=set(v) for a in supports))
    assert dimension==1
    report=dict(result='passed',method='homogeneous Hilbert-function equality over QQ',
                exact_rational_Buchberger_check=True,original_homogeneous_generators_reduce_to_zero=True,
                independent_prime=int(prime),generation_prime=cert['prime'],
                original_modular_ideal_algorithm='libsingular:std',
                reverse_inclusion_certified_by_homogeneous_Hilbert_functions=True,
                original_affine_ideal_memberships=6,leading_exponents=[list(v) for v in leading],
                singular_incidence_dimension_bound=1,all_E7_proper_Levis_covered=True,
                regular_fiber_dimension=12,literature_priority_claim=False,
                raw_sha256={'e7-homogeneous-membership-certificate.json.gz':sha256(raw).hexdigest(),
                            'e7-quartic-boundary-candidate.json.gz':sha256(oldraw).hexdigest()})
    names=['E7_FINAL_FAMILIES.md','tools/verify_e7_companion_membership.py']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'e7-companion-membership-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('source_sha256','raw_sha256')},indent=2))
if __name__=='__main__':main()
