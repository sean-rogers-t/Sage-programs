"""Homogeneous Hilbert-function certificate for exact original membership."""
from pathlib import Path
import gzip,json,time
from sage.all import QQ,GF,PolynomialRing,ZZ,next_prime
from sage.interfaces.singular import singular
ROOT=Path(__file__).resolve().parents[1]
d=json.loads(gzip.decompress((ROOT/'e7-quartic-boundary-candidate.json.gz').read_bytes()))
P=PolynomialRing(QQ,names=d['variables']+['h'],order='degrevlex');h=P.gen(6)
decode=lambda v:P({tuple(e)+(0,):QQ(c) for e,c in v})
original=[decode(v) for v in d['original_generators']]
IH=[f.homogenize(h) for f in original]
t=time.time();singular.lib('modstd.lib')
print('homogeneous modular basis starting',flush=True)
cache=ROOT/'e7-homogeneous-basis-candidate.json.gz'
encode=lambda f:[[list(map(int,e)),str(c)] for e,c in f.dict().items()]
complete=ROOT/'e7-homogeneous-membership-certificate.json.gz'
if complete.exists():
    G=[P({tuple(e):QQ(c) for e,c in f}) for f in json.loads(gzip.decompress(complete.read_bytes()))['homogeneous_basis']]
elif cache.exists():
    G=[P({tuple(e):QQ(c) for e,c in f}) for f in json.loads(gzip.decompress(cache.read_bytes()))]
else:
    G=list(singular.modStd(singular(P.ideal(IH)),0).sage().gens())
    cache.write_bytes(gzip.compress((json.dumps([encode(f) for f in G],separators=(',',':'))+'\n').encode(),mtime=0))
print('homogeneous basis',len(G),'seconds',time.time()-t,flush=True)
# Independent a posteriori certification: homogeneous Hilbert functions,
# exact rational Buchberger test and one finite-field reverse inclusion.
assert all(f.is_homogeneous() for f in G)
assert all(not f.reduce(G) for f in IH)
Gsing=singular(P.ideal(G))
assert singular.eval('system("verifyGB",%s)'%Gsing.name()).strip()=='1'
print('rational standard-basis check passed',flush=True)
prime=ZZ(32003)
while any(c.denominator()%prime==0 for f in G+IH for c in f.coefficients()):prime=next_prime(prime)
Q=P.change_ring(GF(prime));Gmod=[Q(f) for f in G]
assert all(tuple(f.lm().exponents()[0])==tuple(g.lm().exponents()[0]) for f,g in zip(G,Gmod))
Bmod=list(Q.ideal([Q(f) for f in IH]).groebner_basis(algorithm='libsingular:slimgb'))
assert all(not f.reduce(Bmod) for f in Gmod)
assert all(not f.reduce(Gmod) for f in Bmod)
print('homogeneous Hilbert-function equality certified',flush=True)
dehom=[P(f.subs({h:1})) for f in G]
chosen=[decode(v) for v in d['selected_members']]
assert all(not f.reduce(dehom) for f in chosen)
encode=lambda f:[[list(map(int,e)),str(c)] for e,c in f.dict().items()]
payload=dict(variables=list(P.variable_names()),term_order='degrevlex',
             homogeneous_basis=[encode(f) for f in G],prime=int(prime),
             modular_basis=[[[list(map(int,e)),int(c)] for e,c in f.dict().items()] for f in Bmod])
(ROOT/'e7-homogeneous-membership-certificate.json.gz').write_bytes(gzip.compress((json.dumps(payload,separators=(',',':'))+'\n').encode(),mtime=0))
print('six exact memberships saved',flush=True)
