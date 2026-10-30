"""Finite-field exploration of the remaining E8 A2+A1^3 singular ideal.

This deliberately reports evidence, not a characteristic-zero theorem.
The full basis is saved outside the repository to avoid versioning a
large exploratory artifact. Sage and libSingular are required.
"""
from sage.all import *
from pathlib import Path
from sage.interfaces.singular import singular
import time,json,gzip
ROOT=Path(__file__).resolve().parents[1]
P=PolynomialRing(GF(32003),names='x,y,w,a,A,B,C');x,y,w,a,A,B,C=P.gens()
f=x*x+y*y+w*w+a*a-x*x*y*y-x*x*w*w-y*y*w*w-x*x*a*a-y*y*a*a-w*w*a*a-C*(x*x*y*w+x*y*y*a+x*w*w*a+y*w*a*a)-B*(x*y*y*w+x*x*y*a+y*w*w*a+x*w*a*a)-A*(x*y*w*w+x*x*w*a+y*y*w*a+x*y*a*a)-5*x*y*w*a-B*C*(x*y+w*a)-A*C*(x*w+y*a)-A*B*(y*w+x*a)
g=[f.derivative(z) for z in (x,y,w)]+[A*A+B*B+C*C-A*B*C-4];J=matrix(P,[[q.derivative(z) for z in P.gens()] for q in g]);gens=list(set(q for q in g+J.minors(4) if q));I=P.ideal(gens);print('12357 affine start',len(gens),flush=True);t=time.time();G=I.groebner_basis(algorithm='libsingular:slimgb');print('dim',I.dimension(),'basis',len(G),'secs',time.time()-t,flush=True)
payload=dict(variables=list(P.variable_names()),prime=32003,basis=[[[list(e),int(c)] for e,c in f.dict().items()] for f in G],dimension=int(I.dimension()))
cache=ROOT.parent.parent/'e8-five-cases-cache-2026-10-26';cache.mkdir(exist_ok=True)
(cache/'e8-12357-affine-candidate.json.gz').write_bytes(gzip.compress(json.dumps(payload).encode(),mtime=0))
