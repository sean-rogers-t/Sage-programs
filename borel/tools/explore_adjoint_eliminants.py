"""Seek non-monic eliminants to certify fiberwise dimension, retaining sources."""
from pathlib import Path
import json
import time
from sage.all import QQ, PolynomialRing

ROOT=Path(__file__).resolve().parents[1]


def main(generic=False):
    data=json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())
    if generic:
        T=PolynomialRing(QQ,names='t,r'); K=T.fraction_field()
        text=data['potential']
    else:
        K=QQ; text=data['potential'].replace('t','(2)').replace('r','(3)')
    P=PolynomialRing(K,names='u,v,y,x')
    f=P(text); J=P.ideal([f.derivative(a) for a in P.gens()])
    start=time.monotonic(); print('degree reverse lex started',generic,flush=True)
    J.groebner_basis(); print('lex conversion started',time.monotonic()-start,flush=True)
    L=PolynomialRing(K,names='u,v,y,x',order='lex')
    out=J.transformed_basis('fglm',other_ring=L)
    print('lex done',time.monotonic()-start,flush=True)
    univariate=[p for p in out if set(p.variables())<=set([L.gen(3)])]
    for p in univariate:
        print('eliminant factors',p.factor(),flush=True)
    report=dict(generic=generic,elapsed_seconds=time.monotonic()-start,
                basis=list(map(str,out)),univariate=list(map(str,univariate)))
    (ROOT/('adjoint-eliminants-'+('generic' if generic else 'specialized')+'.json')).write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    import sys,resource
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    main('--generic' in sys.argv)
