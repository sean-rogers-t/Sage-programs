"""Exact Schur reduction; no determinant-zero branch is discarded."""
from pathlib import Path
import json,time
from sage.all import QQ,PolynomialRing,matrix,vector
ROOT=Path(__file__).resolve().parents[1]

def main():
    T=PolynomialRing(QQ,names='t,r');K=T.fraction_field();t,r=map(K,T.gens())
    P=PolynomialRing(K,names='x,a,b,y');x,a,b,y=P.gens()
    data=json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())
    text=data['potential'].replace('u','(a/t)').replace('v','(t*b)')
    f=P(text)
    H=matrix(P,[[f.derivative(z,w) for w in (a,b)] for z in (a,b)])
    L=vector(P,[f.derivative(z).subs({a:0,b:0}) for z in (a,b)])
    Q=f.subs({a:0,b:0});D=H.det()
    N=Q*D-(L*H.adjugate()*L)/2
    print('degrees D,N',D.total_degree(),N.total_degree(),flush=True)
    print('N factor',N.factor(),flush=True)
    E=[(N.derivative(z)*D-N*D.derivative(z)) for z in (x,y)]
    for z,e in zip((x,y),E): print('critical numerator',z,'degree',e.total_degree(),'factor',e.factor(),flush=True)
    report=dict(potential=str(f),matrix=[[str(a) for a in row] for row in H.rows()],
                linear=list(map(str,L)),constant=str(Q),determinant=str(D),numerator=str(N),
                critical_numerators=list(map(str,E)))
    (ROOT/'adjoint-quadratic-reduction.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
