"""Independent Sage root-lattice evaluation and differentiated linear solves."""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import RootSystem, QQ, PolynomialRing, matrix, vector

ROOT=Path(__file__).resolve().parents[1]


def root_checks(kind):
    lattice=RootSystem([kind[0],int(kind[1:])]).root_lattice()
    roots=list(lattice.roots());positive=list(lattice.positive_roots());rank=int(kind[1:])
    ambient={a:a.to_ambient() for a in roots}
    scale=QQ(2)/max(v.scalar(v) for v in ambient.values())
    dot=lambda a,b:scale*ambient[a].scalar(ambient[b])
    theta=lattice.highest_root()
    eigen=dot(theta,theta)+sum(dot(theta,a) for a in positive)
    simple=lattice.simple_roots();added=set()
    if kind=='E6':
        for a in roots:
            labels=tuple(a.scalar(simple[i].associated_coroot()) for i in (1,2,3,5,6))
            charge=sum(n*a.scalar(simple[i].associated_coroot()) for i,n in enumerate((2,3,4,6,4,2),1))
            if labels in ((0,1,0,0,0),(0,-1,0,0,0)) and abs(charge)==3:added.add(a)
        assert len(added)==4
    defect=None
    for weighted in (False,True):
        e={a:QQ(2)**sum((i if weighted else 1)*n for i,n in a.monomial_coefficients().items()) for a in roots}
        assert all(e[a]!=1 for a in positive)
        def residual(altered):
            mult={a:1+int(altered and a in added) for a in roots}
            flat=sum(dot(a,a)*mult[a]*e[a] for a in roots)
            # Direct rational evaluation of the radial differential operator;
            # no root-string division or native root coordinates are used.
            drift=sum((e[a]+1)/(e[a]-1)*sum(dot(a,b)*mult[b]*e[b] for b in roots) for a in positive)
            return flat+drift-eigen*(rank+sum(mult[a]*e[a] for a in roots))
        assert residual(False)==0
        if kind=='E6' and not weighted:
            defect=residual(True);assert defect!=0
    report=dict(root_count=len(roots),casimir_eigenvalue=str(eigen),
        root_length_squares=sorted({str(dot(a,a)) for a in roots}),exact_torus_evaluations=2)
    if defect is not None:report['altered_Levi_radial_defect']=str(defect)
    return report


def schur_checks():
    T=PolynomialRing(QQ,'t');t=T.gen();F=T.fraction_field()
    trace=t+F(1)/t
    Vtrace=(t*t+1)/(t*t-1)*t*trace.derivative(t)
    assert Vtrace==trace
    C=PolynomialRing(QQ,'x');cx=C.gen()
    assert cx*(cx*cx-1).derivative(cx)==2*(cx*cx-1)+2
    # Specialize only the central coordinates, then reconstruct the block and
    # differentiate before substituting x,y. No stored Schur matrix is read.
    P=PolynomialRing(QQ,names='x,a,b,y');x,a,b,y=P.gens()
    raw=json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())['potential']
    raw=raw.replace('u','(a/2)').replace('v','(2*b)').replace('t','(2)').replace('r','(3)')
    f=P(raw);H=matrix(P,[[f.derivative(i,j) for j in (a,b)] for i in (a,b)])
    L=vector(P,[f.derivative(i)(a=0,b=0) for i in (a,b)])
    Q=f(a=0,b=0);D=H.det();star=-H.change_ring(P.fraction_field()).solve_right(L)
    R=f(a=star[0],b=star[1]);N=P(R*D)
    wedge=R.derivative(x)*D.derivative(y)-R.derivative(y)*D.derivative(x)
    assert wedge(x=1,y=0)==QQ(515)/6
    assert (Q.derivative(x)*D.derivative(y)-Q.derivative(y)*D.derivative(x))(x=1,y=0)==-QQ(525)/2
    assert D(x=1,y=0)==-QQ(9)/4 and N(x=1,y=0)==QQ(27)/4
    assert D(x=-QQ(1)/2,y=0)==0 and N(x=-QQ(1)/2,y=0)==12
    assert D.derivative(x)(x=-QQ(1)/2,y=0)==6
    hh=H.apply_map(lambda p:QQ(p(x=-QQ(1)/2,y=0)))
    ll=vector(QQ,[p(x=-QQ(1)/2,y=0) for p in L])
    assert hh.rank()==1 and hh.augment(ll).rank()==2
    P=PolynomialRing(QQ,names='h,j,k,l,m,q,dh,dj,dk,dl,dm,dq,tau')
    h,j,k,l,m,q,dh,dj,dk,dl,dm,dq,tau=P.gens();F=P.fraction_field()
    H=matrix(F,[[h,j],[j,k]]);Hd=matrix(F,[[dh,dj],[dj,dk]])
    L=vector(F,[l,m]);Ld=vector(F,[dl,dm])
    star=-H.solve_right(L)
    Ltau=L+tau*Ld;Rtau=q+tau*dq-(Ltau*(H+tau*Hd).solve_right(Ltau))/2
    assert Rtau.derivative(tau)(tau=0)==dq+star*Ld+(star*Hd*star)/2
    return dict(E8_R_D_wedge='515/6',E8_Q_D_wedge='-525/2',
        E8_simple_pole_verified=True,E8_pole_point_quadratic_equations_inconsistent=True,
        universal_Schur_envelope_derivative_verified=True,
        SL2_radial_derivation_regular_on_quotient_discriminant=True)


def main():
    native=json.loads((ROOT/'adjoint-radial-schur-audit.json').read_text())
    records={kind:root_checks(kind) for kind in ('A3','B3','C3','D4','G2','F4','E6','E7','E8')}
    for kind,row in records.items():
        for key in ('root_count','casimir_eigenvalue','root_length_squares'):
            assert row[key]==native['root_identities'][kind][key]
    assert records['E6']['altered_Levi_radial_defect']==native['root_identities']['E6']['altered_Levi_radial_defect']
    schur=schur_checks()
    for key in ('E8_R_D_wedge','E8_Q_D_wedge'):assert schur[key]==native['schur'][key]
    report=dict(result='passed',engine='Sage ambient root lattices, rational evaluation and linear solves',
        root_identities=records,schur=schur,independent_of_native_root_strings=True,
        general_singular_image_bound_proved=False,new_Levi_coverage=False,
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
            for p in (Path(__file__),ROOT/'adjoint-radial-schur-audit.json',
                      ROOT/'adjoint-nilpotent-source-screen.json',ROOT/'UNIFORM_ADJOINT_PROOF_TARGET.md')})
    (ROOT/'adjoint-radial-schur-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
