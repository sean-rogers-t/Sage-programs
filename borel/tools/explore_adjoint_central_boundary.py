"""Exploration of the remaining central boundary of ONE E8 candidate.

Denominators of a generic Groebner basis are leads, not yet a certified
exceptional locus for the original incidence ideal. No target assertion
is inferred from a specialized central source fiber.
"""
from pathlib import Path
import json
import time
import argparse
import gzip
from hashlib import sha256
from itertools import combinations
from sage.all import QQ, PolynomialRing, matrix, NumberField

ROOT=Path(__file__).resolve().parents[1]


def potential():
    T=PolynomialRing(QQ,names='t,r'); K=T.fraction_field()
    P=PolynomialRing(K,names='x,u,v,y')
    data=json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())
    return T,P,P(data['potential'])


def main():
    T,P,f=potential()
    start=time.monotonic()
    J=P.ideal([f.derivative(x) for x in P.gens()])
    gb=J.groebner_basis()
    denominators={c.denominator() for g in gb for c in g.coefficients()}
    print('basis ready; distinct denominators',len(denominators),flush=True)
    factors=set()
    for d in denominators:
        for h,m in d.factor():
            factors.add(h/h.lc())
    report=dict(status='exploratory_denominators_only',
        denominator_factors=sorted(map(str,factors)),
        original_ideal_boundary_certified=False,
        elapsed_seconds=time.monotonic()-start)
    (ROOT/'adjoint-central-boundary-exploration.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


def slice_test():
    T=PolynomialRing(QQ,names='t'); K=T.fraction_field();t=K.gen()
    P=PolynomialRing(K,names='x,u,v,y,r,z'); x,u,v,y,r,z=P.gens()
    data=json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())
    f=P.fraction_field()(data['potential'])
    generators=[P(f.derivative(a).numerator()) for a in (x,u,v,y)]
    J=P.ideal(generators+[z*r-1])
    start=time.monotonic()
    print('full r-torus slice over QQ(t) started',flush=True)
    gb=J.groebner_basis()
    dim=int(J.dimension())
    print('slice basis',len(gb),'dimension',dim,'elapsed',time.monotonic()-start,flush=True)
    report=dict(status='incomplete',base='QQ(t)',central_r_retained=True,
                dimension=dim,groebner_generators=len(gb),target_equations_imposed=False)
    path=ROOT/'adjoint-central-boundary-slice.json'
    path.write_text(json.dumps(report,indent=2)+'\n')
    jac=matrix(P,[[g.derivative(a) for a in (x,u,v,y,r)] for g in generators])
    minors=[J.reduce(jac.matrix_from_columns(list(cols)).det()) for cols in combinations(range(5),4)]
    print('slice Jacobian minors reduced',flush=True)
    singular=J+P.ideal(minors)
    sd=int(singular.dimension())
    report.update(status='completed',singular_locus_dimension=sd,
                  slice_radical_by_complete_intersection_and_R0=(dim==1 and sd<1),
                  elapsed_seconds=time.monotonic()-start)
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


def diagonal_test():
    P=PolynomialRing(QQ,names='x,u,v,y,q,z');x,u,v,y,q,z=P.gens()
    F=P.fraction_field()
    data=json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())
    # Character restriction first is legitimate only for the AFFINE
    # derivatives. Central derivatives must still use the original character.
    text=data['potential'].replace('t','q').replace('r','q')
    f=F(text)
    J=P.ideal([P(f.derivative(a).numerator()) for a in (x,u,v,y)]+[z*q-1])
    start=time.monotonic()
    print('entire diagonal t=r source boundary started',flush=True)
    gb=J.groebner_basis()
    dim=int(J.dimension())
    report=dict(status='completed',boundary='t=r',dimension=dim,
                target_image_dimension_at_most=dim,
                cannot_dominate_two_dimensional_target=(dim<2),
                groebner_generators=len(gb),elapsed_seconds=time.monotonic()-start)
    (ROOT/'adjoint-central-boundary-diagonal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


def full_test():
    P=PolynomialRing(QQ,names='x,u,v,y,t,r,z');x,u,v,y,t,r,z=P.gens()
    data=json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())
    f=P.fraction_field()(data['potential'])
    generators=[P(f.derivative(a).numerator()) for a in (x,u,v,y)]
    J=P.ideal(generators+[z*t*r-1])
    start=time.monotonic()
    print('full relative ideal over QQ, both central coordinates retained',flush=True)
    gb=J.groebner_basis();dim=int(J.dimension())
    print('full basis',len(gb),'dimension',dim,'elapsed',time.monotonic()-start,flush=True)
    path=ROOT/'adjoint-central-boundary-full.json'
    report=dict(status='basis_only',dimension=dim,groebner_generators=len(gb),
                central_coordinates_retained=['t','r'])
    path.write_text(json.dumps(report,indent=2)+'\n')
    H=matrix(P,[[g.derivative(a) for a in (x,u,v,y)] for g in generators])
    det=J.reduce(H.det())
    print('full relative Hessian determinant reduced',flush=True)
    sd=int((J+P.ideal(det)).dimension())
    report.update(status='completed',relative_hessian_zero_locus_dimension=sd,
                  radical_by_complete_intersection_and_R0=(dim==2 and sd<2),
                  elapsed_seconds=time.monotonic()-start)
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


def simple_boundaries():
    data=json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())
    cases=[('t=r','q','q'),('t=-r','(-q)','q'),
           ('t*r=1','q','(1/q)'),('t*r=-1','q','(-1/q)'),
           ('t^2+1=0','ii','q'),('r^2+1=0','q','ii')]
    records=[]
    for label,tt,rr in cases:
        start=time.monotonic()
        C=QQ
        if 'ii' in (tt,rr):
            U=PolynomialRing(QQ,names='T');C=NumberField(U.gen()**2+1,'ii')
        P=PolynomialRing(C,names='x,u,v,y,q,z');x,u,v,y,q,z=P.gens()
        f=P.fraction_field()(data['potential'].replace('t',tt).replace('r',rr))
        J=P.ideal([P(f.derivative(a).numerator()) for a in (x,u,v,y)]+[z*q-1])
        gb=J.groebner_basis();dim=int(J.dimension())
        item=dict(boundary=label,dimension=dim,groebner_generators=len(gb),
                  cannot_dominate_two_dimensional_target=(dim<2),
                  elapsed_seconds=time.monotonic()-start)
        records.append(item);print(json.dumps(item),flush=True)
        (ROOT/'adjoint-simple-boundaries.json').write_text(json.dumps(
            dict(status='partial' if len(records)<len(cases) else 'completed',records=records),indent=2)+'\n')


def lift_test(reverse=False,tracked=False,order=None):
    from sage.libs.singular.function_factory import singular_function
    T,P,f=potential()
    if reverse:
        P=PolynomialRing(T.fraction_field(),names='y,v,u,x');f=P(str(f))
    if order:
        assert sorted(order.split(','))==['u','v','x','y']
        P=PolynomialRing(T.fraction_field(),names=order);f=P(str(f))
    J=P.ideal([f.derivative(x) for x in P.gens()])
    start=time.monotonic()
    if tracked:
        from sage.interfaces.singular import singular
        names=','.join(P.variable_names())
        singular.eval('ring boundary_ring=(0,t,r),('+names+'),dp;')
        singular.eval('short=0;')
        singular.eval('ideal boundary_ideal='+','.join(map(str,J.gens()))+';')
        singular.eval('matrix boundary_lift;')
        print('tracked standard basis started',flush=True)
        singular.eval('ideal boundary_basis=liftstd(boundary_ideal,boundary_lift);')
        size=int(singular.eval('size(boundary_basis)'))
        gb=[P(singular.eval('boundary_basis['+str(j+1)+']')) for j in range(size)]
        lifted=matrix(P,4,size,[P(singular.eval('boundary_lift['+str(i+1)+','+str(j+1)+']'))
                              for i in range(4) for j in range(size)])
        for j,g in enumerate(gb):
            leading=g.lc()
            gb[j]=g/leading
            for i in range(4):lifted[i,j]/=leading
    else:
        gb=J.groebner_basis();print('generic basis ready; lifting membership',flush=True)
        lifted=singular_function('lift')(J,P.ideal(gb))
    print('membership matrix returned',flush=True)
    assert lifted.nrows()==4 and lifted.ncols()==len(gb)
    for j,g in enumerate(gb):
        assert sum(J.gens()[i]*lifted[i,j] for i in range(4))==g
    assert all(g.reduce(gb)==0 for g in J.gens())
    print('every membership identity checked',flush=True)
    denominators={c.denominator() for f in list(gb)+list(lifted.list()) for c in f.coefficients()}
    factors=set()
    for d in denominators:
        for h,m in d.factor():factors.add(h/h.lc())
    report=dict(status='membership_certified',variable_order=list(P.variable_names()),
        denominator_factors=sorted(map(str,factors)),
        original_ideal_equals_monic_model_off_denominators=True,
        elapsed_seconds=time.monotonic()-start)
    tag='reverse' if reverse else 'forward'
    if order:tag=order.replace(',','')
    payload=dict(variable_order=list(P.variable_names()),generators=list(map(str,J.gens())),
                 basis=list(map(str,gb)),
                 membership=[[str(lifted[i,j]) for i in range(4)] for j in range(len(gb))])
    archive=ROOT/('adjoint-central-membership-'+tag+'.json.gz')
    archive.write_bytes(gzip.compress(json.dumps(payload,sort_keys=True).encode(),mtime=0))
    report.update(certificate_file=archive.name,certificate_sha256=sha256(archive.read_bytes()).hexdigest(),
                  both_ideal_inclusions_checked=True,monic_basis_generators=len(gb))
    (ROOT/('adjoint-central-boundary-lift-'+tag+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


def extra_boundary(which):
    data=json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())
    if which==1:
        U=PolynomialRing(QQ,names='T');C=NumberField(U.gen()**2+1,'ii')
    else:C=QQ
    P=PolynomialRing(C,names='x,u,v,y,t,r,z');x,u,v,y,t,r,z=P.gens()
    f=P.fraction_field()(data['potential'])
    if which==1:
        boundary=t*t*r+t*r*r-r-t-2*C.gen()*t*r
        label='t-t^-1+r-r^-1=2i'
    else:
        record=json.loads((ROOT/'adjoint-central-boundary-exploration.json').read_text())
        boundary=P(next(a for a in record['denominator_factors'] if '64/9' in a))
        label='degree_12_denominator'
    J=P.ideal([P(f.derivative(a).numerator()) for a in (x,u,v,y)]+[z*t*r-1,boundary])
    start=time.monotonic();print('extra boundary started',label,flush=True)
    gb=J.groebner_basis();dim=int(J.dimension())
    report=dict(status='completed',boundary=label,dimension=dim,groebner_generators=len(gb),
                cannot_dominate_two_dimensional_target=(dim<2),
                elapsed_seconds=time.monotonic()-start)
    (ROOT/('adjoint-extra-boundary-'+str(which)+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


def square_full():
    from sage.all import prod
    T,A,f=potential();t,r=map(T.fraction_field(),T.gens());x,u,v,y=A.gens()
    transformed=f.subs({x:t*x,y:r*y,u:r*u,v:r*t*t*v})
    P=PolynomialRing(QQ,names='x,u,v,y,p,q,z');x,u,v,y,p,q,z=P.gens()
    phi=P.zero()
    for ex,c in transformed.dict().items():
        numerator=T(c)
        for (a,b),value in numerator.dict().items():
            assert a%2==b%2==0
            phi+=value*p**(a//2)*q**(b//2)*prod(xx**ee for xx,ee in zip((x,u,v,y),ex))
    generators=[phi.derivative(a) for a in (x,u,v,y)]
    J=P.ideal(generators+[z*p*q-1])
    start=time.monotonic();print('full polynomial square-coordinate ideal started',flush=True)
    gb=J.groebner_basis();dim=int(J.dimension())
    print('square-coordinate basis',len(gb),'dimension',dim,flush=True)
    path=ROOT/'adjoint-central-boundary-square-full.json'
    report=dict(status='basis_only',potential=str(phi),dimension=dim,
                coordinate_change='old x=t*x, old y=r*y, old u=r*u, old v=r*t^2*v; p=t^2,q=r^2')
    path.write_text(json.dumps(report,indent=2)+'\n')
    H=matrix(P,[[g.derivative(a) for a in (x,u,v,y)] for g in generators])
    det=J.reduce(H.det());print('square Hessian reduced',flush=True)
    sd=int((J+P.ideal(det)).dimension())
    report.update(status='completed',hessian_zero_dimension=sd,
                  full_relative_ideal_radical=(dim==2 and sd<2),
                  elapsed_seconds=time.monotonic()-start)
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    import resource
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    parser=argparse.ArgumentParser()
    parser.add_argument('--slice',action='store_true')
    parser.add_argument('--diagonal',action='store_true')
    parser.add_argument('--full',action='store_true')
    parser.add_argument('--simple-boundaries',action='store_true')
    parser.add_argument('--lift',action='store_true')
    parser.add_argument('--reverse',action='store_true')
    parser.add_argument('--tracked',action='store_true')
    parser.add_argument('--order')
    parser.add_argument('--extra-boundary',type=int,choices=(1,2))
    parser.add_argument('--square-full',action='store_true')
    args=parser.parse_args()
    if args.slice:
        slice_test()
    elif args.diagonal:
        diagonal_test()
    elif args.full:
        full_test()
    elif args.simple_boundaries:
        simple_boundaries()
    elif args.lift:
        lift_test(args.reverse,args.tracked,args.order)
    elif args.extra_boundary:
        extra_boundary(args.extra_boundary)
    elif args.square_full:
        square_full()
    else:
        main()
