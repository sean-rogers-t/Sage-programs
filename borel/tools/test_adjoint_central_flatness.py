"""Test the proposed central flatness route on the remaining E8 candidate.

Source specializations do not by themselves make generic-target assertions.
"""
from pathlib import Path
from hashlib import sha256
import json
import time
from sage.all import QQ, PolynomialRing, NumberField

ROOT = Path(__file__).resolve().parents[1]


def main():
    data = json.loads((ROOT / 'adjoint-nilpotent-source-screen.json').read_text())
    P = PolynomialRing(QQ, names='x,u,v,y')
    records = []
    for tt, rr in [(1, 2), (-1, 2), (2, 1), (2, -1), (2, 3)]:
        start = time.monotonic()
        b = (tt**2-rr**2)*((tt*rr)**2-1)*(tt**2+1)*(rr**2+1)
        assert b != 0
        phi = P(data['potential'].replace('t', '('+str(tt)+')').replace('r', '('+str(rr)+')'))
        J = P.ideal([phi.derivative(a) for a in P.gens()])
        gb = J.groebner_basis()
        dim = int(J.dimension())
        item = dict(t=tt, r=rr, b=b, dimension=dim,
                    length=int(J.vector_space_dimension()) if dim == 0 else None,
                    groebner_generators=len(gb),
                    basis=list(map(str, gb)),
                    elapsed_seconds=time.monotonic()-start)
        print(json.dumps({k:v for k,v in item.items() if k != 'basis'}), flush=True)
        records.append(item)
        (ROOT / 'adjoint-central-flatness-screen.json').write_text(json.dumps(
            dict(status='source_specializations_only', records=records,
                 source_sha256=source_hashes()), indent=2)+'\n')


def shared_curves():
    data = json.loads((ROOT / 'adjoint-nilpotent-source-screen.json').read_text())
    forward = json.loads((ROOT / 'adjoint-central-boundary-lift-forward.json').read_text())
    reverse = json.loads((ROOT / 'adjoint-central-boundary-lift-reverse.json').read_text())
    T = PolynomialRing(QQ, names='t,r'); t, r = T.gens()
    b = (t*t-r*r)*((t*r)**2-1)*(t*t+1)*(r*r+1)*t*r
    common = [T(s) for s in forward['denominator_factors']
              if s in reverse['denominator_factors'] and b % T(s) != 0]
    U = PolynomialRing(QQ, 'q'); q = U.gen()
    records = []
    for h in common:
        for hp, _ in U(h(t=q, r=2)).factor():
            C = NumberField(hp, 'a'); tt = C.gen(); rr = C(2)
            bb = (tt*tt-rr*rr)*((tt*rr)**2-1)*(tt*tt+1)*(rr*rr+1)
            assert bb != 0
            P = PolynomialRing(C, names='x,u,v,y')
            phi = P(data['potential'].replace('t', '(a)').replace('r', '(2)'))
            J = P.ideal([phi.derivative(a) for a in P.gens()])
            start = time.monotonic()
            print('special curve started', h, hp, flush=True)
            gb = J.groebner_basis(); dim = int(J.dimension())
            item = dict(curve=str(h), r=2, t_minimal_polynomial=str(hp),
                        dimension=dim, length=int(J.vector_space_dimension()) if dim == 0 else None,
                        groebner_generators=len(gb), basis=list(map(str, gb)),
                        elapsed_seconds=time.monotonic()-start)
            records.append(item); print(json.dumps({k:v for k,v in item.items() if k != 'basis'}), flush=True)
            (ROOT / 'adjoint-central-flatness-curve-screen.json').write_text(json.dumps(
                dict(status='source_specializations_only', records=records,
                     source_sha256=source_hashes()), indent=2)+'\n')


def source_hashes():
    return {p.name:sha256(p.read_bytes()).hexdigest() for p in
            (Path(__file__), ROOT/'adjoint-nilpotent-source-screen.json')}


def whole_basis():
    from sage.interfaces.singular import singular
    data = json.loads((ROOT / 'adjoint-nilpotent-source-screen.json').read_text())
    P = PolynomialRing(QQ, names='x,u,v,y,t,r,z')
    x,u,v,y,t,r,z = P.gens()
    phi = P.fraction_field()(data['potential'])
    generators = [P(phi.derivative(a).numerator()) for a in (x,u,v,y)] + [z*t*r-1]
    singular.eval('ring full_ring=0,(x,u,v,y,t,r,z),dp; short=0;')
    singular.eval('ideal full_ideal='+','.join(map(str,generators))+';')
    start = time.monotonic()
    print('whole original ideal slimgb started', flush=True)
    singular.eval('ideal full_basis=slimgb(full_ideal);')
    size = int(singular.eval('size(full_basis)'))
    dim = int(singular.eval('dim(full_basis)'))
    report = dict(status='whole_basis_only', dimension=dim, groebner_generators=size,
                  elapsed_seconds=time.monotonic()-start,
                  basis=[singular.eval('full_basis['+str(i+1)+']') for i in range(size)])
    (ROOT/'adjoint-central-flatness-whole-basis.json').write_text(json.dumps(report, indent=2)+'\n')
    print('whole basis', size, dim, time.monotonic()-start, flush=True)


def comprehensive():
    from sage.interfaces.singular import singular
    data = json.loads((ROOT / 'adjoint-nilpotent-source-screen.json').read_text())
    P = PolynomialRing(QQ, names='x,u,v,y,t,r')
    x,u,v,y,t,r = P.gens(); phi = P.fraction_field()(data['potential'])
    generators = [P(phi.derivative(a).numerator()) for a in (x,u,v,y)]
    singular.eval('LIB "grobcov.lib"; ring parametric_ring=(0,t,r),(x,u,v,y),dp; short=0;')
    singular.eval('ideal parametric_ideal='+','.join(map(str,generators))+';')
    singular.eval('ideal boundary=t*r*(t^2-r^2)*((t*r)^2-1)*(t^2+1)*(r^2+1);')
    start = time.monotonic(); print('comprehensive basis over b open started', flush=True)
    singular.eval('list cover=cgsdr(parametric_ideal,"nonnull",boundary,"can",2);')
    size = int(singular.eval('size(cover)'))
    (ROOT/'adjoint-central-flatness-comprehensive.txt').write_text(singular.eval('cover;')+'\n')
    report = dict(status='cover_computed_needs_audit', segments=size, elapsed_seconds=time.monotonic()-start)
    (ROOT/'adjoint-central-flatness-comprehensive.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


def cover_intersections():
    from sage.all import prod
    data = [json.loads((ROOT/('adjoint-central-boundary-lift-'+tag+'.json')).read_text())
            for tag in ('forward','reverse')]
    P = PolynomialRing(QQ, names='t,r'); t,r = P.gens()
    common = set(data[0]['denominator_factors']) & set(data[1]['denominator_factors'])
    common_product = prod(P(h) for h in common)
    extras = [[P(h) for h in d['denominator_factors'] if h not in common] for d in data]
    records = []
    for i, f in enumerate(extras[0]):
        for j, g in enumerate(extras[1]):
            J = P.ideal([f,g])
            if J.is_one(): continue
            K, exponent = J.saturation(P.ideal(common_product))
            if K.is_one(): continue
            item = dict(forward=str(f), reverse=str(g), dimension=int(K.dimension()),
                        length=int(K.vector_space_dimension()), basis=list(map(str,K.groebner_basis())))
            records.append(item)
            print('uncovered intersection',i,j,item['length'],item['basis'],flush=True)
    report=dict(status='base_cover_intersections_only', records=records,
                common_factors=sorted(common), residual_pairs=len(records))
    (ROOT/'adjoint-central-flatness-cover-intersections.json').write_text(json.dumps(report,indent=2)+'\n')


def generic_curve(index):
    data = json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())
    cover = json.loads((ROOT/'adjoint-central-flatness-cover-intersections.json').read_text())
    U = PolynomialRing(QQ, 'r'); K = U.fraction_field()
    P = PolynomialRing(K, names='x,u,v,y,t'); x,u,v,y,t=P.gens(); r=P(K.gen())
    b = (t*t-r*r)*((t*r)**2-1)*(t*t+1)*(r*r+1)*t*r
    curves = [P(s) for s in cover['common_factors'] if b % P(s) != 0]
    h = curves[index]
    phi = P.fraction_field()(data['potential'])
    generators = [P(phi.derivative(a).numerator()) for a in (x,u,v,y)] + [h]
    # The nonzero constant term of h in t makes t invertible here.
    assert h(t=0) in K and h(t=0) != 0
    J = P.ideal(generators)
    print('generic curve basis started',h,flush=True); start=time.monotonic()
    gb=J.groebner_basis(); dim=int(J.dimension())
    report=dict(status='generic_curve_only',curve=str(h),dimension=dim,
                length=int(J.vector_space_dimension()) if dim==0 else None,
                basis=list(map(str,gb)),elapsed_seconds=time.monotonic()-start)
    (ROOT/('adjoint-central-flatness-generic-curve-'+str(index)+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    print('generic curve result',index,dim,report['length'],report['elapsed_seconds'],flush=True)


if __name__ == '__main__':
    import sys
    import resource
    resource.setrlimit(resource.RLIMIT_AS, (4*1024**3, 4*1024**3))
    if '--whole-basis' in sys.argv: whole_basis()
    elif '--generic-curve' in sys.argv: generic_curve(int(sys.argv[-1]))
    elif '--cover-intersections' in sys.argv: cover_intersections()
    elif '--comprehensive' in sys.argv: comprehensive()
    elif '--shared-curves' in sys.argv: shared_curves()
    else: main()
