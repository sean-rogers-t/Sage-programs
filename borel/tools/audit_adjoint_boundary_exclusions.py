"""Sage-free verification of explicit original-ideal dimension certificates."""
from pathlib import Path
from itertools import combinations
from hashlib import sha256
import gzip
import json
import sympy as s
import re
from functools import lru_cache
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]


def main(require_all=True):
    archive=ROOT/'adjoint-boundary-dimension-certificates.json.gz'
    records=json.loads(gzip.decompress(archive.read_bytes()))
    source=json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())
    variables=s.symbols('x u v y q z');x,u,v,y,q,z=variables
    t,r=s.symbols('t r')
    namespace=dict(zip(('x','u','v','y','q','z'),variables),t=t,r=r,ii=s.I)
    parse=lambda a:s.sympify(a.replace('^','**'),locals=namespace)
    phi=parse(source['potential'])
    sparse_domain=s.QQ_I
    sparse_ring,_x,_u,_v,_y,_q,_z=ring('x,u,v,y,q,z',sparse_domain)
    pattern=re.compile(r'\*?([xuvyqz])(?:\^(\d+))?')
    positions={v:i for i,v in enumerate('xuvyqz')}
    @lru_cache(maxsize=20000)
    def coefficient(txt):
        if txt in ('','+'):txt='1'
        elif txt=='-':txt='-1'
        if sparse_domain==s.QQ:
            numerator,separator,denominator=txt.strip('()').partition('/')
            return s.QQ(int(numerator),int(denominator) if separator else 1)
        return sparse_domain.from_sympy(s.sympify(txt,locals={'ii':s.I}))
    def sparse(txt):
        txt=txt.replace(' ','');parts=[];start=0;depth=0
        for i,c in enumerate(txt):
            if c=='(':depth+=1
            elif c==')':depth-=1
            elif c in '+-' and depth==0 and i>start:
                parts.append(txt[start:i]);start=i
        parts.append(txt[start:]);answer={}
        for part in parts:
            exponent=[0]*6
            for match in pattern.finditer(part):exponent[positions[match[1]]]+=int(match[2] or 1)
            cc=coefficient(pattern.sub('',part));exponent=tuple(exponent)
            answer[exponent]=answer.get(exponent,sparse_domain.zero)+cc
        return sparse_ring.from_dict(answer)
    checked=[]
    for record in records:
        if record.get('certificate_kind')=='potential_involution':
            reference=record['reference_boundary']
            assert any(row['boundary']==reference for row in checked)
            original_vars=(x,u,v,y,t,r)
            values=tuple(parse(a) for a in record['substitution'])
            substitution=dict(zip(original_vars,values))
            assert s.cancel(phi.subs(substitution,simultaneous=True)-phi)==0
            assert all(s.cancel(a.subs(substitution,simultaneous=True)-b)==0
                       for a,b in zip(values,original_vars))
            # All coordinate functions and the inverse are Laurent regular on
            # the central torus; no affine or central divisor is removed.
            for value in values:
                numerator,denominator=s.fraction(s.cancel(value))
                assert not (denominator.free_symbols-{t,r})
                assert len(s.Poly(denominator,t,r).terms())==1
            jac=s.Matrix(values[:4]).jacobian((x,u,v,y))
            assert s.cancel(jac.det()) in (1,-1)
            transformed_grad=s.Matrix([s.diff(phi,a).subs(substitution,simultaneous=True)
                                       for a in (x,u,v,y)])
            assert all(s.cancel(a)==0 for a in jac.T*transformed_grad-s.Matrix([s.diff(phi,a) for a in (x,u,v,y)]))
            if record['boundary']=='t*r=-1':
                assert reference=='t*r=1' and s.cancel(values[4]*values[5]-1+t*r+1)==0
            else:
                assert record['boundary']=='r^2+1=0' and reference=='t^2+1=0'
                assert s.cancel(values[4]**2-r*r)==0
            checked.append(dict(boundary=record['boundary'],method='original-potential involution and full gradient chain rule',
                                reference_boundary=reference,dimension_upper_bound=1,
                                target_image_dimension_upper_bound=1))
            print(record['boundary'],'independent original-scheme transfer verified',flush=True)
            continue
        sparse_domain=s.QQ if record['field']=='QQ' else s.QQ_I
        sparse_ring,*_=ring('x,u,v,y,q,z',sparse_domain)
        coefficient.cache_clear()
        def poly(f):return s.Poly(f,*variables,extension=s.I)
        f=phi.subs({t:parse(record['t_substitution']),r:parse(record['r_substitution'])},simultaneous=True)
        generators=[poly(parse(a)) for a in record['generators']]
        sparse_generators=[sparse(a) for a in record['generators']]
        integer_ring=sparse_ring.clone(domain=s.ZZ) if sparse_domain==s.QQ else None
        if integer_ring is not None:
            integer_generators=[(d,a.set_ring(integer_ring)) for d,a in
                                (g.clear_denoms() for g in sparse_generators)]
        for i,a in enumerate((x,u,v,y)):
            ratio=s.cancel(generators[i].as_expr()/s.diff(f,a))
            assert not (ratio.free_symbols-set([q]))
            num,den=s.fraction(ratio)
            assert len(s.Poly(num,q,extension=s.I).terms())==1
            assert len(s.Poly(den,q,extension=s.I).terms())==1
        assert generators[4]==poly(z*q-1)
        supports=[]
        for witness in record['witnesses']:
            g=sparse(witness['polynomial'])
            coeffs=[sparse(a) for a in witness['membership']]
            if integer_ring is not None:
                dg,gg=g.clear_denoms();gg=gg.set_ring(integer_ring)
                lifted=[(d,a.set_ring(integer_ring)) for d,a in (f.clear_denoms() for f in coeffs)]
                common=dg
                for (da,a),(db,b) in zip(integer_generators,lifted):
                    common=s.ZZ.lcm(common,da*db)
                value=sum((a*b*(common//(da*db)) for (da,a),(db,b) in zip(integer_generators,lifted)),integer_ring.zero)
                assert value==gg*(common//dg),record['boundary']
            else:
                value=sum((a*b for a,b in zip(sparse_generators,coeffs)),sparse_ring.zero)
                assert value==g,record['boundary']
            ex=max(g.keys(),key=lambda e:(sum(e),tuple(-a for a in reversed(e))))
            assert list(ex)==witness['leading_exponent']
            supports.append({i for i,e in enumerate(ex) if e})
        assert all(any(support<=set(pair) for support in supports)
                   for pair in combinations(range(6),2))
        checked.append(dict(boundary=record['boundary'],membership_identities=len(supports),
                            dimension_upper_bound=1,target_image_dimension_upper_bound=1))
        print(record['boundary'],'independently verified',flush=True)
    if require_all:assert len(checked)==6
    report=dict(result='passed',uses_Sage=False,boundaries=checked,
                independently_replayed_count=len(checked),partial_audit=(len(checked)<6),
                direct_membership_boundaries=sum('membership_identities' in row for row in checked),
                potential_involution_boundaries=sum('reference_boundary' in row for row in checked),
                original_gradients_reconstructed=True,
                rational_memberships_checked_after_exact_denominator_clearing=True,
                every_membership_identity_verified=True,
                dimension_bounds_do_not_assume_witnesses_are_groebner_basis=True,
                target_nondominance_proved_for_six_boundaries=(len(checked)==6),
                entire_remaining_boundary_resolved=False,
                source_sha256={p.name:sha256(p.read_bytes()).hexdigest() for p in
                    (Path(__file__),archive,ROOT/'adjoint-nilpotent-source-screen.json')})
    (ROOT/'adjoint-boundary-exclusions-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    import sys
    main(require_all='--partial' not in sys.argv)
