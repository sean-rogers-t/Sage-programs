"""Independent modular corroboration of exact characteristic-zero fiber lengths.

This does not certify flatness on a whole parameter open.
"""
from pathlib import Path
from hashlib import sha256
from itertools import product
import json
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def main():
    variables = s.symbols('x u v y'); t,r,q,a=s.symbols('t r q a')
    namespace = dict(zip('xuvy', variables), t=t, r=r, q=q, a=a)
    parse = lambda f: s.sympify(f.replace('^','**'), locals=namespace)
    source = ROOT/'adjoint-nilpotent-source-screen.json'
    screen = ROOT/'adjoint-central-flatness-curve-screen.json'
    phi = parse(json.loads(source.read_text())['potential'])
    records = json.loads(screen.read_text())['records']; checked=[]
    for item in records:
        hp = parse(item['t_minimal_polynomial'])
        assert s.Poly(hp,q,domain=s.QQ).is_irreducible
        b_at_r2=(q*q-4)*(4*q*q-1)*(q*q+1)*5
        assert s.degree(s.gcd(hp,b_at_r2),q)==0
        assert s.rem(parse(item['curve']).subs({r:2,t:q}),hp,q)==0
        original_basis=[s.Poly(parse(g),*variables) for g in item['basis']]
        assert all(g.LC(order='grevlex')==1 for g in original_basis)
        prime=65521
        def mod(f, vv=variables):
            return sum((int(c.p)*pow(int(c.q),-1,prime)%prime)*s.prod(v**e for v,e in zip(vv,ex))
                       for ex,c in s.Poly(f,*vv,domain=s.QQ).terms())
        roots=[]
        while not roots:
            factors=s.factor_list(mod(hp,(q,)),q,modulus=prime)[1]
            for f,m in factors:
                fp=s.Poly(f,q,modulus=prime)
                if fp.degree()==1:
                    aa,bb=map(int,fp.all_coeffs()); roots.append((-bb*pow(aa,-1,prime))%prime)
            if not roots: prime=int(s.nextprime(prime))
        tt=roots[0]; rr=2
        assert ((tt*tt-rr*rr)*((tt*rr)**2-1)*(tt*tt+1)*(rr*rr+1)*tt*rr)%prime
        f=phi.subs({t:tt,r:rr})
        partials=[mod(s.diff(f,v)) for v in variables]
        basis=s.groebner(partials,*variables,modulus=prime,order='grevlex')
        stored=[s.Poly(mod(parse(g).subs(a,tt)),*variables,modulus=prime) for g in item['basis']]
        assert set(basis.polys)==set(stored)
        assert [g.LM(order='grevlex').exponents for g in stored]==[
            g.LM(order='grevlex').exponents for g in original_basis]
        leading=[g.LM(order='grevlex').exponents for g in basis.polys]
        bounds=[min(e[i] for e in leading if e[i] and sum(e)==e[i]) for i in range(4)]
        standard=[e for e in product(*(range(n) for n in bounds))
                  if not any(all(a<=b for a,b in zip(lm,e)) for lm in leading)]
        assert len(standard)==item['length']
        checked.append(dict(prime=prime,t=tt,r=rr,length=len(standard),
                            curve=item['curve'],basis_recomputed_and_equal=True))
        print('independently corroborated length',len(standard),'at',tt,rr,'mod',prime,flush=True)
    report=dict(result='passed',independent_engine='SymPy',checks=checked,
                generic_relative_rank=38,finite_over_entire_b_open_ruled_out=True,
                characteristic_zero_claim_source='Sage number-field Groebner bases',
                independent_check_scope='Modular recomputation and exact basis comparison',
                whole_open_flatness_proved=False,generic_target_reducedness_proved=False,
                source_sha256={p.name:sha256(p.read_bytes()).hexdigest() for p in
                               (Path(__file__),source,screen)})
    (ROOT/'adjoint-central-flatness-audit.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__': main()
