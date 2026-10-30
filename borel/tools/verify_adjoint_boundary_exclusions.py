"""Exact dimension-upper-bound certificates for six selected source boundaries.

Every witness is lifted to the original affine-gradient ideal (with the
torus inverse). A monomial in each coordinate two-plane bounds dimension
by one, regardless of whether the selected witnesses form a Groebner basis.
"""
from pathlib import Path
from itertools import combinations
from hashlib import sha256
import gzip
import json
import sys
import re
from functools import lru_cache
from sage.all import QQ, PolynomialRing, NumberField, sage_eval
from sage.interfaces.singular import singular

ROOT=Path(__file__).resolve().parents[1]
CASES=[('t=r','q','q'),('t=-r','(-q)','q'),
       ('t*r=1','q','(1/q)'),('t*r=-1','q','(-1/q)'),
       ('t^2+1=0','ii','q'),('r^2+1=0','q','ii')]
TRANSFERS={'t*r=-1':('t*r=1','sign_t'),
           'r^2+1=0':('t^2+1=0','swap')}


def transfer_record(label,source):
    reference,kind=TRANSFERS[label]
    P=PolynomialRing(QQ,names='x,u,v,y,t,r');x,u,v,y,t,r=P.gens()
    F=P.fraction_field();phi=F(source)
    values=(-x,u,v,y,-t,r) if kind=='sign_t' else (y,t*u/r,r*v/t,x,r,t)
    substitution=dict(zip(P.gens(),values))
    assert phi.subs(substitution)==phi
    assert all(F(v).subs(substitution)==g for v,g in zip(values,P.gens()))
    # The transformations are Laurent polynomial automorphisms, affine-linear
    # and invertible in x,u,v,y over the central torus. Chain rule preserves
    # the entire original affine-gradient scheme, including nilpotents.
    from sage.all import matrix
    assert matrix(F,[[F(value).derivative(g) for g in (x,u,v,y)] for value in values[:4]]).det() in (1,-1)
    if kind=='sign_t':assert values[4]*values[5]-1==-(t*r+1)
    else:assert values[4]**2+1==r*r+1
    return dict(boundary=label,reference_boundary=reference,certificate_kind='potential_involution',
                involution=kind,substitution=list(map(str,values)),
                dimension_upper_bound=1,target_image_dimension_upper_bound=1)


def sparse_parser(P):
    """Parse expanded Singular polynomials without a giant symbolic AST."""
    C=P.base_ring(); pattern=re.compile(r'\*?([xuvyqz])(?:\^(\d+))?')
    positions={v:i for i,v in enumerate('xuvyqz')}
    @lru_cache(maxsize=20000)
    def coefficient(txt):
        if txt in ('','+'):return C.one()
        if txt=='-':return -C.one()
        return C(sage_eval(txt,locals={'ii':C.gen()} if C!=QQ else {}))
    def parse(txt):
        txt=txt.replace(' ','').replace('\n','');parts=[];start=0;depth=0
        for i,c in enumerate(txt):
            if c=='(':depth+=1
            elif c==')':depth-=1
            elif c in '+-' and depth==0 and i>start:
                parts.append(txt[start:i]);start=i
        parts.append(txt[start:]);answer={}
        for part in parts:
            exponent=[0]*6
            for match in pattern.finditer(part):exponent[positions[match[1]]]+=int(match[2] or 1)
            exponent=tuple(exponent);cc=coefficient(pattern.sub('',part))
            answer[exponent]=answer.get(exponent,C.zero())+cc
        return P(answer)
    return parse


def main():
    sys.setrecursionlimit(30000)
    data=json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())
    archive=ROOT/'adjoint-boundary-dimension-certificates.json.gz'
    records=json.loads(gzip.decompress(archive.read_bytes())) if archive.exists() else []
    # Replace only newly generated redundant lifts by exact small isomorphism
    # certificates; the original two saved membership proofs are preserved.
    records=[transfer_record(record['boundary'],data['potential'])
             if record['boundary'] in TRANSFERS else record for record in records]
    for label,tt,rr in CASES:
        if any(record['boundary']==label for record in records):
            print(label,'reusing saved exact certificate (independent audit replays it)',flush=True)
            continue
        if label in TRANSFERS:
            assert any(record['boundary']==TRANSFERS[label][0] for record in records)
            records.append(transfer_record(label,data['potential']))
            print(label,'exact original-potential involution checked',flush=True)
            continue
        C=QQ
        if 'ii' in (tt,rr):
            U=PolynomialRing(QQ,names='T');C=NumberField(U.gen()**2+1,'ii')
        P=PolynomialRing(C,names='x,u,v,y,q,z');x,u,v,y,q,z=P.gens()
        parse=sparse_parser(P)
        f=P.fraction_field()(data['potential'].replace('t',tt).replace('r',rr))
        generators=[P(f.derivative(a).numerator()) for a in (x,u,v,y)]+[z*q-1]
        coefficient='0' if C==QQ else '(0,ii)'
        singular.eval('ring certificate_ring='+coefficient+',(x,u,v,y,q,z),dp;')
        if C!=QQ:singular.eval('minpoly=ii^2+1;')
        singular.eval('short=0;')
        singular.eval('ideal certificate_ideal='+','.join(map(str,generators))+';')
        singular.eval('matrix certificate_lift;')
        singular.eval('ideal certificate_basis=liftstd(certificate_ideal,certificate_lift);')
        size=int(singular.eval('size(certificate_basis)'))
        gb=[parse(singular.eval('certificate_basis['+str(j+1)+']')) for j in range(size)]
        witnesses=[];costs={}
        def lift_size(g):
            j=gb.index(g)+1
            if j not in costs:
                costs[j]=sum(int(singular.eval('size(certificate_lift['+str(i+1)+','+str(j)+'])')) for i in range(5))
            return costs[j]
        for pair in combinations(range(6),2):
            possible=[g for g in gb if set(i for i,e in enumerate(g.lm().exponents()[0]) if e)<=set(pair)]
            assert possible,label
            selected=min(possible,key=lift_size)
            if selected not in witnesses:witnesses.append(selected)
        print(label,'selected',len(witnesses),'membership witnesses',flush=True)
        columns=gb
        lift=[[parse(singular.eval('certificate_lift['+str(i+1)+','+str(columns.index(g)+1)+']'))
               for g in witnesses] for i in range(5)]
        for j,g in enumerate(witnesses):
            assert sum(generators[i]*lift[i][j] for i in range(5))==g
        records.append(dict(boundary=label,t_substitution=tt,r_substitution=rr,
            field='QQ(i)' if C!=QQ else 'QQ',variables=list(P.variable_names()),
            generators=list(map(str,generators)),
            witnesses=[dict(polynomial=str(g),leading_exponent=list(map(int,g.lm().exponents()[0])),
                            membership=[str(lift[i][j]) for i in range(5)]) for j,g in enumerate(witnesses)],
            dimension_upper_bound=1,target_image_dimension_upper_bound=1))
        print(label,'exact memberships checked',flush=True)
        (ROOT/'adjoint-boundary-dimension-certificates.json.gz').write_bytes(
            gzip.compress(json.dumps(records,sort_keys=True).encode(),mtime=0))
    archive=ROOT/'adjoint-boundary-dimension-certificates.json.gz'
    archive.write_bytes(gzip.compress(json.dumps(records,sort_keys=True).encode(),mtime=0))
    report=dict(result='passed',boundary_count=len(records),
                direct_membership_boundaries=sum('witnesses' in r for r in records),
                potential_involution_boundaries=sum('involution' in r for r in records),
                every_dimension_upper_bound=1,all_target_images_nondominant=True,
                full_remaining_boundary_resolved=False,uniform_conjecture_proved=False,
                certificate_file=archive.name,
                source_sha256={p.name:sha256(p.read_bytes()).hexdigest() for p in
                    (Path(__file__),ROOT/'adjoint-nilpotent-source-screen.json',archive)})
    (ROOT/'adjoint-boundary-exclusions-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    import resource
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    main()
