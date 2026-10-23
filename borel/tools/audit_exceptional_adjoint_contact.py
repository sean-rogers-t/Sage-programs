"""Sage-free root and module reconstruction for the contact-fiber proof.

Shares Bourbaki numbering, not Sage root or character implementations.
The explicit Cartan matrices and reflection routines come from the earlier
native certificate auditor; expected modules are built here combinatorially.
"""
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json

import sympy as sp
from audit_affine_branching_certificates import cartan, root_data, inverse, transpose

ROOT=Path(__file__).resolve().parents[1]


def local_cartan(kind):
    if kind.startswith('E'):return cartan(kind)
    n=int(kind[1:]);A=[[2*int(i==j) for j in range(n)] for i in range(n)]
    if kind[0] in ('A','C'):
        for i in range(n-1):A[i][i+1]=A[i+1][i]=-1
        if kind[0]=='C':A[n-2][n-1]=-2
    elif kind[0]=='D':
        for i in range(n-2):A[i][i+1]=A[i+1][i]=-1
        A[n-3][n-1]=A[n-1][n-3]=-1
    else:raise ValueError(kind)
    return tuple(tuple(row) for row in A)


def all_roots(A):
    positive=root_data(A)[1]
    return list(positive)+[tuple(-v for v in w) for w in positive]


def orbit(A,highest):
    found={highest};todo=deque([highest])
    while todo:
        w=todo.popleft()
        for i in range(len(A)):
            v=tuple(w[j]-w[i]*A[j][i] for j in range(len(A)))
            if v not in found:found.add(v);todo.append(v)
    return Counter({w:1 for w in found})


def model_weights(kind,highest):
    n=int(kind[1:])
    if kind=='A1':return Counter({(3,):1,(1,):1,(-1,):1,(-3,):1})
    if kind=='A5':
        out=Counter()
        for selected in combinations(range(6),3):
            a=[int(i in selected) for i in range(6)]
            out[tuple(a[i]-a[i+1] for i in range(5))]+=1
        return out
    if kind=='C3':
        vectors=[tuple(sign*int(i==j) for i in range(3)) for j in range(3) for sign in (-1,1)]
        labels=lambda a:(a[0]-a[1],a[1]-a[2],a[2])
        out=Counter(labels(tuple(sum(v[i] for v in subset) for i in range(3)))
                    for subset in combinations(vectors,3))
        out.subtract(Counter(labels(v) for v in vectors))
        assert all(m>=0 for m in out.values())
        return +out
    if kind=='D6':
        parity=0 if highest[-1]==1 else 1
        out=Counter()
        for signs in product((-1,1),repeat=6):
            if signs.count(-1)%2!=parity:continue
            a=[F(v,2) for v in signs]
            w=tuple(int(a[i]-a[i+1]) for i in range(5))+(int(a[4]+a[5]),)
            out[w]+=1
        return out
    assert kind=='E7' and highest==(0,0,0,0,0,0,1)
    # This fundamental representation is minuscule: its weights are its orbit.
    return orbit(cartan('E7'),highest)


def reconstruct(record):
    kind=record['type'];A=cartan(kind);n=len(A);j=record['omitted_node']-1
    order=[i-1 for i in record['derived_standard_nodes']]
    local=local_cartan(record['derived_type'])
    assert tuple(tuple(A[i][k] for k in order) for i in order)==local
    H=tuple(row[j] for row in inverse(transpose(A)))
    assert list(H)==record['grading_cocharacter'] and H[j]==2
    roots=all_roots(A);grades={k:Counter() for k in (-2,-1,0,1,2)}
    for w in roots:
        grade=sum(w[i]*H[i] for i in range(n));assert grade.denominator==1
        grades[int(grade)][tuple(w[i] for i in order)]+=1
    grades[0][(0,)*(n-1)]+=n
    for k,weights in grades.items():
        reported=Counter({tuple(r['weight']):r['multiplicity'] for r in record['grades'][str(k)]})
        assert weights==reported
    expected=model_weights(record['derived_type'],tuple(record['V_highest']))
    assert grades[1]==grades[-1]==expected
    adj=Counter(all_roots(local));adj[(0,)*(n-1)]+=n
    assert grades[0]==adj
    assert grades[2]==grades[-2]==Counter({(0,)*(n-1):1})
    trace=sum(sum(w[i]*H[i] for i in range(n))**2 for w in roots)
    assert trace==2*sum(expected.values())+8==record['trace_normalization']
    if kind=='E6':
        c=[int(H[i])%2 for i in order]
        diagonal=[c[0]]+[c[i]-c[i-1] for i in range(1,5)]+[-c[-1]]
        assert all(v%2==1 for v in diagonal)
    return dict(type=kind,roots=len(roots),V_dimension=sum(expected.values()),
                all_graded_weights_reconstructed=True,trace_normalization=int(trace),
                cocharacter_cover_index=int(H[j]))


def native_e6_point():
    roots=all_roots(cartan('E6'));H=[1,2,2,3,2,1];order=[0,2,3,4,5]
    # Monomials 2^a * omega^b * i^c, omega a primitive cube root.
    eig=[(1,0),(1,1),(1,2),(-1,0),(-1,1),(-1,2)]
    params=[[0,0,h] for h in H]
    for k,pos in enumerate(order):
        params[pos][0]=sum(a for a,b in eig[:k+1])
        params[pos][1]=sum(b for a,b in eig[:k+1])
    omega=(-1+sp.sqrt(3)*sp.I)/2
    phases={(b,c):sp.expand(omega**b*sp.I**c) for b in range(3) for c in range(4)}
    values=[]
    for w in roots:
        a,b,c=[sum(w[i]*params[i][j] for i in range(6)) for j in range(3)]
        value=sp.expand(sp.Rational(2)**a*phases[b%3,c%4])
        assert value!=1;values.append(value)
    B=sp.Matrix(6,6,lambda i,j:sum(w[i]*w[j] for w in roots))
    grad=sp.Matrix([sp.expand(sum(value*w[i] for value,w in zip(values,roots))) for i in range(6)])
    target=sp.Rational(65,192)*sp.I*sp.Matrix(H)
    assert all(sp.expand(e)==0 for e in grad-B*target)
    Hess=sp.Matrix(6,6,lambda i,j:sp.expand(sum(value*w[i]*w[j] for value,w in zip(values,roots))))
    assert Hess.rank()==4
    assert (sp.Matrix(H).T*B*sp.Matrix(H))[0]==48
    root_rank=sum(sum(w[i]*H[i] for i in range(6))!=0 for w in roots)
    assert root_rank==42 and 78-4-root_rank==32
    return dict(normalized_target_reconstructed=True,source_regular=True,
                Cartan_differential_rank=4,root_direction_rank=root_rank,
                full_differential_rank=46,kernel_dimension=32,
                arithmetic='exact Q(i,sqrt(3)); independently enumerated roots')


def transfer_controls():
    v,y,z,u,c=sp.symbols('v y z u c');field=sp.QQ.frac_field(c)
    phi=v**3+(z+u)*y+z*z+u*u
    actual=sp.groebner([sp.diff(phi,v),sp.diff(phi,y),z*sp.diff(phi,z)-u*sp.diff(phi,u)-c,z*u-1],v,y,z,u,domain=field)
    expected=sp.groebner([v*v,z*z+1,u+z,2*z*y-c],v,y,z,u,domain=field)
    assert actual==expected and actual.reduce(v*v)[1]==0 and actual.reduce(v)[1]!=0
    x=sp.symbols('x');D=3*x*x-2;U=-2*x/D
    mu_squared=(U*U-4)*(2*U+x**3-2*x)**2
    formula=-4*x**6*(3*x*x-8)**2*(9*x*x-4)*(x*x-1)/D**4
    assert sp.cancel(mu_squared-formula)==0 and sp.diff(sp.cancel(mu_squared),x)!=0
    return dict(nonreduced_residual_structure_preserved=True,
                G2_target_square_identity=True,
                nonreduced_control_is_actual_character=False)


def main():
    source=json.loads((ROOT/'exceptional-adjoint-contact-verification.json').read_text())
    assert source['result']=='passed'
    for name,expected in source['source_sha256'].items():
        assert sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()==expected,name
    rows=[reconstruct(row) for row in source['contact_gradings']]
    point=native_e6_point();controls=transfer_controls()
    paths=['EXCEPTIONAL_ADJOINT_CONTACT_FIBERS.md','tools/audit_exceptional_adjoint_contact.py',
           'tools/audit_affine_branching_certificates.py','exceptional-adjoint-contact-verification.json']
    report=dict(result='passed',date='2026-10-23',uses_Sage=False,
                independent_module_models=['SL2 symmetric cube','SL6 exterior cube',
                  'Sp6 primitive exterior cube','Spin12 signed half-spin weights','E7 minuscule Weyl orbit'],
                contact_gradings=rows,E6_independent_point=point,controls=controls,
                shared_conventions='Bourbaki numbering and standard contact module names; no Sage root or character arithmetic',
                all_exceptional_Levis_classified=False,literature_priority_claim=False,
                source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in paths})
    (ROOT/'exceptional-adjoint-contact-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
