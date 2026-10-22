"""Symbolic integral checks for the symplectic projective K-image determinant."""
from hashlib import sha256
import json
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]


def check(r):
    L=sp.Symbol('L');ss=sp.symbols('s1:'+str(r+1))
    P=sp.Poly(sp.prod(L*L-s*L+1 for s in ss),L)
    inverse=-sum(P.nth(k)*L**(k-1) for k in range(1,2*r+1))
    b=sp.Poly(L-inverse,L);p=sp.Poly(1,L);columns=[]
    for _ in range(2*r):
        columns.append([p.nth(i) for i in range(2*r)])
        p=sp.rem(p*b,P)
    M=sp.Matrix(columns).T
    determinant=DomainMatrix.from_Matrix(M).det().as_expr()
    expected=sp.Integer(2)**r*sp.prod((ss[i]+ss[j])**2 for i in range(r) for j in range(i+1,r))
    assert sp.expand(determinant-expected)==0
    assert int(expected.subs({s:2 for s in ss}))==2**(r*(2*r-1))
    # The characteristic polynomial of b in this target is the monic
    # product(z^2-(s_i^2-4)), independently of the determinant formula.
    z=sp.Symbol('z');B=sp.Matrix(2*r,2*r,lambda i,j:sp.rem(b*sp.Poly(L**j,L),P).nth(i))
    actual=B.charpoly(z).as_expr()
    expected_char=sp.prod(z*z-(s*s-4) for s in ss)
    assert sp.expand(actual-expected_char)==0
    return dict(symplectic_rank=r,matrix_size=2*r,
                determinant=str(sp.factor(determinant)),
                full_source_characteristic_relation_verified=True,
                integral_augmentation_index=str(2**(r*(2*r-1))))


def main():
    cases=[]
    for r in range(1,4):
        cases.append(check(r));print('rank',r,'symbolic determinant and characteristic relation checked',flush=True)
    report=dict(result='passed',date='2026-10-22',cases=cases,
                scope='Exact symbolic projective-bundle K-ring multiplication, independent of root-discriminant or fixed-point proof; finite corroboration of the written all-rank integral formula.',
                source_sha256={'experiments/other_theories_symplectic_equivariant_pilot.py':sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
    (ROOT/'experiments/other-theories-symplectic-equivariant-pilot.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
