"""Check the formal-inverse height mechanism on symplectic projective flags.

Uses p-typical logarithms for heights 1--5 over Q, checks 2-integrality,
then reduces the inverse series modulo two. An independent supersingular
elliptic formal group supplies the height-two check.
"""
from hashlib import sha256
import json
from pathlib import Path
from sage.all import EllipticCurve,GF,PowerSeriesRing,QQ,matrix

ROOT=Path(__file__).resolve().parents[1]


def bit_product(a,b,precision):
    out=0
    while a:
        bit=a & -a;shift=bit.bit_length()-1;out ^= b << shift;a ^= bit
    return out & ((1 << precision)-1)


def binary_rank(columns):
    piv={}
    for a in columns:
        while a:
            i=a.bit_length()-1
            if i in piv:a ^= piv[i]
            else:piv[i]=a;break
    return len(piv)


def main():
    precision=130;Q=PowerSeriesRing(QQ,'t',default_prec=precision);t=Q.gen();cases=[]
    inverse_h2=None
    for height in range(1,6):
        q=2**height;log=Q.zero();power=1;k=0
        while power<precision:
            log+=t**power/(2**k);power*=q;k+=1
        log=log.add_bigoh(precision)
        inv=log.reverse()(-log)
        assert all(c.denominator()%2 for c in inv.list())
        inverse_bits=sum((int(c.numerator()) & 1) << j for j,c in enumerate(inv.list()))
        b=(1 << 1) ^ inverse_bits
        assert (b & -b).bit_length()-1==q
        if height==2:inverse_h2=inverse_bits
        rows=[]
        for m in range(2,130,2):
            value=1;columns=[]
            for _ in range(m):columns.append(value);value=bit_product(value,b,m)
            rank=binary_rank(columns);predicted=(m-1)//q+1
            assert rank==predicted
            first_zero=next(j for j,c in enumerate(columns+[value]) if not c)
            assert first_zero==predicted
            rows.append(dict(symplectic_rank=m//2,projective_dimension=m-1,image_rank=rank,nilpotence_order=first_zero))
        cases.append(dict(height=height,formal_inverse_first_difference_degree=q,
                          inverse_mod2_exponents=[j for j in range(precision) if (inverse_bits >> j)&1],checks=rows))
    E=EllipticCurve(GF(2),[0,0,1,0,0]);inv=E.formal_group().inverse(prec=precision)
    elliptic_bits=sum(int(c) << j for j,c in enumerate(inv.list()))
    assert (elliptic_bits ^ 2) & 15==0 and ((elliptic_bits >> 4)&1)==1
    # Agreement of height, not an assumption that the two coordinates or full
    # formal group laws are identical.
    elliptic_checks=[]
    for m in range(2,130,2):
        b=elliptic_bits ^ 2;value=1;cols=[]
        for _ in range(m):cols.append(value);value=bit_product(value,b,m)
        assert binary_rank(cols)==(m-1)//4+1
        elliptic_checks.append(dict(projective_dimension=m-1,image_rank=binary_rank(cols)))
    report=dict(result='passed',date='2026-10-22',precision=precision,cases=cases,
                supersingular_elliptic_curve='y^2+y=x^3 over F2',independent_height_two_checks=elliptic_checks,
                scope='Finite formal-series corroboration of the written valuation/nilpotence proof. The general theorem assumes a unit leading coefficient of the 2-series at degree 2^h.',
                source_sha256={'experiments/other_theories_height_pilot.py':sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
    (ROOT/'experiments/other-theories-height-pilot.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASSED',len(cases)*64,'height/rank cases and',len(elliptic_checks),'independent elliptic cases')


if __name__=='__main__':main()
