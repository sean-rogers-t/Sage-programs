"""Retain both monic quadratic relations as a free rank-four algebra.

Avoids substituting chi3=chi1^2-... into large powers before elimination.
The five remaining equations act on a four-element basis over five variables.
"""
from pathlib import Path
import json
import re
import sys
from sage.all import GF, PolynomialRing, matrix, identity_matrix, zero_matrix, block_matrix
from sage_borel import e7_adjoint_quotient
from exceptional_degree_fibers import prepare_fiber
from e7_adjoint_module_fiber import module_length, preflight


def tower_presentation(R,equations,targets):
    F=R.base_ring()
    B=PolynomialRing(F,5,names=['x%d'%(i+1) for i in range(5)])
    u,b,v,d,z=B.gens()
    S=PolynomialRing(B,'s')
    s=S.gen()
    original=e7_adjoint_quotient(representation='adjoint').character_polynomials[0]
    P=PolynomialRing(F,7,names=['h%d'%(i+1) for i in range(7)])
    first=P({tuple(e):F(c) for e,c in original.dict().items()})-targets[0]
    assert all(first.degree(P.gen(i))==0 for i in (1,3,4,6))
    qraw=P.hom([s,0,b,0,0,d,0],S)(first)
    lc=qraw.leading_coefficient()
    assert qraw.degree()==2 and lc.is_constant() and lc
    qraw=S(qraw/lc)
    shift=qraw[1]/2
    assert shift.is_constant()
    q1=qraw(s-shift)
    assert q1[1]==0 and q1[2]==1
    k=-q1[0]
    assert k.total_degree()==1
    I2=identity_matrix(B,2)
    Z2=zero_matrix(B,2)
    A2=matrix(B,[[0,k],[1,0]])
    original_A2=A2-shift*I2
    W=PolynomialRing(S,'c')
    c=W.gen()
    hom=R.hom([s-shift,u,c,v,d,z],W)
    f6=hom(equations[0])
    lc6=f6.leading_coefficient()
    assert f6.degree()==2 and lc6.degree()==0 and lc6[0].is_constant() and lc6
    f6=W(f6/lc6)
    assert f6[2]==1
    def evaluate_coefficient(f):
        answer=zero_matrix(B,2)
        for coefficient in reversed(S(f).list()):
            answer=answer*A2+coefficient*I2
        return answer
    P2=evaluate_coefficient(f6[1])
    R2=evaluate_coefficient(f6[0])
    D2=P2*P2/4-R2
    A4=block_matrix([[original_A2,Z2],[Z2,original_A2]])
    C4=block_matrix([[Z2,D2],[I2,Z2]])-block_matrix([[P2/2,Z2],[Z2,P2/2]])
    assert A4*C4==C4*A4
    I4=identity_matrix(B,4)
    assert A4*A4+qraw[1]*A4+qraw[0]*I4==zero_matrix(B,4)
    apowers=[I4]
    cpowers=[I4]
    for _ in range(max(f.degree(R.gen(0)) for f in equations)):
        apowers.append(apowers[-1]*A4)
    for _ in range(max(f.degree(R.gen(2)) for f in equations)):
        cpowers.append(cpowers[-1]*C4)
    products={}
    def evaluate(f):
        grouped={}
        for e,coefficient in f.dict().items():
            pair=(e[0],e[2])
            exponent=(e[1],0,e[3],e[4],e[5])
            coefficients=grouped.setdefault(pair,{})
            coefficients[exponent]=coefficients.get(exponent,F.zero())+coefficient
        answer=zero_matrix(B,4)
        for pair,coefficients in grouped.items():
            if pair not in products:
                products[pair]=apowers[pair[0]]*cpowers[pair[1]]
            answer+=B(coefficients)*products[pair]
        return answer
    assert evaluate(equations[0])==zero_matrix(B,4)
    columns=[]
    degrees=[]
    for f in equations[1:]:
        M=evaluate(f)
        assert M*A4==A4*M and M*C4==C4*M
        columns.extend(tuple(M.column(i)) for i in range(4))
        degrees.append(max(int(g.total_degree()) for g in M.list() if g))
    info=dict(tower_ranks=[2,2],module_basis=['1','s','c_prime','s*c_prime'],
        first_character_shift=str(shift),first_square_rhs_degree=1,
        residual_matrix_coefficient_degrees=degrees,
        commuting_regular_representation=True,first_target_matrix_identity=True,
        second_target_matrix_identity=True,remaining_multiplication_commutators=10)
    return B,columns,info


def main():
    prime=int(sys.argv[1]) if len(sys.argv)>1 else 32003
    seed=int(sys.argv[2]) if len(sys.argv)>2 else 20261020
    order=sys.argv[3] if len(sys.argv)>3 else 'weighted'
    assert order in ('dp','weighted')
    preflight()
    smoke=PolynomialRing(GF(32003),2,names=('x','y'))
    columns=[]
    for variable in smoke.gens():
        for i in range(4):
            column=[smoke.zero()]*4
            column[i]=variable
            columns.append(tuple(column))
    assert module_length(smoke,columns,weights=(2,2,0,1,1,2),
                         label='rank-four-nilpotent')['chart_fiber_length']==4
    R,equations,targets,eliminated,chart=prepare_fiber('E7','adjoint',prime,seed,mode='eliminate-center')
    B,columns,tower=tower_presentation(R,equations,targets)
    print('Tower matrix coefficient degrees:',tower['residual_matrix_coefficient_degrees'],flush=True)
    weights=(1,4,2,2,2,0,2,5,7) if order=='weighted' else None
    report=module_length(B,columns,weights=weights,label='e7-adjoint-tower')
    report.update(status='specialized-fiber-only',cartan_type='E7',representation='adjoint',
        characteristic=prime,seed=seed,targets=list(map(str,targets)),
        fiber_length=chart['cover_degree']*report['chart_fiber_length'],
        eliminated_coordinate=eliminated,center_chart=chart,tower_reduction=tower,
        preflight_lengths=[1,4],rank_four_preflight_length=4)
    path=Path(__file__).with_name('e7-adjoint-fiber-tower-%d-%d-%s.json'%(prime,seed,order))
    encoded=json.dumps(report,indent=2)
    encoded=re.sub(r'\[\n\s+(-?\d+(?:,\n\s+-?\d+)*)\n\s+\]',
                   lambda m:'['+','.join(re.findall(r'-?\d+',m.group(1)))+']',encoded)
    path.write_text(encoded+'\n')
    print({k:v for k,v in report.items() if k!='module_leading_terms'},flush=True)


if __name__=='__main__':
    main()
