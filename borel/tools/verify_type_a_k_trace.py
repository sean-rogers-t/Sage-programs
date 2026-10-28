"""Independent companion-matrix and symbolic checks of the trace lemma.

This does not import the recurrence builder or the integer witness auditor.
Matrices are multiplied over an integral polynomial ring, without choosing
roots, discarding nilpotents, or dividing by a discriminant.
"""
from pathlib import Path
from hashlib import sha256
from math import comb
import json
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import ZZ
from sympy.matrices.normalforms import smith_normal_form

ROOT=Path(__file__).resolve().parents[1]


def main():
    cases=[]
    for m in range(2,13):
        result=ring(['L']+['c%d'%j for j in range(1,m+1)],ZZ)
        domain=result[0];L=result[1];c=[domain.one]+list(result[2:])
        zero=domain.zero;one=domain.one
        # X times the basis 1,X,...,X^(m-1), using F(X)=0.
        M=[[zero for _ in range(m)] for _ in range(m)]
        for j in range(m-1):M[j+1][j]=one
        for k in range(1,m+1):M[m-k][m-1]=(-1)**(k+1)*c[k]
        D=[[(L if i==j else zero)-M[i][j] for j in range(m)] for i in range(m)]
        support=[[(k,D[k][j]) for k in range(m) if D[k][j]] for j in range(m)]
        B=[[one if i==j else zero for j in range(m)] for i in range(m)]
        f=[one]+[sum(((-1)**k*comb(m-k,j-k)*c[k]*L**(j-k)
                      for k in range(j+1)),zero) for j in range(1,m+1)]
        moments=[domain(m)]
        for r in range(1,2*m+3):
            B=[[sum((B[i][k]*value for k,value in support[j]),zero)
                for j in range(m)] for i in range(m)]
            actual=sum((B[i][i] for i in range(m)),zero)
            predicted=sum(((-1)**(k-1)*f[k]*moments[r-k]
                           for k in range(1,min(r-1,m)+1)),zero)
            if r<=m:predicted+=(-1)**(r-1)*r*f[r]
            assert actual==predicted,(m,r)
            moments.append(predicted)
        cases.append(dict(m=m,moments_checked=2*m+2,
                          coefficient_ring='Z[L,c_1,...,c_m]',
                          f_m_retained_until_root_quotient=True))
        print('companion trace passed:',m,flush=True)

    # Verify Hasse translation with rank itself indeterminate. The finite
    # index range tests the implementation; the note proves the identity
    # for all indices by the factorial/binomial identity.
    n=sp.Symbol('n');symbolic=0
    def choose(poly,k):
        return sp.prod(poly-i for i in range(k))/sp.factorial(k)
    for j in range(11):
        for k in range(j+1):
            for r in range(j-k+1):
                lhs=choose(n-k,j-k)*comb(j-k,r)
                rhs=choose(n-j+r,r)*choose(n-k,j-k-r)
                assert sp.expand(lhs-rhs)==0
                symbolic+=1

    # The invariant-ideal shortcut has larger torsion already in A2.
    a,b=sp.symbols('a b');q=a*a-a*b+b*b
    cubic=(a+b)*(2*a-b)*(a-2*b)
    rows=[]
    for value in (a*cubic,b*cubic):
        remainder=sp.rem(value,q,a)
        rows.append([remainder.coeff(a).coeff(b,3),remainder.coeff(b,4)])
    assert rows==[[-3,6],[-6,3]]
    matrix=sp.Matrix(rows)
    assert abs(matrix.det())==27
    assert sp.gcd_list(list(matrix))==3
    control=json.loads((ROOT/'type-a-k-flag-transfer-verification.json').read_text(encoding='utf-8'))
    assert control['result']=='passed' and control['full_integral_closure']
    a,b=sp.symbols('a b')
    parse=lambda rows:sp.Matrix([[sp.sympify(x,locals={'a':a,'b':b}) for x in row] for row in rows])
    B=parse(control['presentation']);R=parse(control['syzygies'])
    assert (B*R).applyfunc(sp.expand)==sp.zeros(B.rows,R.cols)
    B0=B.subs({a:3,b:3});R0=R.subs({a:3,b:3})
    assert B0.rank()==6 and R0.rank()==B.cols-6
    source_smith=smith_normal_form(R0,domain=ZZ)
    ordinary_smith=smith_normal_form(B0,domain=ZZ)
    assert [abs(int(source_smith[i,i])) for i in range(min(source_smith.shape)) if abs(source_smith[i,i])>1]==[3]*4
    assert [abs(int(ordinary_smith[i,i])) for i in range(6)]==[1,1,3,3,3,3]
    names=['TYPE_A_K_INTEGRAL_RECURRENCE.md','tools/verify_type_a_k_trace.py',
           'type-a-k-flag-transfer-verification.json','tools/verify_type_a_k_flag_transfer.py']
    report=dict(result='passed',method='Independent integral companion-matrix powers and symbolic-rank binomial checks.',
        cases=cases,symbolic_rank_hasse_identities=symbolic,
        invariant_shortcut_counterexample_smith=[3,9],
        independent_full_flag_abelian_factor_replay=[3,3,3,3],
        written_proof_scope=dict(machine_formalized=False,external_review=False,
            actual_image_used=True,no_reduced_fiber_replacement=True,
            no_base_change_injectivity_assumed=True,arbitrary_words_covered_by_written_proof=True,
            arbitrary_coefficient_family_requires_augmentation_to_Z_and_coalesced_polynomial=True),
        source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in names})
    (ROOT/'type-a-k-trace-verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('passed',sum(c['moments_checked'] for c in cases),'matrix moments and',symbolic,'symbolic-rank identities')


if __name__=='__main__':main()
