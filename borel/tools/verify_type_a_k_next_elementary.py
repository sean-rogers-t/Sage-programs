"""Actual elementary images and exhaustive next-layer Witt comparison."""
from hashlib import sha256
from math import comb
from pathlib import Path
import json
from sage.all import ZZ,QQ,GF,matrix,vector,lcm,identity_matrix
from verify_type_a_k_first_correction import build_image,primary
from verify_type_a_k_stable_multiplicities import expected

ROOT=Path(__file__).resolve().parents[1]


def compare(m,p):
    d=p+3
    R,basis,J,C,source=build_image(m,d,'moment')
    keep=[i for i,e in enumerate(basis) if e[0]>0]
    J=J.matrix_from_columns(keep).hermite_form(include_zero_rows=False)
    C=C.matrix_from_columns(keep).hermite_form(include_zero_rows=False)
    basis=[basis[i] for i in keep]
    index={e:i for i,e in enumerate(basis)}
    K=R.change_ring(QQ);L,*aa=K.gens();a=[K(m)]+aa;q=[K(m)]+list(map(K,source))
    b=(a[p]-a[1]**p)/p;alpha=(m-m**p)//p
    H=alpha*L**p+sum((comb(p,i)//p)*(a[i]-m**(p-i)*a[1]**i)*L**(p-i)
                    for i in range(1,p))
    external=[b*(q[3]-a[3]),a[1]*b*(q[2]-a[2])]
    vA,vB,vC=m*b*a[2]*L,m*b*a[1]**2*L,m*b*L*q[2]
    w3,w12,w111=a[3]*H,a[1]*a[2]*H,a[1]**3*H
    w1q2,w2q1,w11q1=a[1]*q[2]*H,a[2]*q[1]*H,a[1]**2*q[1]*H
    E3,E12=q[3]*H-alpha*q[d],q[1]*q[2]*H-m*alpha*q[d]
    ideal=[vA,vB,vC,w3,w12,w111,w1q2,w2q1,w11q1]
    new=[vA,vC,w3,w12,E3]
    R3=q[3]*q[p]-m*q[d]
    def vec(f):
        out=[QQ(0)]*len(basis)
        for e,c in f.dict().items():
            if e[0]:out[index[tuple(e)]]=c
        return vector(QQ,out)
    def member(f,A):
        sol=A.transpose().solve_right(vec(f))
        return all(x.denominator().valuation(p)==0 for x in sol)
    for f in [vB,w111,w11q1,vA+w2q1,QQ(m)/p*R3+3*(w1q2+w2q1),p*E3-R3]:
        assert member(f,C)
    assert member(E12+vC,J)
    for f in new[:4]:
        assert member(p*f,C)
    coordinates=J.transpose().solve_right(matrix(QQ,[vec(f) for f in new]).transpose()).transpose()
    assert all((p*x).denominator().valuation(p)==0 for x in coordinates.list())
    assert matrix(GF(p),p*coordinates).rank()==5
    rawJ=matrix(ZZ,[p*r for r in J.rows()]+[p*vec(f) for f in external+ideal+[E3,E12]])
    rawC=matrix(ZZ,[p*r for r in C.rows()]+[p*vec(f) for f in external+ideal])
    JJ=rawJ.hermite_form(include_zero_rows=False)
    CC=rawC.hermite_form(include_zero_rows=False)
    factors=primary(JJ,CC,p)
    want=[int(p**ZZ(m).valuation(p))]*expected(d)
    assert factors==want
    # Check the reduced list defines the same p-local source and ideal.
    reducedJ=matrix(ZZ,[p*r for r in J.rows()]+[p*vec(f) for f in external+new])
    reducedC=matrix(ZZ,[p*r for r in C.rows()]+[p*vec(f) for f in external+new[:4]+[w1q2+w2q1]])
    for A,B in [(JJ,reducedJ),(CC,reducedC)]:
        reduced=B.hermite_form(include_zero_rows=False)
        for X,Y in [(A,reduced),(reduced,A)]:
            solution=Y.transpose().solve_right(X.transpose())
            assert all(x.denominator().valuation(p)==0 for x in solution.list())
    # Actual elementary calculation is independently built from elementary jets.
    _,_,JE,CE,_=build_image(m,d,'actual' if m>p else 'auxiliary')
    assert primary(JE,CE,p)==want
    return dict(m=m,p=p,degree=d,factors=factors,five_new_classes_independent=True,
                all_ideal_class_identities=True,exhaustive_and_reduced_lattices_agree=True,
                elementary_model='actual' if m>p else 'auxiliary_unrooted')


def main():
    cases=[]
    for m,p in [(5,5),(10,5),(25,5),(14,7),(49,7)]:
        cases.append(compare(m,p));print(json.dumps(cases[-1]),flush=True)
    boundaries=[]
    for p in [5,7]:
        stages={}
        for model in ['auxiliary','unrooted','actual']:
            ring,basis,J,C,fs=build_image(p,p+3,model)
            stages[model]=primary(J,C,p)
            if model=='unrooted':
                L,*cs=ring.gens();F=fs[p-1]
                index={e:i for i,e in enumerate(basis)}
                def vec(f):
                    out=[ZZ(0)]*len(basis)
                    for e,c in f.dict().items():out[index[tuple(e)]]=c
                    return vector(ZZ,out)
                K=matrix(ZZ,[vec(F*L**3),vec(F*cs[0]*L**2),vec(F*cs[1]*L),vec(F*cs[0]**2*L)])
                coordinates=J.transpose().solve_right(K.transpose()).transpose()
                intersection=identity_matrix(ZZ,4)
                for column in coordinates.columns():
                    v=intersection*column
                    denominator=lcm(x.denominator() for x in v)
                    if denominator==1:continue
                    constraint=matrix(ZZ,[list(denominator*v)+[-denominator]])
                    preimage=constraint.right_kernel_matrix().matrix_from_columns(range(4))
                    intersection=(preimage*intersection).hermite_form(include_zero_rows=False)
                predicted=matrix(ZZ,[[p,-3,-6,0],[0,p,0,0],[0,0,p,0],[0,0,0,1]]).hermite_form()
                assert intersection==predicted
                target=K[3]
                assert all(x.denominator().valuation(p)==0 for x in J.transpose().solve_right(target))
                phi=vector(QQ,[0]*len(basis))
                for monomial,value in [(cs[0]**2*L**(p+1),QQ(1)/p),
                                       (cs[0]*L**(p+2),QQ(p+2)/p**2)]:
                    phi[index[tuple(monomial.exponents()[0])]]=value
                assert all(x.denominator()==1 for x in C*phi)
                assert target*phi==QQ(1)/p
                root_intersection=[list(map(int,row)) for row in intersection]
        assert [len(stages[x]) for x in ['auxiliary','unrooted','actual']]==[
            expected(p+3),expected(p+3)-1,expected(p+3)-2]
        boundaries.append(dict(p=p,stages=stages,root_intersection=root_intersection,
                               literal_source_membership=True,root_class_functional_value='1/%d'%p,
                               full_coefficient_ideal_integrality=True))
    names=['TYPE_A_K_NEXT_ELEMENTARY.md','TYPE_A_K_NEXT_LAYER_OBSTRUCTION.md',
           'TYPE_A_K_FIRST_CORRECTION.md','tools/verify_type_a_k_next_elementary.py',
           'tools/verify_type_a_k_first_correction.py','tools/verify_type_a_k_stable_multiplicities.py']
    report=dict(result='passed',comparisons=cases,root_boundary_controls=boundaries,
                all_rank_claim_from_finite_checks=False,root_boundary_theorem=True,
                source_sha256={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names})
    (ROOT/'type-a-k-next-elementary-verification.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
