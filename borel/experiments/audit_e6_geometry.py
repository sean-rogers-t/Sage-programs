"""Independent E6 matrix audit and exact unisolvent certificates.

Build the 27 and 78 dimensional representations independently, compare all
78 differential blocks coefficientwise, and check every quotient coordinate
at a provably unisolvent set of exact rational torus points.
"""
from pathlib import Path
from time import perf_counter
from collections import Counter
from itertools import combinations
import json
from sage.all import QQ, GF, matrix, vector, prod, WeylCharacterRing
from sage_borel import (springer_character_map, SpringerGeometry,
                       e6_adjoint_quotient, e6_cayley_determinant)
from e6_matrix_oracle import matrices, adjoint_matrices, C
from e6_geometry_support import candidates

ROOT=Path(__file__).resolve().parents[1]


def characters(point,oracle,field):
    values=[prod(field(t)**a for t,a in zip(point,w)) for w in oracle['weights']]
    dual=[1/v for v in values]
    p=[sum(v**k for v in values) for k in (1,2,3)]
    d=[sum(v**k for v in dual) for k in (1,2)]
    roots=[prod(field(t)**a for t,a in zip(point,w)) for w in oracle['labels'][6:]]
    return (p[0],6+sum(roots),(p[0]**2-p[1])/2,
            (p[0]**3-3*p[0]*p[1]+2*p[2])/6,(d[0]**2-d[1])/2,d[0])


def monomials(chars,exponents):
    powers=[[c**k for k in range(max(e[i] for e in exponents)+1)] for i,c in enumerate(chars)]
    return [prod(powers[i][a] for i,a in enumerate(e)) for e in exponents]


def cartan_trace_data(basis,labels,ring):
    diagonal=[[h[k,k] for k in range(basis[0].nrows())] for h in basis[:6]]
    def laurent(coefficients):
        D={}
        for w,c in zip(labels,coefficients):
            D[w]=D.get(w,0)+c
        return ring(D)
    gram=matrix(QQ,6,6,lambda i,j:sum(a*b for a,b in zip(diagonal[i],diagonal[j])))
    A=matrix(ring,6,6,lambda i,j:laurent([a*b for a,b in zip(diagonal[i],diagonal[j])]))
    return diagonal,gram,gram.inverse()*A,laurent


def audit_representation(name,oracle,basis,labels):
    start=perf_counter()
    highest=(1,0,0,0,0,0) if name=='minimal' else (0,1,0,0,0,0)
    s=springer_character_map('E6',highest)
    g=SpringerGeometry(s)
    diagonal,gram,J,laurent=cartan_trace_data(basis,labels,s.torus_ring)
    assert gram==s.trace_gram and J==g.torus_jacobian()
    lookup={a:i for i,a in enumerate(oracle['labels']) if i>=6}
    root_factors=dict(g.root_factors())
    root_blocks=0
    root_polynomials=[]
    for alpha,index in lookup.items():
        opposite=lookup[tuple(-a for a in alpha)]
        product=basis[index]*basis[opposite]
        coefficients=[product[k,k] for k in range(product.nrows())]
        norm=sum(coefficients)
        assert norm
        actual=laurent([a/norm for a in coefficients])
        if alpha in root_factors:
            expected=prod(t**a for t,a in zip(s.torus_ring.gens(),alpha))*root_factors[alpha]
        else:
            expected=root_factors[tuple(-a for a in alpha)]
        assert actual==expected
        root_polynomials.append(actual)
        root_blocks+=1
    quotient=e6_adjoint_quotient(representation=name)
    determinant=e6_cayley_determinant(representation=name)
    polynomials=dict(zip(quotient.degrees,quotient.character_polynomials))
    certificates=[]
    for mode in ('quotient','determinant'):
        data=json.loads((ROOT/('experiments/e6-%s-%s-interpolation.json'%(name,mode))).read_text())
        for block in data['cases']:
            exponents=candidates(highest,max(block['degrees']))
            assert len(exponents)==block['candidate_count']
            allowed=set(exponents)
            inputs=([polynomials[d] for d in block['degrees']] if mode=='quotient'
                    else [determinant.torus_character_polynomial])
            assert all(tuple(e) in allowed for f in inputs for e in f.dict())
            finite=GF(block['modulus'])
            modular_rows=[]
            equations=0
            for point in block['points']:
                cs=characters(point,oracle,QQ)
                modular_rows.append(monomials(characters(point,oracle,finite),exponents))
                T=[prod(QQ(t)**a for t,a in zip(point,w)) for w in labels]
                if mode=='quotient':
                    moments=vector(QQ,[sum(a*v for a,v in zip(row,T)) for row in diagonal])
                    h=gram.inverse()*moments
                    # Direct diagonal of the independently constructed 27d
                    # Lie matrix rho_27(h), then trace its powers.
                    eigenvalues=[sum(a*b for a,b in zip(w,h)) for w in oracle['weights']]
                    for degree in block['degrees']:
                        normalization=12 if degree==2 else 1
                        assert polynomials[degree](*cs)==sum(v**degree for v in eigenvalues)/normalization
                        equations+=1
                else:
                    direct=matrix(QQ,6,6,lambda i,j:sum(diagonal[i][k]*diagonal[j][k]*T[k]
                                   for k in range(len(T)))).det()/gram.det()
                    assert determinant.torus_character_polynomial(*cs)==direct
                    equations+=1
            rank=int(matrix(finite,modular_rows).rank())
            assert rank==len(exponents)==len(block['points'])
            certificates.append(dict(mode=mode,degrees=block['degrees'],rank=rank,
                                     exact_scalar_identities=equations))
    assert determinant.evaluate([1]*6)==1
    finite=GF(65537)
    point=tuple(map(finite,[2,3,5,7,11,13]))
    def evaluate(f):
        return sum(finite(c)*prod(t**a for t,a in zip(point,e)) for e,c in f.dict().items())
    independent=matrix([[evaluate(f) for f in row] for row in J]).det()*prod(evaluate(f) for f in root_polynomials)
    assert independent==determinant.evaluate(point)
    return dict(representation=name,dimension=basis[0].nrows(),
                coefficientwise_cartan_block=True,coefficientwise_root_blocks=root_blocks,
                full_differential_dimension=78,normalization_at_identity=1,
                coordinate_terms={str(d):len(f.dict()) for d,f in polynomials.items()},
                cartan_determinant_terms=len(determinant.torus_character_polynomial.dict()),
                finite_field_determinant=dict(prime=65537,point=[2,3,5,7,11,13],value=int(independent)),
                unisolvent_certificates=certificates,seconds=perf_counter()-start)


def invariant_certificate(oracle):
    # Basicness uses the fundamental degrees and one exact nonzero Jacobian.
    # Invariance follows from a sum over the entire minuscule weight orbit.
    h=C.change_ring(QQ).inverse()*vector(QQ,[1,2,3,4,5,6])
    degrees=(2,5,6,8,9,12)
    J=matrix(QQ,[[sum(d*w[i]*sum(a*b for a,b in zip(w,h))**(d-1)
                     for w in oracle['weights'])/(12 if d==2 else 1)
                  for i in range(6)] for d in degrees])
    determinant=J.det()
    assert determinant
    positive=[a for a in oracle['labels'][6:]
              if all(v>=0 for v in C.change_ring(QQ).inverse()*vector(a))]
    discriminant=prod(sum(a*b for a,b in zip(alpha,h)) for alpha in positive)
    assert len(positive)==36 and discriminant
    return dict(degrees=degrees,point=list(map(str,h)),jacobian=str(determinant),
                jacobian_over_positive_root_product=str(determinant/discriminant))


def main():
    start=perf_counter()
    oracle=matrices()
    adjoint,brackets=adjoint_matrices(oracle)
    WCR=WeylCharacterRing('E6',style='coroots')
    minimal=WCR((1,0,0,0,0,0))
    assert minimal.exterior_power(2)==WCR((0,0,1,0,0,0))
    assert minimal.exterior_power(3)==WCR((0,0,0,1,0,0))
    coroots=WCR.space().simple_coroots()
    def weight_counter(character):
        return Counter({tuple(int(w.inner_product(coroots[i+1])) for i in range(6)):int(m)
                        for w,m in character.weight_multiplicities().items()})
    assert weight_counter(WCR((0,0,0,0,0,1)))==Counter(tuple(-a for a in w) for w in oracle['weights'])
    assert weight_counter(WCR((0,0,0,0,1,0)))==Counter(
        tuple(-a-b for a,b in zip(x,y)) for x,y in combinations(oracle['weights'],2))
    assert weight_counter(WCR((0,1,0,0,0,0)))==Counter(oracle['labels'])
    cases=[]
    for name,basis,labels in (('minimal',oracle['basis'],oracle['weights']),
                              ('adjoint',adjoint,oracle['labels'])):
        result=audit_representation(name,oracle,basis,labels)
        cases.append(result)
        print(name,result,flush=True)
    report=dict(result='passed',simple_generator_relations=oracle['relations'],
                independently_verified_commutators=brackets,
                fundamental_exterior_character_identifications=True,
                invariant_certificate=invariant_certificate(oracle),cases=cases,
                seconds=perf_counter()-start)
    (ROOT/'experiments/e6-geometry-audit.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
