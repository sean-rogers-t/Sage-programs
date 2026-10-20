"""Independent 56/133-dimensional matrices and exact E7 certificates.

The representation oracle uses no character-ring or interpolation code.
Freudenthal independently reconstructs the 912-dimensional character.
Every coordinate is checked on a proved unisolvent set, not just probes.
"""
from pathlib import Path
from time import perf_counter
from collections import Counter
import json
from sage.all import QQ, ZZ, GF, matrix, vector, prod, WeylCharacterRing
from sage.misc.cython import cython_import
from sage_borel import (springer_character_map, SpringerGeometry,
                       e7_adjoint_quotient, e7_cayley_determinant)
from e7_matrix_oracle import matrices, adjoint_matrices, C
from e7_geometry_support import candidates, parent_plan, MARKS
from audit_trace_source_survey import weights_from_strings

ROOT=Path(__file__).resolve().parents[1]


def character_values(point,oracle,extra):
    powers=[{k:QQ(t)**k for k in range(-4,5)} for t in point]
    monomial=lambda w:prod(powers[i][int(a)] for i,a in enumerate(w))
    values=[monomial(w) for w in oracle['weights']]
    p1,p2,p3,p4=(sum(v**k for v in values) for k in (1,2,3,4))
    roots=[monomial(w) for w in oracle['labels'][7:]]
    a1=7+sum(roots)
    a2=7+sum(v*v for v in roots)
    e2=(p1*p1-p2)/2
    e3=(p1**3-3*p1*p2+2*p3)/6
    e4=(p1**4-6*p1*p1*p2+3*p2*p2+8*p1*p3-6*p4)/24
    return (a1,sum(m*monomial(w) for w,m in extra.items()),
            (a1*a1-a2)/2-a1,e4-e2,e3-p1,e2-1,p1)


def cleared_row(point,chars,exponents,plan):
    B=ZZ(prod(abs(t) for t in point))
    numerators=[ZZ(c*B**a) for c,a in zip(chars,MARKS)]
    parents,factors,costs=plan
    K=max(costs)
    powers=[B**k for k in range(K+1)]
    values=[]
    for e,p,(i,j) in zip(exponents,parents,factors):
        # Independently verify the exponent identity in the multiplication plan.
        base=(0,)*7 if p<0 else exponents[p]
        expected=list(base)
        if i>=0:
            expected[i]+=1
        if j>=0:
            expected[j]+=1
        assert tuple(expected)==e
        value=ZZ.one() if p<0 else values[p]
        if i>=0:
            value*=numerators[i]
        if j>=0:
            value*=numerators[j]
        values.append(value)
    return [v*powers[K-c] for v,c in zip(values,costs)],B**K


def cartan_trace_data(basis,labels,ring):
    diagonal=[[h[k,k] for k in range(basis[0].nrows())] for h in basis[:7]]
    def laurent(coefficients):
        D={}
        for w,c in zip(labels,coefficients):
            D[w]=D.get(w,0)+c
        return ring(D)
    gram=matrix(QQ,7,7,lambda i,j:sum(a*b for a,b in zip(diagonal[i],diagonal[j])))
    A=matrix(ring,7,7,lambda i,j:laurent([a*b for a,b in zip(diagonal[i],diagonal[j])]))
    return diagonal,gram,gram.inverse()*A,laurent


def audit_representation(name,oracle,basis,labels,extra):
    start=perf_counter()
    highest=(0,0,0,0,0,0,1) if name=='minimal' else (1,0,0,0,0,0,0)
    s=springer_character_map('E7',highest)
    g=SpringerGeometry(s)
    diagonal,gram,J,laurent=cartan_trace_data(basis,labels,s.torus_ring)
    assert gram==s.trace_gram and J==g.torus_jacobian()
    lookup={a:i for i,a in enumerate(oracle['labels']) if i>=7}
    root_factors=dict(g.root_factors())
    root_polynomials=[]
    for alpha,index in lookup.items():
        opposite=lookup[tuple(-a for a in alpha)]
        product=basis[index]*basis[opposite]
        coefficients=[product[k,k] for k in range(product.nrows())]
        norm=sum(coefficients)
        assert norm
        actual=laurent([a/norm for a in coefficients])
        expected=(prod(t**a for t,a in zip(s.torus_ring.gens(),alpha))*root_factors[alpha]
                  if alpha in root_factors else root_factors[tuple(-a for a in alpha)])
        assert actual==expected
        root_polynomials.append(actual)
    print(name,'all 49 Cartan entries and 126 root blocks agree',flush=True)
    quotient=e7_adjoint_quotient(representation=name)
    determinant=e7_cayley_determinant(representation=name)
    polynomials=dict(zip(quotient.degrees,quotient.character_polynomials))
    certificates=[]
    kernel=cython_import(str(ROOT/'experiments/e7_interpolation_kernels.pyx'),use_cache=True)
    for block_name in ('low','high','determinant'):
        data=json.loads((ROOT/('experiments/e7-%s-%s-interpolation.json'%(name,block_name))).read_text())
        block=data['cases'][0]
        exponents=candidates(highest,max(block['degrees']))
        assert len(exponents)==block['candidate_count']
        allowed=set(exponents)
        plan=parent_plan(exponents)
        inputs=([polynomials[d] for d in block['degrees']] if block_name!='determinant'
                else [determinant.torus_character_polynomial])
        coefficient_vectors=[]
        for f,record in zip(inputs,block['coordinates']):
            normalization=24 if block_name!='determinant' and record['degree']==2 else 1
            assert all(tuple(e) in allowed for e in f.dict())
            integers={tuple(e):ZZ(c) for e,c in record['terms']}
            assert f==f.parent()({e:QQ(c)/(record['scale']*normalization) for e,c in integers.items()})
            coefficient_vectors.append([integers.get(e,ZZ.zero()) for e in exponents])
        finite=GF(block['rank_prime'])
        modular_chars=[]
        equations=0
        for k,point in enumerate(block['points']):
            cs=character_values(point,oracle,extra)
            modular_chars.append([int(finite(c)) for c in cs])
            row,denominator=cleared_row(point,cs,exponents,plan)
            T=[prod(QQ(t)**a for t,a in zip(point,w)) for w in labels]
            if block_name!='determinant':
                moments=vector(QQ,[sum(a*v for a,v in zip(diag,T)) for diag in diagonal])
                h=gram.inverse()*moments
                # Direct eigenvalues of the independently certified minuscule
                # Lie matrix rho_56(h), retaining exactly scaled integral traces.
                scale=24 if name=='minimal' else 36
                eigenvalues=[scale*sum(a*b for a,b in zip(w,h)) for w in oracle['weights']]
                expected=[sum(v**d for v in eigenvalues) for d in block['degrees']]
            else:
                expected=[matrix(QQ,7,7,lambda i,j:sum(diagonal[i][u]*diagonal[j][u]*T[u]
                                  for u in range(len(T)))).det()]
            for coefficients,value in zip(coefficient_vectors,expected):
                assert sum(a*b for a,b in zip(row,coefficients))==denominator*value
                equations+=1
            if (k+1)%1000==0:
                print(name,block_name,'exact matrix-trace equations at',k+1,'points',flush=True)
        array=kernel.evaluation_rows(modular_chars,plan[0],plan[1],int(finite.characteristic()))
        A=matrix(finite,array)
        del array
        rank=int(A.rank())
        del A
        assert rank==len(exponents)==len(block['points'])
        certificates.append(dict(block=block_name,degrees=block['degrees'],rank=rank,
                                 exact_scalar_identities=equations))
        print(name,block_name,'independent exact certificate passed',flush=True)
    assert determinant.evaluate([1]*7)==1
    finite=GF(65537)
    point=tuple(map(finite,[2,3,5,7,11,13,17]))
    def evaluate(f):
        return sum(finite(c)*prod(t**a for t,a in zip(point,e)) for e,c in f.dict().items())
    independent=matrix([[evaluate(f) for f in row] for row in J]).det()*prod(evaluate(f) for f in root_polynomials)
    assert independent==determinant.evaluate(point)
    return dict(representation=name,dimension=basis[0].nrows(),
                coefficientwise_cartan_entries=49,coefficientwise_root_blocks=126,
                full_differential_dimension=133,normalization_at_identity=1,
                coordinate_terms={str(d):len(f.dict()) for d,f in polynomials.items()},
                cartan_determinant_terms=len(determinant.torus_character_polynomial.dict()),
                finite_field_determinant=dict(prime=65537,point=[2,3,5,7,11,13,17],value=int(independent)),
                unisolvent_certificates=certificates,seconds=perf_counter()-start)


def invariant_certificate(oracle):
    h=C.change_ring(QQ).inverse()*vector(QQ,[1,2,3,4,5,6,7])
    degrees=(2,6,8,10,12,14,18)
    J=matrix(QQ,[[sum(d*w[i]*sum(a*b for a,b in zip(w,h))**(d-1)
                     for w in oracle['weights'])/(24 if d==2 else 1)
                  for i in range(7)] for d in degrees])
    determinant=J.det()
    assert determinant
    positive=[a for a in oracle['labels'][7:]
              if all(v>=0 for v in C.change_ring(QQ).inverse()*vector(a))]
    discriminant=prod(sum(a*b for a,b in zip(alpha,h)) for alpha in positive)
    assert len(positive)==63 and discriminant
    return dict(degrees=degrees,point=list(map(str,h)),jacobian=str(determinant),
                jacobian_over_positive_root_product=str(determinant/discriminant))


def main():
    start=perf_counter()
    oracle=matrices()
    adjoint,brackets=adjoint_matrices(oracle)
    extra=weights_from_strings('E7',(0,1,0,0,0,0,0))
    data=json.loads((ROOT/'src/sage_borel/data/e7_geometry.json').read_text())
    assert dict((tuple(w),m) for w,m in data['weights_912'])==extra
    assert sum(extra.values())==912 and len(extra)==632
    assert Counter(map(tuple,data['weights_56']))==Counter(oracle['weights'])
    assert Counter(map(tuple,data['roots']))==Counter(oracle['labels'][7:])
    R=WeylCharacterRing('E7',style='coroots')
    fund=lambda i:R(tuple(int(j==i) for j in range(7)))
    minimal=fund(6)
    assert minimal.exterior_power(2)==fund(5)+1
    assert minimal.exterior_power(3)==fund(4)+fund(6)
    assert minimal.exterior_power(4)==fund(3)+fund(5)+1
    assert fund(0).exterior_power(2)==fund(2)+fund(0)
    print('independent matrices, 912 weights, and exterior identities passed',flush=True)
    cases=[]
    for name,basis,labels in (('minimal',oracle['basis'],oracle['weights']),
                              ('adjoint',adjoint,oracle['labels'])):
        cases.append(audit_representation(name,oracle,basis,labels,extra))
    report=dict(result='passed',simple_generator_relations=oracle['relations'],
                independently_verified_commutators=brackets,
                independently_reconstructed_912_weights=len(extra),
                fundamental_exterior_character_identifications=True,
                invariant_certificate=invariant_certificate(oracle),cases=cases,
                seconds=perf_counter()-start)
    (ROOT/'experiments/e7-geometry-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    main()
