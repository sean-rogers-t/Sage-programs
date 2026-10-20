"""Literal E7/P7 Levi-character form: all-degree certificate and full image.

Uses integral A5 permutation coordinates to compress root-lattice coefficients.
Run from borel/: PYTHONPATH=src sage -python experiments/e7_character_form.py
"""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
import json
from sage.all import (QQ,ZZ,CartanMatrix,WeylCharacterRing,RootSystem,matrix,vector,
    PolynomialRing,LaurentPolynomialRing,prod,binomial,identity_matrix,diagonal_matrix,gcd)
from sage_borel import exceptional_borel
from e6_character_form import integral_preimage,certify_preimage
from exceptional_minuscule_images import monomials,minor_smith
from exceptional_localization import FixedPointBorel
from exceptional_published_rings import PublishedMinusculeRing

ROOT=Path(__file__).resolve().parents[1]
NAMES=('omega7','psi2','psi5','psi6','psi8','psi9','psi12')
DEGREES=(1,2,5,6,8,9,12)


def compressed_integral_preimage(M):
    """Exact row lattice, preconditioned and accumulated in small batches.

    Column contents are accounted for in the returned basis. Early return
    requires checking every original row against the candidate basis, not
    just the accumulated subset; the later minor certificate is independent.
    """
    scales=diagonal_matrix(QQ,[1/gcd(col) for col in M.columns()])
    A=matrix(ZZ,M*scales)
    rows=sorted(A.rows(),key=lambda row:(max(abs(c) for c in row),sum(abs(c) for c in row)))
    H=matrix(ZZ,0,A.ncols())
    for i in range(0,len(rows),16):
        H=H.stack(matrix(ZZ,rows[i:i+16])).hermite_form(algorithm='pari4',include_zero_rows=False)
        if H.nrows()==A.ncols() and i%128==0:
            B=scales*H.change_ring(QQ).inverse()
            # A subset of constraints gives a superlattice. If every basis
            # vector satisfies all original constraints, it is the full lattice.
            if (M*B).denominator()==1:
                return B
    assert H.nrows()==M.ncols()
    B=scales*H.change_ring(QQ).inverse()
    assert (M*B).denominator()==1
    return B


def torus():
    C=CartanMatrix('E7');highest=(0,0,0,0,0,0,1)
    weights={highest};pending=[highest]
    while pending:
        w=pending.pop()
        for i in range(7):
            v=tuple(w[j]-w[i]*C[i,j] for j in range(7))
            if v not in weights:weights.add(v);pending.append(v)
    roots=RootSystem('E7').ambient_space().simple_coroots()
    char=WeylCharacterRing('E7',style='coroots')(*highest)
    independent={tuple(w.inner_product(roots[i+1]) for i in range(7)):int(m)
                 for w,m in char.weight_multiplicities().items()}
    assert independent=={w:1 for w in weights} and len(weights)==56
    S=sorted(weights,reverse=True)
    M=matrix(ZZ,[[w[i] for w in S] for i in range(7)])
    assert M*M.transpose()==12*C
    U=matrix(ZZ,[[-1,-2,-3,-4,-3,-2,-1],[0,0,-1,0,0,0,0],
       [0,-1,-1,-1,-1,-1,-1],[0,-1,-1,-1,-1,-1,0],
       [0,-1,-1,-1,-1,0,0],[0,-1,-1,-1,0,0,0],[-1,-3,-4,-4,-3,-2,-1]])
    assert U.det()==1
    H=U*M;leading=[]
    for p in (0,2):
        pos=[next(i for i,c in enumerate(row) if c and (not p or c%p)) for row in H]
        coeff=[int(H[j,i]) for j,i in enumerate(pos)]
        assert coeff==([1]*6+[2] if not p else [1]*6+[-3])
        E=matrix(ZZ,[S[i] for i in pos]);assert E.det()!=0
        leading.append(dict(characteristic=p,positions=pos,coefficients=coeff,
                            exponents=[list(map(int,row)) for row in E],determinant=int(E.det())))
    L=LaurentPolynomialRing(QQ,['t'+str(i) for i in range(1,8)])
    ts=L.gens();chars=[prod(t**e for t,e in zip(ts,w)) for w in S]
    moments=M*vector(L,chars)
    theta=(12*C).inverse()*moments
    assert all(f(*([1]*7))==0 for f in theta)
    assert matrix(QQ,7,7,lambda i,j:theta[i].derivative(ts[j])(*([1]*7)))==identity_matrix(7)
    return C,theta,dict(cartan=[list(map(int,row)) for row in C],trace_index=12,
      ordered_weights=[list(map(int,w)) for w in S],unimodular_change=[list(map(int,row)) for row in U],
      transformed_coefficients=[list(map(int,row)) for row in H],leading=leading)


def compressed_inputs(b,C):
    R=PolynomialRing(QQ,['E'+str(i) for i in range(1,7)]+['u'])
    es=R.gens()[:6];u=R.gens()[6];E=es[0];z=3*E-6*u
    power=[R(6)]
    for k in range(1,13):
        value=sum((-1)**(i-1)*es[i-1]*power[k-i] for i in range(1,min(k,7)))
        if k<=6:value+=(-1)**(k-1)*k*es[k-1]
        power.append(R(value))
    X=[R(sum(binomial(k,j)*6**j*(-E)**(k-j)*power[j] for j in range(k+1))) for k in range(13)]
    psi={}
    for d in DEGREES[1:]:
        first=2*sum(binomial(d,j)*X[d-j]*z**j for j in range(0,d+1,2))
        pairs=QQ((-1)**d)/2*(sum(binomial(d,j)*X[j]*X[d-j] for j in range(d+1))-2**d*X[d])
        psi[d]=R(first+pairs)
    inputs=((3*u-E)/2,)+tuple(psi[d] for d in DEGREES[1:])
    A=PolynomialRing(QQ,['a'+str(i) for i in range(1,8)]);a=A.gens()
    chain=[0,2,3,4,5]
    vs=[a[6]+sum(a[j] for j in chain[i:]) for i in range(6)]
    uu=a[1]+sum(vs[:3])
    T=matrix(ZZ,[[v.monomial_coefficient(x) for x in a] for v in vs+[uu]])
    assert abs(T.det())==1
    # Every A5 simple reflection permutes v's and fixes u, over the exact lattice.
    reflection_checks=[]
    for i in chain:
        reflect=lambda f:f(*[a[j]-C[j,i]*a[i] for j in range(7)])
        perm=[vs.index(reflect(v)) for v in vs]
        assert reflect(uu)==uu
        reflection_checks.append(dict(node=i+1,permutation=perm))
    from itertools import combinations
    # Compare the defining linear forms first, then verify Newton's formula
    # in independent permutation variables. This avoids expanding a dense
    # change of basis in a degree-12 fundamental-weight polynomial.
    definition=json.loads((ROOT/'src/sage_borel/data/exceptional_generators.json').read_text())['conventions']['psi_e7']
    total=sum(vs)
    for v,coeff in zip(vs,definition['x_simple_root_coefficients']):
        assert 6*v-total==sum(c*x for c,x in zip(coeff,a))
    assert 3*total-6*uu==sum(c*x for c,x in zip(definition['z_simple_root_coefficients'],a))
    assert (3*uu-total)/2==(C.inverse()*vector(A,a))[6]
    V=PolynomialRing(QQ,['v'+str(i) for i in range(1,7)]+['u'])
    vv=V.gens()[:6];vu=V.gens()[6];ve=sum(vv)
    subst=[sum(prod(vv[i] for i in I) for I in combinations(range(6),k)) for k in range(1,7)]+[vu]
    xx=[6*v-ve for v in vv];zz=3*ve-6*vu
    for d,f in zip(DEGREES[1:],inputs[1:]):
        direct=sum((x+zz)**d+(x-zz)**d for x in xx)+sum((-xx[i]-xx[j])**d for i in range(6) for j in range(i+1,6))
        assert f(*subst)==direct,d
    return R,inputs,dict(root_basis_change=[list(map(int,row)) for row in T],
       determinant=int(T.det()),permutation_checks=reflection_checks,
       compressed_generators={n:str(f) for n,f in zip(NAMES,inputs)},coefficient_weights=[1,2,3,4,5,6,1])


def direct_laurent(b,theta):
    names=NAMES[:2];ds=DEGREES[:2]
    polys=[b.generator_polynomial(n) for n in names]
    vals=[f(*theta) for f in polys]
    powers=[(d,p) for d in range(4) for p in monomials(ds,d)]
    expansions=[prod(f**e for f,e in zip(vals,p)) for d,p in powers]
    dicts=[f.dict() for f in expansions];support=sorted(set().union(*(set(D) for D in dicts)))
    M=matrix(QQ,[[D.get(e,0) for D in dicts] for e in support]);rows=[]
    for d in range(4):
        for mode in ('homogeneous','filtered'):
            indices=[i for i,(degree,p) in enumerate(powers) if degree==d or (mode=='filtered' and degree<d)]
            A=M.matrix_from_columns(indices);B=integral_preimage(A)
            rows.append(dict(degree=d,mode=mode,monomial_exponents=[list(powers[i][1]) for i in indices],
                source_basis_matrix=[[str(v) for v in row] for row in B],
                completeness_certificate=certify_preimage(A,B)))
    return rows


def lattice_compare(H,oldH,form,old_factors):
    to_old=oldH.change_ring(QQ).inverse()*H
    to_new=H.change_ring(QQ).inverse()*oldH
    S=H.augment(oldH).transpose().hermite_form(include_zero_rows=False).transpose()
    in_sum=matrix(ZZ,S.change_ring(QQ).inverse()*H)
    old_in_sum=matrix(ZZ,S.change_ring(QQ).inverse()*oldH)
    K=H.augment(-oldH).right_kernel().basis_matrix().transpose()
    I=H*K[:H.ncols(),:]
    I=I.transpose().hermite_form(include_zero_rows=False).transpose()
    i_new=matrix(ZZ,H.change_ring(QQ).inverse()*I)
    i_old=matrix(ZZ,oldH.change_ring(QQ).inverse()*I)
    assert abs(I.det()*S.det())==abs(H.det()*oldH.det())
    return dict(form=form,old_contains_character=to_old.denominator()==1,
       character_contains_old=to_new.denominator()==1,old_smith_factors=old_factors,
       old_basis=[list(map(int,row)) for row in oldH],sum_basis=[list(map(int,row)) for row in S],
       intersection_basis=[list(map(int,row)) for row in I],
       sum_mod_character=list(map(int,in_sum.elementary_divisors())),
       sum_mod_old=list(map(int,old_in_sum.elementary_divisors())),
       old_mod_intersection=list(map(int,i_old.elementary_divisors())),
       character_mod_intersection=list(map(int,i_new.elementary_divisors())))


def main():
    start=perf_counter();C,theta,torus_cert=torus();b=exceptional_borel('E7',7)
    print('Torus certificate checked',round(perf_counter()-start,2),flush=True)
    R,inputs,compression=compressed_inputs(b,C)
    print('Integral compression checked',round(perf_counter()-start,2),flush=True)
    low=direct_laurent(b,theta)
    print('Direct low-degree Laurent check completed',round(perf_counter()-start,2),flush=True)
    result=dict(scope='Full E7/P7 image of Sym_Z(12Q)^(W_L), with Q the E7 root lattice.',
       lattice_algorithm='Exact column-content preconditioning and PARI4 Hermite batches of 16. Early return only when every original coefficient row is integral on the candidate basis. Independent gcd-one maximal-minor completeness certificate.',
       torus_certificate=torus_cert,compression_certificate=compression,direct_laurent_rows=low,degrees=[])
    output=ROOT/'experiments/e7-character-form-results.json'
    loc=FixedPointBorel('E7',7);pub=PublishedMinusculeRing('E7')
    vals=[loc.polynomial_values(b.generator_polynomial(n)) for n in NAMES]
    old=json.loads((ROOT/'experiments/exceptional-minuscule-results.json').read_text())
    oldcases=[c for c in old['cases'] if c['kind']=='E7']
    for d in range(28):
        powers=list(monomials(DEGREES,d))
        polys=[prod(f**e for f,e in zip(inputs,p)) for p in powers]
        dicts=[f.dict() for f in polys];exponents=sorted(set().union(*(set(D) for D in dicts)))
        coefficient_rows=sorted(set(tuple(D.get(e,0) for D in dicts) for e in exponents))
        M=matrix(QQ,coefficient_rows);B0=compressed_integral_preimage(M);cert=certify_preimage(M,B0);B=12**d*B0
        if d<=3:
            direct=next(r for r in low if r['degree']==d and r['mode']=='homogeneous')
            assert [list(p[:2]) for p in powers]==direct['monomial_exponents']
            transition=matrix(QQ,direct['source_basis_matrix']).inverse()*B
            assert transition.denominator()==1 and abs(transition.det())==1
        target=[w for w in loc.basis if w.length()==d];columns=[]
        for p in powers:
            values=[prod(vals[j][i]**e for j,e in enumerate(p)) for i in range(len(loc.basis))]
            im=loc.coefficients(values,d)
            rp=prod(pub.generator_expressions[n]**e for n,e in zip(NAMES,p))
            check={loc.W.from_reduced_word(w.reduced_word()):v for w,v in pub.expand(rp,d).items()}
            assert im==check
            columns.append([im.get(w,0) for w in target])
        im=matrix(ZZ,matrix(QQ,columns).transpose()*B)
        H=im.transpose().hermite_form(include_zero_rows=False).transpose()
        factors=tuple(im.elementary_divisors());assert factors==minor_smith(im)
        comparisons=[]
        for case in oldcases:
            row=case['degrees'][d]
            oldM=matrix(ZZ,[row['matrix'][row['schubert_words'].index(list(w.reduced_word()))] for w in target])
            oldH=oldM.transpose().hermite_form(include_zero_rows=False).transpose()
            comparisons.append(lattice_compare(H,oldH,case['form'],row['smith_factors']))
        result['degrees'].append(dict(degree=d,source_rank=len(powers),names=list(NAMES),
          monomial_exponents=[list(p) for p in powers],root_coefficient_rows=len(coefficient_rows),
          source_basis_matrix=[[str(v) for v in row] for row in B],completeness_certificate=cert,
          schubert_words=[list(w.reduced_word()) for w in target],image_matrix=[list(map(int,row)) for row in im],
          image_basis=[list(map(int,row)) for row in H],smith_factors=list(map(int,factors)),comparisons=comparisons))
        output.write_text(json.dumps(result,indent=2)+'\n')
        print('degree',d,'rank',len(powers),'rows',len(coefficient_rows),'Smith',factors,
              'contained',[c['old_contains_character'] for c in comparisons],'seconds',round(perf_counter()-start,2),flush=True)
    result.update(result='passed',seconds=round(perf_counter()-start,3))
    files=[Path(__file__).resolve(),ROOT/'experiments/e6_character_form.py',
      ROOT/'experiments/exceptional_localization.py',ROOT/'experiments/exceptional_published_rings.py',
      ROOT/'experiments/exceptional_minuscule_images.py',ROOT/'experiments/exceptional-minuscule-results.json',
      ROOT/'src/sage_borel/data/exceptional_generators.json',ROOT/'tests/test_exceptional.py']
    files+=sorted((ROOT/'src/sage_borel').glob('*.py'))
    result['source_hash_semantics']='SHA256 with CRLF normalized to LF.'
    result['source_sha256']={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in files}
    output.write_text(json.dumps(result,indent=2)+'\n');print('PASSED',result['seconds'],flush=True)


if __name__=='__main__':main()
