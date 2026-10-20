"""Literal Levi-character intersection through a specified polynomial degree.

Run with PYTHONPATH=src sage -python experiments/e6_character_form.py.
The full intersection is filtered, not assumed graded or finitely generated.
"""
from pathlib import Path
import json
from hashlib import sha256
from time import perf_counter
from sage.all import (QQ, ZZ, RootSystem, WeylCharacterRing, LaurentPolynomialRing,
                      PolynomialRing, matrix, vector, prod, lcm, identity_matrix)
from sage_borel import exceptional_borel
from exceptional_minuscule_images import monomials
from exceptional_localization import FixedPointBorel
from exceptional_published_rings import PublishedMinusculeRing

ROOT = Path(__file__).resolve().parents[1]


def springer():
    domain = RootSystem('E6').ambient_space()
    coroots = domain.simple_coroots()
    roots = domain.simple_roots()
    cartan = matrix(ZZ,6,6,lambda i,j:roots[i+1].inner_product(coroots[j+1]))
    highest = (1,0,0,0,0,0)
    weights, pending = {highest}, [highest]
    while pending:
        mu = pending.pop()
        for i in range(6):
            nu = tuple(mu[j]-mu[i]*cartan[i,j] for j in range(6))
            if nu not in weights:
                weights.add(nu); pending.append(nu)
    character = WeylCharacterRing('E6',style='coroots')(*highest)
    independent = {tuple(w.inner_product(coroots[i+1]) for i in range(6)):int(m)
                   for w,m in character.weight_multiplicities().items()}
    assert independent == {w:1 for w in weights}
    weights = sorted(weights)
    gram = matrix(QQ,6,6,lambda i,j:sum(w[i]*w[j] for w in weights))
    L = LaurentPolynomialRing(QQ,['t'+str(i) for i in range(1,7)])
    ts = L.gens()
    chars = {w:prod(t**e for t,e in zip(ts,w)) for w in weights}
    theta = gram.inverse()*vector(L,[sum(w[i]*chars[w] for w in weights) for i in range(6)])
    assert all(f(*([1]*6))==0 for f in theta)
    assert matrix(QQ,6,6,lambda i,j:theta[i].derivative(ts[j])(*([1]*6)))==identity_matrix(6)
    return L,theta,cartan


def integral_preimage(M):
    """Columns B span exactly {q in Q^n : M*q is integral}."""
    denominator = M.denominator()
    H = matrix(ZZ,denominator*M).hermite_form(include_zero_rows=False)
    assert H.nrows()==H.ncols()==M.ncols()
    B = denominator*H.change_ring(QQ).inverse()
    assert (M*B).denominator()==1
    return B


def certify_preimage(M,B):
    """Independently certify completeness by a gcd-one set of maximal minors.

    Modular pivots select small nonsingular submatrices; exact determinants
    certify primitive columns. No Smith or Hermite algorithm is used here.
    """
    from sage.all import GF, gcd
    N=matrix(ZZ,M*B)
    minors=[];divisor=ZZ.zero()
    for p in (2,3,5,7,11,13,17,19,23,29):
        selected=N.change_ring(GF(p)).transpose().pivots()
        assert len(selected)==N.ncols()
        determinant=abs(N.matrix_from_rows(selected).det())
        assert determinant
        minors.append(dict(rows=list(map(int,selected)),determinant=str(determinant),selection_prime=p))
        divisor=gcd(divisor,determinant)
        if divisor==1:break
    assert divisor==1
    return dict(maximal_minors=minors,minor_gcd=1,rank=int(N.ncols()))


def low_degree(L,theta,b):
    names=('y1','e1','e2')
    degrees=(1,2,4)
    P=PolynomialRing(QQ,names)
    # Build from the six linear coordinates, avoiding a large polynomial substitution.
    coords=b.generator_info('y1').coordinates
    ys=[sum(QQ(c)*t for c,t in zip(coords[str(i)],theta)) for i in range(1,7)]
    squares=[v*v for v in ys[1:]]
    inputs=(ys[0],sum(squares),sum(squares[i]*squares[j] for i in range(5) for j in range(i+1,5)))
    polynomials=tuple(b.generator_polynomial(n) for n in names)
    powers=[(d,p) for d in range(5) for p in monomials(degrees,d)]
    labels=[prod(x**e for x,e in zip(P.gens(),p)) for d,p in powers]
    expansions=[prod(g**e for g,e in zip(inputs,p)) for d,p in powers]
    dictionaries=[{tuple(e):c for e,c in f.dict().items()} for f in expansions]
    support=sorted(set().union(*(set(d) for d in dictionaries)))
    M=matrix(QQ,[[D.get(w,0) for D in dictionaries] for w in support])
    print('Direct Laurent coefficient matrix',M.dimensions(),flush=True)
    report=[]
    for d in range(5):
        for mode in ('homogeneous','filtered'):
            indices=[i for i,(degree,p) in enumerate(powers) if degree==d or (mode=='filtered' and degree<d)]
            A=M.matrix_from_columns(indices)
            B=integral_preimage(A)
            certificate=certify_preimage(A,B)
            source=[sum(B[i,j]*labels[k] for i,k in enumerate(indices)) for j in range(B.ncols())]
            target_labels=[(e,w) for e in range(d+1) for w in b.basis(e)] if mode=='filtered' else [(d,w) for w in b.basis(d)]
            images=[]
            for k in indices:
                e,p=powers[k]
                f=prod(g**a for g,a in zip(polynomials,p))
                im=b.expand(f)
                images.append([im.get(w,0) if e==degree else 0 for degree,w in target_labels])
            C=matrix(QQ,images).transpose()*B
            if C.denominator()==1:
                smith=list(map(int,matrix(ZZ,C).elementary_divisors()))
            else:
                smith=None
            report.append(dict(degree=d,mode=mode,source_basis=list(map(str,source)),
                               source_basis_matrix=[[str(v) for v in row] for row in B.rows()],
                               monomials=[str(labels[i]) for i in indices],
                               monomial_exponents=[list(powers[i][1]) for i in indices],
                               image_matrix=[[str(v) for v in row] for row in C.rows()],
                               smith_factors=smith,completeness_certificate=certificate))
    return report


def main():
    start=perf_counter()
    L,theta,cartan=springer()
    b=exceptional_borel('E6',1)
    # An integer change of the six trace-moment coordinates.
    U=matrix(ZZ,[[0,-2,-1,-3,-2,-1],[1,-1,1,0,0,0],
                 [0,-2,-1,-2,-2,-1],[0,-2,-1,-2,-1,-1],
                 [0,-2,-1,-2,-1,0],[0,-3,-1,-3,-2,-1]])
    assert U.det()==-1
    moments=6*cartan*theta
    assert all(c.denominator()==1 for f in moments for c in f.dict().values())
    transformed=U*moments
    support=sorted(set().union(*(set(f.dict()) for f in moments)),reverse=True)
    coefficient_rows=[[int(f.dict().get(e,0)) for e in support] for f in transformed]
    leading=[]
    for prime in (0,3):
        positions=[next(i for i,c in enumerate(row) if c and (not prime or c%prime)) for row in coefficient_rows]
        exponents=[list(map(int,support[i])) for i in positions]
        coefficients=[coefficient_rows[j][i] for j,i in enumerate(positions)]
        determinant=matrix(ZZ,exponents).det()
        assert determinant!=0
        if not prime:
            assert coefficients==[1,1,1,1,1,3]
        else:
            assert coefficients==[1,1,1,1,1,-2]
        leading.append(dict(characteristic=prime,positions=positions,exponents=exponents,
                            coefficients=coefficients,exponent_determinant=int(determinant)))
    low=low_degree(L,theta,b)
    names=('y1','e1','e2','euler','e3','e4')
    degrees=(1,2,4,5,6,8)
    P=PolynomialRing(QQ,names)
    A=PolynomialRing(QQ,['a'+str(i) for i in range(1,7)])
    weights=cartan.inverse()*vector(A,A.gens())
    # Work first with the root lattice Q; multiply a homogeneous basis by 6^d.
    coordinates=b.generator_info('y1').coordinates
    ys=[sum(QQ(c)*t for c,t in zip(coordinates[str(i)],weights)) for i in range(1,7)]
    squares=[y*y for y in ys[1:]]
    from itertools import combinations
    elementary=lambda k:sum(prod(squares[i] for i in I) for I in combinations(range(5),k))
    roots=(ys[0],elementary(1),elementary(2),prod(ys[1:]),elementary(3),elementary(4))
    loc=FixedPointBorel('E6',1)
    pub=PublishedMinusculeRing('E6')
    values=[loc.polynomial_values(b.generator_polynomial(n)) for n in names]
    old=json.loads((ROOT/'experiments/exceptional-minuscule-results.json').read_text())
    oldcases=[c for c in old['cases'] if c['kind']=='E6']
    result=dict(scope='All degrees 0..16 for the literal Levi-character intersection, using the proved full-torus integral form Sym_Z(6Q).',
                torus_certificate=dict(cartan=[list(map(int,row)) for row in cartan],
                  unimodular_change=[list(map(int,row)) for row in U],leading=leading,
                  ordered_weights=[list(map(int,e)) for e in support],coefficient_rows=coefficient_rows),
                direct_laurent_rows=low,degrees=[])
    for d in range(17):
        powers=list(monomials(degrees,d))
        labels=[prod(x**e for x,e in zip(P.gens(),p)) for p in powers]
        polys=[prod(f**e for f,e in zip(roots,p)) for p in powers]
        dicts=[f.dict() for f in polys]
        exponents=sorted(set().union(*(set(D) for D in dicts)))
        # Repeated coefficient rows impose identical integrality constraints.
        rows=sorted(set(tuple(D.get(e,0) for D in dicts) for e in exponents))
        M=matrix(QQ,rows)
        root_basis=integral_preimage(M)
        certificate=certify_preimage(M,root_basis)
        B=6**d*root_basis
        if d<=4:
            direct=next(r for r in low if r['degree']==d and r['mode']=='homogeneous')
            assert direct['monomials']==list(map(str,labels))
            B0=matrix(QQ,direct['source_basis_matrix'])
            transition=B0.inverse()*B
            assert transition.denominator()==1 and abs(transition.det())==1
        target=[w for w in loc.basis if w.length()==d]
        columns=[]
        for p in powers:
            vals=[prod(values[j][i]**e for j,e in enumerate(p)) for i in range(len(loc.basis))]
            im=loc.coefficients(vals,d)
            ringpoly=prod(pub.generator_expressions[n]**e for n,e in zip(names,p))
            check={loc.W.from_reduced_word(w.reduced_word()):v for w,v in pub.expand(ringpoly,d).items()}
            assert im==check
            columns.append([im.get(w,0) for w in target])
        image=matrix(ZZ,matrix(QQ,columns).transpose()*B)
        H=image.transpose().hermite_form(include_zero_rows=False).transpose()
        # Smith factors independently via all minors of a <=3 by 3 basis.
        from exceptional_minuscule_images import minor_smith
        factors=tuple(image.elementary_divisors())
        assert factors==minor_smith(image)
        comparisons=[]
        for c in oldcases:
            oldrow=c['degrees'][d]
            oldlabels=[tuple(w) for w in oldrow['schubert_words']]
            reordered=matrix(ZZ,[oldrow['matrix'][oldlabels.index(tuple(w.reduced_word()))] for w in target])
            oldH=reordered.transpose().hermite_form(include_zero_rows=False).transpose()
            coordinates=oldH.change_ring(QQ).inverse()*H
            assert coordinates.denominator()==1
            quotient=list(map(int,matrix(ZZ,coordinates).elementary_divisors()))
            comparisons.append(dict(form=c['form'],old_image_contains_character_image=True,
                                    old_smith_factors=oldrow['smith_factors'],quotient_factors=quotient,
                                    containment_matrix=[list(map(int,row)) for row in coordinates]))
        sources=[str(sum(B[i,j]*labels[i] for i in range(len(labels)))) for j in range(B.ncols())]
        result['degrees'].append(dict(degree=d,source_rank=len(powers),monomials=list(map(str,labels)),monomial_exponents=[list(p) for p in powers],
          root_coefficient_rows=len(rows),source_basis=sources,
          source_basis_matrix=[[str(v) for v in row] for row in B],
          completeness_certificate=certificate,schubert_words=[list(w.reduced_word()) for w in target],
          image_matrix=[list(map(int,row)) for row in image],image_basis=[list(map(int,row)) for row in H],
          smith_factors=list(map(int,factors)),comparisons=comparisons))
        print('Degree',d,'source rank',len(powers),'coefficient rows',len(rows),'Smith',factors,
              'elapsed',round(perf_counter()-start,2),flush=True)
        (ROOT/'experiments/e6-character-form-results.json').write_text(json.dumps(result,indent=2)+'\n')
    result.update(result='passed',seconds=round(perf_counter()-start,3))
    files=[Path(__file__).resolve(),ROOT/'experiments/exceptional-minuscule-results.json',
           ROOT/'experiments/exceptional_localization.py',ROOT/'experiments/exceptional_published_rings.py',
           ROOT/'experiments/exceptional_minuscule_images.py',ROOT/'src/sage_borel/data/exceptional_generators.json']
    files+=sorted((ROOT/'src/sage_borel').glob('*.py'))
    result['source_hash_semantics']='SHA256 with CRLF normalized to LF.'
    result['source_sha256']={str(p.relative_to(ROOT)):sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in files}
    (ROOT/'experiments/e6-character-form-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASSED',result['seconds'],'seconds',flush=True)


if __name__=='__main__':
    main()
