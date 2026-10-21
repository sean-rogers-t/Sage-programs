"""Exact 3-local root/weight comparison on the Cayley plane and refinements.

Sage BGG identifies the six D5 invariant images. Graded relation matrices
then provide two-sided local lattice certificates and multiplication
closure, rather than relying on agreement with a previous Smith table.
"""
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from time import perf_counter
import json

from sage.all import QQ, ZZ, GF, CartanMatrix, PolynomialRing, WeylGroup, matrix, vector, diagonal_matrix
from exceptional_published_rings import PublishedMinusculeRing

ROOT = Path(__file__).resolve().parents[1]


def powers(weights, total):
    if not weights:
        if total == 0:
            yield ()
        return
    for i in range(total // weights[0] + 1):
        for rest in powers(weights[1:], total - i*weights[0]):
            yield (i,) + rest


def entries(A):
    return [[str(x) for x in row] for row in A]


def local(A):
    return all(x.denominator() % 3 for x in A.list())


def main():
    start = perf_counter()
    pub = PublishedMinusculeRing('E6')
    R = pub.R
    h, s = R.gens()
    names = ('y1', 'e1', 'e2', 'euler', 'e3', 'e4')
    degrees = (1, 2, 4, 5, 6, 8)
    generators = [pub.generator_expressions[name] for name in names]
    expected = [3*h, -3*h**2, 42*h**4-96*s, -21*h**5+48*h*s,
                282*h**6-672*h**2*s, -1371*h**8+3552*h**4*s-768*s**2]
    assert generators == expected
    print('Six invariant images reconstructed by BGG', flush=True)
    coordinates, bases, degrees_data = {}, {}, []
    for d in range(17):
        mon = [(d-4*j, j) for j in range(d//4+1)]
        chosen = [e for e in mon if e[0] < 9 and e[1] < 3]
        B = [h**i*s**j for i,j in chosen]
        rel = [f*h**i*s**j for f,k in zip(pub.relations, (9,12))
               for i,j in powers((1,4), d-k)]
        columns = B + rel
        A = matrix(QQ, [[f.monomial_coefficient(h**i*s**j) for f in columns] for i,j in mon])
        assert A.is_square() and A.det().valuation(3) == 0
        transform = A.inverse()[:len(B),:]
        def coords(f, mon=mon, transform=transform):
            return transform*vector(QQ, [f.monomial_coefficient(h**i*s**j) for i,j in mon])
        coordinates[d] = coords
        # N has basis h^i, 3 h^i s, 3 h^i s^2 over Z_(3).
        nb = [f*(1 if j == 0 else 3) for f,(i,j) in zip(B,chosen)]
        scale = diagonal_matrix(QQ, [1 if j == 0 else 3 for i,j in chosen])
        formal = list(powers(degrees,d))
        polys = [R.prod(g**e for g,e in zip(generators, exps)) for exps in formal]
        T = matrix(QQ, [scale.inverse()*coords(f) for f in polys]).transpose()
        assert local(T)
        D = T.denominator()
        Aint = matrix(ZZ, D*T)
        smith, U, V = Aint.smith_form()
        exponents = [int(smith[i,i].valuation(3)) for i in range(len(B))]
        H = U.change_ring(QQ).inverse()*diagonal_matrix(QQ, [3**a for a in exponents])
        # Choose actual source monomials as a local module basis. A minor
        # attaining the Smith determinant valuation always exists over a DVR.
        selected = next(indices for indices in combinations(range(T.ncols()),len(B))
                        if T.matrix_from_columns(indices).det() and
                        T.matrix_from_columns(indices).det().valuation(3)==sum(exponents))
        H = T.matrix_from_columns(selected)
        witness = matrix(QQ,T.ncols(),len(B))
        for j,i in enumerate(selected):
            witness[i,j]=1
        assert local(witness) and T*witness == H
        assert local(H.inverse()*T)
        mb = [sum(c*f for c,f in zip(H.column(j),nb)) for j in range(len(B))]
        bases[d] = (nb, mb, scale, H)
        row = dict(degree=d, monomial_basis=[list(e) for e in chosen],
                   relation_basis_determinant=str(A.det()), root_monomials=[list(e) for e in formal],
                   image_columns_in_weight_basis=entries(T), root_basis_in_weight_basis=entries(H),
                   root_basis_witness=entries(witness), smith_exponents=exponents,
                   selected_root_monomials=list(selected),
                   nonunit_factors=[3**a for a in exponents if a],
                   root_basis=[str(f) for f in mb])
        degrees_data.append(row)
    # Finite closed module certificate: basis vectors are in the generated
    # ring (two-sided witnesses above), and multiplication by every algebra
    # generator preserves these degree components. N is checked separately.
    closure = 0
    for d,(nb,mb,scale,H) in bases.items():
        for fs, gs, ws, is_root in ((mb,generators,degrees,True),
                                   (nb,[h,3*s,3*s*s],(1,4,8),False)):
            for f in fs:
                for g,k in zip(gs,ws):
                    if d+k <= 16:
                        _,_,ns,nh = bases[d+k]
                        v = ns.inverse()*coordinates[d+k](f*g)
                        if is_root:
                            v = nh.inverse()*v
                        assert all(x.denominator()%3 for x in v)
                        closure += 1
    # Every D5 parabolic refinement, and its E6 diagram-involuted partner.
    W = WeylGroup('D5',prefix='s')
    elements = list(W)
    refinements = []
    polynomials = {}
    for mask in range(32):
        retained = [i+1 for i in range(5) if mask & (1<<i)]
        counts = [0]*21
        for w in elements:
            if not any(w.has_descent(i,side='right') for i in retained):
                counts[w.length()] += 1
        while counts[-1] == 0:
            counts.pop()
        polynomials[str(mask)] = counts
        factors = {str(a): [sum(counts[k]*sum(x==a for x in degrees_data[d-k]['smith_exponents'])
                                         for k in range(min(d,len(counts)-1)+1)
                                         if d-k <= 16)
                                  for d in range(17+len(counts)-1)] for a in (1,2,3)}
        refinements.append(dict(d5_retained=retained, fiber_betti=counts,
                                e6_retained=[(6,5,4,3,2)[i-1] for i in retained],
                                multiplicities=factors))
    # Match prior finite tables after removing the representation scaling.
    old = json.loads((ROOT/'experiments/minuscule-character-api-results.json').read_text())
    adj = next(c for c in old['cases'] if c['kind']=='E6' and c['representation']=='adjoint')
    for row, prior in zip(degrees_data,adj['degrees']):
        assert row['nonunit_factors'] == prior['comparison']['sum_mod_minimal']['torsion_factors']
    # A common fiber basis also permits deconvolution from full flags
    # whenever the adjoint Levi has no 3-primary isogeny obstruction.
    P=PolynomialRing(ZZ,'q'); q=P.gen()
    def length_poly(ds):return prod(sum(q**j for j in range(k)) for k in ds)
    from sage.all import prod
    D5=length_poly((2,4,5,6,8))
    full=[P(sum(q**r['degree']*r['smith_exponents'].count(a) for r in degrees_data))*D5 for a in (1,2,3)]
    C=matrix(ZZ,CartanMatrix('E6'))
    all_parabolics=[]
    for mask in range(63):
        I=[i for i in range(6) if mask&(1<<i)]
        remaining=set(I); ds=[]
        while remaining:
            comp={remaining.pop()}; pending=list(comp)
            for i in pending:
                near={j for j in remaining if C[i,j]}
                remaining-=near;comp|=near;pending+=list(near)
            n=len(comp)
            branch=max(sum(C[i,j]==-1 for j in comp) for i in comp)
            ds+=([2,4,4,6] if n==4 and branch==3 else
                 [2,4,5,6,8] if n==5 and branch==3 else list(range(2,n+2)))
        fiber=length_poly(ds)
        good=C.matrix_from_columns(I).change_ring(GF(3)).rank()==len(I)
        quotients=[f.quo_rem(fiber) for f in full]
        assert good==all(rem==0 for quo,rem in quotients)
        if good:
            assert all(c>=0 for quo,rem in quotients for c in quo.list())
        all_parabolics.append(dict(e6_retained=[i+1 for i in I],
            common_fiber_basis=good, levi_length_polynomial=list(map(int,fiber.list())),
            multiplicities={str(a):list(map(int,quo.list())) for a,(quo,rem) in zip((1,2,3),quotients)} if good else {},
            nonzero_remainders={str(a):str(rem) for a,(quo,rem) in zip((1,2,3),quotients)} if not good else {}))
    assert sum(r['common_fiber_basis'] for r in all_parabolics)==60
    files = [Path(__file__).resolve(),ROOT/'experiments/exceptional_published_rings.py',
             ROOT/'experiments/minuscule-character-api-results.json']
    result = dict(result='passed', coefficient_ring='Z_(3)',
        scope='E6/P1 base, 48 minuscule refinements by transfer, and 60 proper standard parabolics by good-Levi deconvolution. Formula outputs are not 60 direct BGG recomputations. Three proper parabolics remain outside the theorem.',
        invariant_images={n:str(f) for n,f in zip(names,generators)}, degrees=degrees_data,
        generator_schubert={n:[dict(word=list(w.reduced_word()),coefficient=str(c))
                              for w,c in pub.borel.generator_image(n).items()] for n in names},
        schubert_labels={str(d):[list(w.reduced_word()) for w in pub.labels[d]] for d in range(17)},
        e6_square_anchor=pub.e6_square_anchor,
        module_closure_checks=closure, d5_refinements=refinements,
        proper_parabolics=all_parabolics,
        seconds=round(perf_counter()-start,3),
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in files})
    (ROOT/'experiments/e6-root-weight-comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASSED',result['seconds'],'closure',closure,
          [(r['degree'],r['nonunit_factors']) for r in degrees_data],flush=True)


if __name__ == '__main__':
    main()
