"""Complete E6/P1 and E7/P7 image matrices with independent ring checks.

Run with PYTHONPATH=src sage -python experiments/exceptional_minuscule_images.py.
All catalogue monomials are retained, including redundant generators. No
high-degree fundamental-weight polynomial products are expanded.
"""
from datetime import datetime, timezone
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys
from time import perf_counter

from sage.all import QQ, ZZ, PolynomialRing, gcd, lcm, matrix, prod
from sage.env import SAGE_VERSION
from sage_borel import exceptional_borel
from exceptional_localization import FixedPointBorel
from exceptional_published_rings import PublishedMinusculeRing

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from tests.test_exceptional import independent_polynomial, FIXTURE_PATH


def monomials(degrees, total):
    if not degrees:
        if total == 0:
            yield ()
        return
    for exponent in range(total//degrees[0]+1):
        for tail in monomials(degrees[1:],total-exponent*degrees[0]):
            yield (exponent,)+tail


def evaluate_monomial(values, powers):
    return [prod((values[j][i]**e for j,e in enumerate(powers) if e), QQ.one())
            for i in range(len(values[0]))]


def minor_smith(M):
    """A different integer reduction: Hermite lattice basis, then minor gcds."""
    H = M.transpose().hermite_form(include_zero_rows=False).transpose()
    previous, answer = ZZ.one(), []
    for size in range(1,H.rank()+1):
        divisor = gcd([H.matrix_from_rows_and_columns(rows,cols).det()
                       for rows in combinations(range(H.nrows()),size)
                       for cols in combinations(range(H.ncols()),size)])
        answer.append(abs(divisor//previous))
        previous = abs(divisor)
    return tuple(answer)


def regression():
    source = ROOT/'experiments/exceptional-integral-results.json'
    old = json.loads(source.read_text())
    assert old['result'] == 'passed'
    contexts, columns = {}, 0
    for case in old['cases']:
        key = case['kind'],case['parabolic']
        if key not in contexts:
            contexts[key] = (FixedPointBorel(*key),exceptional_borel(*key))
        loc,b = contexts[key]
        R = PolynomialRing(QQ,len(case['retained_generators']),names=case['retained_generators'])
        values = [loc.polynomial_values(b.ring(p)) for p in case['retained_polynomials']]
        for row in case['degrees']:
            for j,m in enumerate(row['monomials']):
                actual = loc.coefficients(evaluate_monomial(values,tuple(R(m).exponents()[0])),row['degree'])
                expected = {loc.W.from_reduced_word(word):QQ(row['matrix'][i][j])
                            for i,word in enumerate(row['schubert_words']) if row['matrix'][i][j]}
                assert actual == expected, (key,row['degree'],m)
                columns += 1
    print('G2/F4 localization regression:',columns,'columns',flush=True)
    return columns


def main():
    started = perf_counter()
    output = ROOT/'experiments/exceptional-minuscule-results.json'
    result = {'sage_version': SAGE_VERSION, 'generated_utc': datetime.now(timezone.utc).isoformat(),
              'scope': 'Every degree of E6/P1 and E7/P7 for the catalogue and individually primitive-image source forms. No character-intersection or novelty claim.',
              'methods': ['Fixed-point Billey interpolation at two regular rational points',
                          'Duan--Zhao integral presentations plus ordinary root/coroot Chevalley multiplication; BGG supplies catalogue images and one E6 square anchor',
                          'Smith elementary divisors checked by Hermite reduction and determinantal divisors'],
              'g2_f4_regression_columns': regression(), 'cases': []}
    fixture = json.loads(FIXTURE_PATH.read_text())['rows']
    columns, unscaled_columns, inputs, ring_reports = 0, 0, 0, []
    for kind,node in (('E6',1),('E7',7)):
        start_case = perf_counter()
        b = exceptional_borel(kind,node)
        loc = FixedPointBorel(kind,node)
        other = FixedPointBorel(kind,node,point_weights=(2,3,5,7,11,13,17)[:b.rank])
        published = PublishedMinusculeRing(kind)
        assert loc.top == published.top
        assert set(loc.basis) == {loc.W.from_reduced_word(w.reduced_word())
                                  for d in range(loc.top+1) for w in b.basis(d)}
        names = tuple(g.name for g in b.generators())
        degrees = tuple(g.polynomial_degree for g in b.generators())
        R = PolynomialRing(QQ,len(names),names=names)
        polynomials = tuple(b.generator_polynomial(name) for name in names)
        audit, cleared, primitive = [], [], []
        for g,p in zip(b.generators(),polynomials):
            source_row = next(row for row in fixture
                              if row['kind']==kind and row['parabolic']==node and row['name']==g.name)
            assert independent_polynomial(source_row,b.ring) == p
            image = b.expand(p)
            assert loc.expand(p) == {loc.W.from_reduced_word(w.reduced_word()):c for w,c in image.items()}
            denominator = lcm([c.denominator() for c in image.values()])
            numerator = gcd([c*denominator for c in image.values()])
            cleared.append(QQ(denominator))
            primitive.append(QQ(denominator)/numerator if numerator else QQ.one())
            audit.append({'name':g.name,'degree':g.polynomial_degree,'image':image.to_terms(),
                          'least_integer_image_multiplier':str(denominator),
                          'primitive_image_multiplier':str(primitive[-1]),
                          'published_ring_expression':str(published.generator_expressions[g.name])})
            inputs += 1
        values = [loc.polynomial_values(p) for p in polynomials]
        other_values = [other.polynomial_values(p) for p in polynomials]
        ring_inputs = tuple(published.generator_expressions[name] for name in names)
        @lru_cache(None)
        def ring_monomial(powers):
            if not any(powers):
                return published.R.one()
            j = next(j for j,e in enumerate(powers) if e)
            lower = list(powers)
            lower[j] -= 1
            return (ring_inputs[j]*ring_monomial(tuple(lower))).reduce(published.groebner)
        cases = []
        for form,scales in (('image_denominators_cleared',cleared),('individually_primitive_images',primitive)):
            cases.append({'kind':kind,'parabolic':node,'form':form,'complex_dimension':loc.top,
                          'fixed_point_count':len(loc.basis),'generators':audit,
                          'generator_names':list(names),'generator_degrees':list(degrees),
                          'multipliers':dict(zip(names,map(str,scales))),'degrees':[]})
        for d in range(loc.top+1):
            powers = list(monomials(degrees,d))
            labels = [w for w in loc.basis if w.length()==d]
            raw_columns, independent_columns = [], []
            for exponents in powers:
                primary = loc.coefficients(evaluate_monomial(values,exponents),d)
                second_point = other.coefficients(evaluate_monomial(other_values,exponents),d)
                assert primary == {loc.W.from_reduced_word(w.reduced_word()):c for w,c in second_point.items()}
                independent = published.expand(ring_monomial(exponents),d)
                independent = {loc.W.from_reduced_word(w.reduced_word()):c for w,c in independent.items()}
                assert primary == independent, (kind,d,exponents)
                independent_columns.append([independent.get(w,0) for w in labels])
                unscaled_columns += 1
                raw_columns.append(vector_column := [primary.get(w,0) for w in labels])
                assert all(QQ(c).denominator()==1 for c in vector_column)
            for case,scales in zip(cases,(cleared,primitive)):
                scalars = [prod(s**e for s,e in zip(scales,p)) for p in powers]
                M = matrix(ZZ,[[raw_columns[j][i]*scalars[j] for j in range(len(powers))]
                               for i in range(len(labels))])
                # Rescaling a homogeneous source monomial multiplies its
                # column by the product of its explicit generator scales.
                assert M == matrix(ZZ,[[independent_columns[j][i]*scalars[j] for j in range(len(powers))]
                                       for i in range(len(labels))])
                assert M.rank() == len(labels)
                factors = tuple(M.elementary_divisors())
                assert factors == minor_smith(M)
                row = {'degree':d,'rank':len(labels),'ambient_dimension':len(labels),
                       'smith_factors':[int(v) for v in factors],
                       'saturation_index':int(prod(factors)),
                       'monomials':[str(prod(x**e for x,e in zip(R.gens(),p))) for p in powers],
                       'exponents':[list(p) for p in powers],
                       'schubert_words':[list(w.reduced_word()) for w in labels],
                       'matrix':[[int(v) for v in row] for row in M.rows()]}
                case['degrees'].append(row)
                columns += len(powers)
                print(kind,case['form'],d,factors,flush=True)
        # The independent ring has no component above the geometric top.
        assert all(not ring_monomial(p) for p in monomials(degrees,loc.top+1))
        ring_reports.append({'kind':kind,'relations':[str(p) for p in published.relations],
                             'e6_square_anchor':getattr(published,'e6_square_anchor',None),
                             'seconds':round(perf_counter()-start_case,3)})
        result['cases'] += cases
        output.write_text(json.dumps(result,indent=2)+'\n')
        ring_monomial.cache_clear()
        published.oracle.covers.cache_clear()
    result.update(result='passed',degree_rows=sum(len(c['degrees']) for c in result['cases']),
                  independently_checked_columns=columns,independently_reconstructed_inputs=inputs,
                  direct_unscaled_column_comparisons=unscaled_columns,
                  normalization_check='Each unscaled column is independently computed; both source forms then use exact rational column scaling, with integrality and equality checked for every resulting column.',
                  published_ring_checks=ring_reports,seconds=round(perf_counter()-started,3))
    paths = [Path(__file__).resolve(),ROOT/'experiments/exceptional_localization.py',
             ROOT/'experiments/exceptional_published_rings.py',ROOT/'experiments/exceptional-integral-results.json',
             ROOT/'tests/test_exceptional.py',FIXTURE_PATH,ROOT/'src/sage_borel/data/exceptional_generators.json']
    paths += sorted((ROOT/'src/sage_borel').glob('*.py'))
    result['source_hash_semantics'] = 'SHA256 of UTF-8 file bytes with CRLF normalized to LF, matching Git text blobs.'
    result['source_sha256'] = {str(p.relative_to(ROOT)):sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths}
    output.write_text(json.dumps(result,indent=2)+'\n')
    print('PASSED',result['degree_rows'],'rows;',columns,'independently checked columns;',result['seconds'],'seconds',flush=True)


if __name__ == '__main__':
    main()
