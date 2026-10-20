"""Standard-library-only audit of E7 character form and all image comparisons."""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from time import perf_counter
from math import gcd
from itertools import combinations
import json
import sys
from audit_e6_character_form import (determinant,inverse,combine,multiply,
    modular_rank,matmul,smith_square,check_image_basis)

ROOT=Path(__file__).resolve().parents[1]
C=[[2,0,-1,0,0,0,0],[0,2,0,-1,0,0,0],[-1,0,2,-1,0,0,0],
   [0,-1,-1,2,-1,0,0],[0,0,0,-1,2,-1,0],[0,0,0,0,-1,2,-1],[0,0,0,0,0,-1,2]]
X=[[5,0,4,3,2,1,0],[-1,0,4,3,2,1,0],[-1,0,-2,3,2,1,0],
   [-1,0,-2,-3,2,1,0],[-1,0,-2,-3,-4,1,0],[-1,0,-2,-3,-4,-5,0]]
Z=[-3,-6,-6,-9,-6,-3,0]


def integral_matrix(A):
    assert all(F(x).denominator==1 for row in A for x in row)


def contained(A,B):
    return all(F(x).denominator==1 for row in matmul(inverse(B),A) for x in row)


def main():
    start=perf_counter();path=ROOT/'experiments/e7-character-form-results.json'
    data=json.loads(path.read_text());complete=data.get('result')=='passed'
    assert complete or sys.argv[1:]==['--prefix'], 'Full calculation must pass before a full audit.'
    assert [row['degree'] for row in data['degrees']]==list(range(len(data['degrees'])))
    if complete:
        assert len(data['degrees'])==28
        assert sum(row['source_rank'] for row in data['degrees'])==1967
        assert sum(len(row['schubert_words']) for row in data['degrees'])==56
    weights={(0,0,0,0,0,0,1)};pending=list(weights)
    while pending:
        w=pending.pop()
        for i in range(7):
            v=tuple(w[j]-w[i]*C[i][j] for j in range(7))
            if v not in weights:weights.add(v);pending.append(v)
    assert len(weights)==56
    gram=[[sum(w[i]*w[j] for w in weights) for j in range(7)] for i in range(7)]
    assert gram==[[12*c for c in row] for row in C]
    moments=[{w:w[i] for w in weights if w[i]} for i in range(7)]
    cert=data['torus_certificate'];U=cert['unimodular_change'];assert determinant(U)==1
    transformed=[combine(moments,row) for row in U];ordered=sorted(weights,reverse=True)
    assert [list(w) for w in ordered]==cert['ordered_weights']
    assert [[p.get(w,0) for w in ordered] for p in transformed]==cert['transformed_coefficients']
    for case in cert['leading']:
        p=case['characteristic']
        leading=[max(w for w,c in f.items() if c and (not p or c%p)) for f in transformed]
        assert [list(w) for w in leading]==case['exponents']
        assert [f[w] for f,w in zip(transformed,leading)]==case['coefficients']
        assert determinant(leading)==case['determinant']!=0
    for i in range(7):
        for j in range(7):
            reflected={tuple(w[k]-w[i]*C[i][k] for k in range(7)):c for w,c in moments[j].items()}
            assert reflected==combine([moments[j],moments[i]],[1,-C[j][i]])
    # Check the integral A5 permutation model independently in integer root coordinates.
    compression=data['compression_certificate'];T=compression['root_basis_change']
    assert abs(determinant(T))==1
    for check in compression['permutation_checks']:
        i=check['node']-1
        def reflect(v):
            out=v[:];out[i]-=sum(v[j]*C[j][i] for j in range(7));return out
        assert [T[:6].index(reflect(v)) for v in T[:6]]==check['permutation']
        assert reflect(T[6])==T[6]
    E=[sum(row[j] for row in T[:6]) for j in range(7)]
    assert [[6*v[j]-E[j] for j in range(7)] for v in T[:6]]==X
    assert [3*E[j]-6*T[6][j] for j in range(7)]==Z
    assert [F(3*T[6][j]-E[j],2) for j in range(7)]==inverse(C)[6]
    # Independent symbolic certificate for the degree-five character witness.
    variables=[{tuple(int(i==j) for j in range(7)):1} for i in range(7)]
    def power(f,n):
        out={(0,)*7:1}
        for _ in range(n):out=multiply(out,f)
        return out
    elementary=[{tuple(int(j in I) for j in range(7)):1 for I in combinations(range(6),k)} for k in range(1,7)]
    E1=elementary[0];uu=variables[6]
    xx=[combine([v,E1],[6,-1]) for v in variables[:6]];zz=combine([E1,uu],[3,-6])
    linear=[combine([x,zz],[1,s]) for x in xx for s in (1,-1)]
    linear += [combine([xx[i],xx[j]],[-1,-1]) for i in range(6) for j in range(i+1,6)]
    direct_five=combine([power(f,5) for f in linear],[1]*27)
    witness_terms=[(4,(5,0,0,0,0,0,0)),(-15,(4,0,0,0,0,0,1)),(15,(3,0,0,0,0,0,2)),
      (-15,(3,1,0,0,0,0,0)),(54,(2,1,0,0,0,0,1)),(-54,(1,1,0,0,0,0,2)),
      (27,(2,0,1,0,0,0,0)),(-81,(1,0,1,0,0,0,1)),(81,(0,0,1,0,0,0,2)),
      (-27,(1,0,0,1,0,0,0)),(81,(0,0,0,0,1,0,0))]
    expanded=[]
    for coefficient,exponents in witness_terms:
        f={(0,)*7:1}
        for g,e in zip(elementary+[uu],exponents):f=multiply(f,power(g,e))
        expanded.append(f)
    assert direct_five==combine(expanded,[5760*c for c,e in witness_terms])
    # Direct Laurent reconstruction of omega7 and the 27-term quadratic psi2.
    theta=[combine(moments,row) for row in inverse(gram)]
    xs=[combine(moments,[F(c,12) for c in row]) for row in X]
    z=combine(moments,[F(c,12) for c in Z])
    terms=[combine([x,z],[1,sign]) for x in xs for sign in (1,-1)]
    terms += [combine([xs[i],xs[j]],[-1,-1]) for i in range(6) for j in range(i+1,6)]
    quadratic=combine([multiply(f,f) for f in terms],[1]*27)
    inputs=[theta[6],quadratic]
    @lru_cache(None)
    def monomial(p):
        if not any(p):return {(0,)*7:F(1)}
        i=next(i for i,e in enumerate(p) if e);q=list(p);q[i]-=1
        return multiply(inputs[i],monomial(tuple(q)))
    last=next(r for r in data['direct_laurent_rows'] if r['degree']==3 and r['mode']=='filtered')
    support=sorted(set().union(*(set(monomial(tuple(p))) for p in last['monomial_exponents'])))
    positions=0
    for row in data['direct_laurent_rows']:
        polys=[monomial(tuple(p)) for p in row['monomial_exponents']]
        B=[[F(x) for x in r] for r in row['source_basis_matrix']];assert determinant(B)
        images=[combine(polys,[B[i][j] for i in range(len(B))]) for j in range(len(B))]
        assert all(c.denominator==1 for f in images for c in f.values())
        N=[[int(f.get(w,0)) for f in images] for w in support];divisor=0
        for minor in row['completeness_certificate']['maximal_minors']:
            det=abs(determinant([N[i] for i in minor['rows']]))
            assert det==int(minor['determinant']);divisor=gcd(divisor,int(det))
            assert modular_rank(N,minor['selection_prime'])==len(B)
        assert divisor==1;positions+=len(N)*len(B)
    oldpath=ROOT/'experiments/exceptional-minuscule-results.json';old=json.loads(oldpath.read_text())
    lattices=0;comparisons=0
    for row in data['degrees']:
        H=row['image_basis'];check_image_basis(row['image_matrix'],H)
        assert smith_square(H)==row['smith_factors'];lattices+=1
        for comparison in row['comparisons']:
            oldH=comparison['old_basis'];S=comparison['sum_basis'];I=comparison['intersection_basis']
            case=next(c for c in old['cases'] if c['kind']=='E7' and c['form']==comparison['form'])
            prior=case['degrees'][row['degree']]
            oldM=[prior['matrix'][prior['schubert_words'].index(w)] for w in row['schubert_words']]
            check_image_basis(oldM,oldH)
            assert smith_square(oldH)==comparison['old_smith_factors']
            check_image_basis([a+b for a,b in zip(H,oldH)],S)
            assert contained(H,oldH)==comparison['old_contains_character']
            assert contained(oldH,H)==comparison['character_contains_old']
            assert contained(I,H) and contained(I,oldH)
            assert abs(determinant(I)*determinant(S))==abs(determinant(H)*determinant(oldH))
            for key,A,B in [('sum_mod_character',H,S),('sum_mod_old',oldH,S),
                            ('old_mod_intersection',I,oldH),('character_mod_intersection',I,H)]:
                coordinates=matmul(inverse(B),A);integral_matrix(coordinates)
                assert smith_square(coordinates)==comparison[key]
            comparisons+=1;lattices+=1
    degree5=data['degrees'][5]
    original=next(c for c in old['cases'] if c['kind']=='E7' and c['form']=='image_denominators_cleared')['degrees'][5]
    j=original['monomials'].index('psi5')
    witness=[[F(216,5)*original['matrix'][original['schubert_words'].index(w)][j]] for w in degree5['schubert_words']]
    assert contained(witness,degree5['image_basis'])
    assert not contained(witness,degree5['comparisons'][0]['old_basis'])
    files=[path,oldpath,Path(__file__).resolve(),ROOT/'experiments/audit_e6_character_form.py']
    out=dict(result='passed' if complete else 'prefix_passed',degrees_checked=len(data['degrees']),
        engine='Python standard library only',weights=56,
        moment_reflection_checks=49,leading_characteristic_cases=2,integral_permutation_reflections=5,
        direct_laurent_lattices=len(data['direct_laurent_rows']),coefficient_positions=positions,
        image_lattices_and_smith_lists=lattices,embedded_sum_and_intersection_comparisons=comparisons,
        degree_five_witness_symbolic_and_lattice_checks=True,
        seconds=round(perf_counter()-start,3),source_hash_semantics='SHA256 with CRLF normalized to LF.',
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in files})
    if complete:
        (ROOT/'experiments/e7-character-form-audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
