"""Bind the proved 60-case E6 formula and three-case completion into 63 rows."""
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from summarize_e7_all_parabolics import divide, poincare, levi_degrees

ROOT = Path(__file__).resolve().parents[1]


def main():
    inputs = ['e6-root-weight-comparison.json', 'e6-remaining-a5-audit.json',
              'e6-remaining-a2a2a1-audit.json', 'e6-remaining-structure.json']
    data = {}
    for name in inputs:
        data[name] = json.loads((ROOT/'experiments'/name).read_text())
        assert data[name]['result'] == 'passed'
    uniform = {tuple(row['e6_retained']): row for row in
               data[inputs[0]]['proper_parabolics'] if row['common_fiber_basis']}
    assert len(uniform) == 60
    special = {}
    for nodes, name in (((1,3,4,5,6), inputs[1]), ((1,2,3,5,6), inputs[2])):
        rows = data[name]['degree_rows']
        special[nodes] = {str(3**a): [row['multiplicities'][a-1] for row in rows]
                          for a in (1,2,3)}
    special[(1,3,5,6)] = {
        factor: [coefficients[d] if d < len(coefficients) else 0 for d in range(31)]
        for factor, coefficients in special[(1,2,3,5,6)].items()}
    for factor, coefficients in special[(1,3,5,6)].items():
        base = special[(1,2,3,5,6)][factor]
        for d in range(1,31):
            coefficients[d] += base[d-1]
    W = poincare([2,5,6,8,9,12])
    rows, families = [], Counter()
    for size in range(6):
        for nodes in combinations(range(1,7), size):
            degrees, typ = levi_degrees(nodes)
            H = divide(W, poincare(degrees))
            dimension = 36-sum(d-1 for d in degrees)
            assert len(H) == dimension+1
            if nodes in uniform:
                family = 'good-Levi-60'
                counts = {str(3**a): uniform[nodes]['multiplicities'][str(a)] for a in (1,2,3)}
            else:
                family = {(1,3,4,5,6): 'A5', (1,2,3,5,6): 'A2+A2+A1',
                          (1,3,5,6): 'A2+A2-transfer'}[nodes]
                counts = special[nodes]
            counts = {factor: list(c)+[0]*(dimension+1-len(c)) for factor, c in counts.items()}
            assert all(len(c) == dimension+1 and min(c) >= 0 for c in counts.values())
            assert all(sum(c[d] for c in counts.values()) <= H[d] for d in range(dimension+1))
            families[family] += 1
            rows.append(dict(retained_nodes=list(nodes), levi_type=typ, dimension=dimension,
                             schubert_rank=sum(H), cohomology_ranks=H,
                             multiplicity_coefficients=counts, proof_family=family))
    paths = [Path(__file__).resolve(), ROOT/'experiments/summarize_e7_all_parabolics.py',
             ROOT/'experiments/e7_even_block_formula.py']+[ROOT/'experiments'/name for name in inputs]
    report = dict(result='passed', complete=True, cartan_type='E6',
                  proper_standard_parabolics=63, nonunit_factors=[3,9,27],
                  polynomial_degree_convention='d means H^(2d)',
                  root_source='Full (Sym_Z Q)^W_I', weight_source='Full (Sym_Z X)^W_I',
                  families=dict(families), parabolics=rows,
                  Springer_normalizations=dict(minimal_27='6^d M_I^d', adjoint_78='24^d N_I^d'),
                  scope='Complete integral relative image comparison. Existing source/image proofs are shared inputs, not 63 new direct computations. Absolute H/N and priority remain separate.',
                  source_sha256={p.relative_to(ROOT).as_posix(): sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'experiments/e6-all-parabolics-summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print('PASSED E6 all 63 proper parabolics', dict(families))


if __name__ == '__main__':
    main()
