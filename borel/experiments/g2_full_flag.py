"""Full literal G2/B character images, exact integral kernel and Smith factors."""
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter
from sage.all import ZZ, matrix
from sage_borel import g2_character_algebra, springer_character_map

ROOT = Path(__file__).resolve().parents[1]


def main():
    start = perf_counter()
    reports = []
    for representation, highest in (('minimal', (1,0)), ('adjoint', (0,1))):
        algebra = g2_character_algebra(representation=representation)
        character_map = springer_character_map('G2', highest)
        scale = 4 if representation == 'adjoint' else 1
        assert character_map.linear_basis_matrix == scale*matrix(ZZ, [[6,0],[0,2]])
        X, Y = algebra.domain.gens()
        relations = [X**2-3*X*Y+3*Y**2, Y**6]
        assert all(not algebra.borel_map.expand(algebra.substitute(f)) for f in relations)
        # Both polynomials are monic, with coprime leading monomials X^2,Y^6.
        # Their integer quotient is free on Y^j, X Y^j (0<=j<6).
        basis = [(0,j) for j in range(6)]+[(1,j) for j in range(6)]
        rows = []
        for d in range(8):
            row = algebra.integral_degree(d, include_ambient=True)
            selected = [p for p in basis if sum(p) == d]
            assert len(selected) == row.rank == row.ambient_dimension
            normal = matrix(ZZ, [[algebra.evaluate(algebra.domain.monomial(*p)).get(w, 0)
                                  for p in selected] for w in row.schubert_indices])
            # The free quotient injects: its homogeneous basis has full rank.
            assert normal.rank() == len(selected)
            smith = tuple(abs(int(c)) for c in normal.elementary_divisors() if c)
            assert smith == tuple(map(int, row.smith_factors))
            rows.append(dict(degree=d, cohomological_degree=2*d,
                             rank=row.rank, ambient_rank=row.ambient_dimension,
                             source_exponents=[list(map(int,p.exponents()[0])) for p in row.monomials],
                             schubert_words=[list(map(int,w.reduced_word())) for w in row.schubert_indices],
                             coefficient_matrix=[list(map(int, r)) for r in row.matrix],
                             quotient_basis_exponents=[list(p) for p in selected],
                             smith_factors=list(map(int, row.smith_factors))))
        reports.append(dict(representation=representation, dimension=7 if scale == 1 else 14,
                            polynomial_source='Z[x,y]', linear_generators=[6*scale, 2*scale],
                            integral_kernel=['x^2-3*x*y+3*y^2', 'y^6'],
                            free_quotient_rank=12, degrees=rows))
    for a, b in zip(reports[0]['degrees'], reports[1]['degrees']):
        assert b['smith_factors'] == [4**a['degree']*c for c in a['smith_factors']]
    paths = [Path(__file__).resolve(), ROOT/'src/sage_borel/characters.py',
             ROOT/'G2_CHARACTER_FORM.md', ROOT/'experiments/g2-character-form-audit.json']
    result = dict(result='passed', complete=True, group='G2', parabolic='B',
                  source='Full literal integral polynomial-character source for normalized trace projection',
                  polynomial_degree_convention='d means H^(2d)',
                  cohomology_ranks=[1,2,2,2,2,2,1], cases=reports,
                  integral_kernel_proof='Monic relations with coprime leading monomials; free rank-12 quotient; all 12 basis images rationally independent. Hence no additional integer kernel. Source completeness is the earlier full torus-source theorem.',
                  scope='Complete absolute embedded images in all degrees for minimal and adjoint G2/B; no literature-priority claim.',
                  seconds=round(perf_counter()-start, 3),
                  source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'experiments/g2-full-flag.json').write_text(json.dumps(result, indent=2)+'\n')
    print('PASSED complete G2/B images', [(r['degree'],r['smith_factors']) for r in reports[0]['degrees']])


if __name__ == '__main__':
    main()
