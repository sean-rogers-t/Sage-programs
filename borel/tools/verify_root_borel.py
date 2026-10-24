"""Audit automatic root data, character integration and exact torus jets."""
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter
import unittest
from sage.all import QQ
from sage.env import SAGE_VERSION
from sage_borel import root_borel_map

ROOT = Path(__file__).resolve().parents[1]


def main():
    sys.path.insert(0, str(ROOT/'tests'))
    from test_root_borel import RootBorelTests, adjoint_data, product_unipotent, trace_pair
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(RootBorelTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)
    e8 = root_borel_map(('E', 8), (0, 0, 0, 0, 0, 0, 0, 1))
    u = tuple(QQ((-1)**i) for i in range(e8.dimension))
    start = perf_counter()
    actual = e8.coordinates((1,)*8, u)
    recurrence_seconds = perf_counter()-start
    T, positive, negative, H = adjoint_data(e8, (QQ.one(),)*8)
    g = product_unipotent(positive, u)
    expected = tuple(trace_pair(g, f)/trace_pair(e, f) for e, f in zip(positive, negative))
    assert actual == expected and len(actual) == 120 and e8.max_height == 29
    paths = ['src/sage_borel/root_borel.py', 'src/sage_borel/characters.py',
             'src/sage_borel/borel_reconstruction.py', 'src/sage_borel/geometry.py',
             'src/sage_borel/__init__.py', 'tests/test_root_borel.py',
             'tests/test_borel_coordinates.py', 'tests/test_borel_reconstruction.py',
             'ROOT_BOREL_API.md', 'tools/verify_root_borel.py']
    report = dict(result='passed', sage_version=SAGE_VERSION, test_methods=result.testsRun,
                  coefficient_field='QQ; theorem in characteristic zero',
                  torus_lattice='simply connected full weight lattice',
                  classical_standard_matrix_cases=12, previous_basis_convention_bridges=2,
                  larger_symmetric_square_trace_cases=2,
                  exceptional_adjoint_trace_cases=15, exceptional_minimal_conjugation_cases=4,
                  exceptional_types=['G2', 'F4', 'E6', 'E7', 'E8'],
                  symbolic_adapter_types=['A3', 'C3', 'G2'],
                  root_string_magnitudes_checked=True, mixed_character_normalization_checked=True,
                  root_hypertorus_evaluation='Finite torus jets, with precision highest-root height plus one',
                  E8_dense_identity=dict(positive_roots=120, adjoint_dimension=248,
                      nonzero_unipotent_parameters=120, jet_precision=30,
                      recurrence_seconds=round(recurrence_seconds, 6),
                      full_adjoint_trace_agrees=True),
                  independent_scope='Classical standard matrices are independently constructed. Exceptional '
                      'full adjoint trace and conjugation computations share Sage Chevalley constants but '
                      'do not use the Borel fixed-vector recurrence. No external review or separate CAS.',
                  limitations='Exact rational evaluation; symbolic expansion separately budgeted. '
                      'No complete singular-fiber, branch-locus, affine or all-isogeny-form classification.',
                  novelty_claim=False,
                  source_sha256={p: sha256((ROOT/p).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                                 for p in paths})
    (ROOT/'root-borel-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'source_sha256'}, indent=2))


if __name__ == '__main__':
    main()
