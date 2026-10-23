"""Bind the torus/root reconstruction proof to independent matrix tests."""
from hashlib import sha256
import json
from pathlib import Path
import sys
import unittest
from sage.env import SAGE_VERSION

ROOT = Path(__file__).resolve().parents[1]


def main():
    sys.path.insert(0, str(ROOT/'tests'))
    from test_borel_reconstruction import BorelReconstructionTests
    result = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(BorelReconstructionTests))
    if not result.wasSuccessful():
        raise SystemExit(1)
    from sage_borel import reconstruct_borel_map
    omega2 = reconstruct_borel_map('Sp4', 'second_fundamental')
    x, y = omega2.torus_parameters
    assert omega2.vertical_jacobian().det() == (x+y)*(x*y+1)*(x*x+1)*(y*y+1)/(16*y)
    paths = ['src/sage_borel/borel_reconstruction.py', 'src/sage_borel/borel_coordinates.py',
             'src/sage_borel/__init__.py', 'tests/test_borel_reconstruction.py',
             'tests/test_borel_coordinates.py', 'BOREL_RECONSTRUCTION.md',
             'tools/verify_borel_reconstruction.py']
    report = dict(result='passed', sage_version=SAGE_VERSION, test_methods=result.testsRun,
                  ready_backends=['SL3', 'Sp4'], characteristic=0,
                  symbolic_full_representation_projection_cases=9,
                  low_level_additional_defining_matrix_checks=['SL2', 'SL4'],
                  independently_built_representations=['defining', 'dual', 'symmetric square',
                                                       'exterior square', 'adjoint', 'direct sum'],
                  independent_arithmetic='Different constructions inside Sage; not a separate CAS or external review.',
                  representation_matrices_used_by_runtime=False,
                  regular_extension_checked=True, changed_collapse_loci_checked=True,
                  printed_Sp4_omega2_determinant_verified=True,
                  general_recurrence_proof='BOREL_RECONSTRUCTION.md',
                  scope='Uniform characteristic-zero recurrence for consistent finite root data and a genuine '
                        'equivariant regular map. Ready coordinate adapters SL3/Sp4. '
                        'No arbitrary-character positivity certification or general singular-fiber classification.',
                  novelty_claim=False,
                  source_sha256={p: sha256((ROOT/p).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                                 for p in paths})
    (ROOT/'borel-reconstruction-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'source_sha256'}, indent=2))


if __name__ == '__main__':
    main()
