"""Source-bound exact symbolic verification of the Borel-coordinate API."""
from hashlib import sha256
import json
from pathlib import Path
import sys
import unittest
from sage.env import SAGE_VERSION

ROOT = Path(__file__).resolve().parents[1]


def main():
    sys.path.insert(0, str(ROOT/'tests'))
    from test_borel_coordinates import BorelCoordinateTests
    result = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(BorelCoordinateTests))
    if not result.wasSuccessful():
        raise SystemExit(1)
    paths = ['src/sage_borel/borel_coordinates.py', 'src/sage_borel/__init__.py',
             'tests/test_borel_coordinates.py', 'BOREL_COORDINATES.md',
             'tools/verify_borel_coordinates.py']
    report = dict(result='passed', sage_version=SAGE_VERSION, test_methods=result.testsRun,
                  groups=['SL3', 'Sp4'], representation='defining', characteristic=0,
                  independent_full_trace_gram_dimensions=[8, 10],
                  symbolic_all_parameter_matrix_identities=True,
                  symbolic_collapse_stratum_ranks=True,
                  exact_empty_and_positive_dimensional_fibers=True,
                  runtime_API_changed=True, novelty_claim=False,
                  scope='Fixed-t vertical maps; independent constructions in Sage, not external review. '
                        'No full B/G critical or target branch-locus classification.',
                  source_sha256={p: sha256((ROOT/p).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                                 for p in paths})
    (ROOT/'borel-coordinates-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'source_sha256'}, indent=2))


if __name__ == '__main__':
    main()
