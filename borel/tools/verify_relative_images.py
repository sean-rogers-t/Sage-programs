"""Run and bind the independent E7 relative-image API checks."""
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from time import perf_counter
import unittest
import sage_borel
from package_e7_relative_images import packaged_data

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def main():
    start = perf_counter()
    package = ROOT/'src/sage_borel/data/e7_relative_images.json'
    assert json.loads(package.read_text()) == packaged_data()
    fixture = ROOT/'tests/data/e7_relative_audited_bases.json'
    for row in json.loads(fixture.read_text()).values():
        assert digest(ROOT/row['source']) == row['source_sha256']
        source = json.loads((ROOT/row['source']).read_text())
        counts = source.get('multiplicity_coefficients')
        if counts is None:
            counts = [0]+[r['multiplicity_2'] for r in source['degree_rows']]
        assert counts == row['coefficients']
    spec = importlib.util.spec_from_file_location('relative_checks', ROOT/'tests/test_relative_images.py')
    checks = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checks)
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(checks))
    if not result.wasSuccessful():
        raise SystemExit(1)
    paths = [package, fixture, ROOT/'tests/test_relative_images.py', Path(__file__).resolve(),
             ROOT/'tools/package_e7_relative_images.py', ROOT/'src/sage_borel/relative_images.py',
             ROOT/'experiments/e7-all-parabolics-summary.json']
    report = dict(result='passed', version=sage_borel.__version__,
                  tests=result.testsRun, proper_standard_parabolics=127,
                  cohomology_polynomials_reconstructed=127,
                  multiplicity_polynomials_reconstructed=125,
                  multiplicity_polynomials_from_independent_image_audits=2,
                  transformed_integer_lattice_cases=28,
                  fresh_e7_p7_character_degrees=[0,1,2,3,6,10],
                  containment_and_four_quotients_checked=True,
                  scope='Independent formula arithmetic, integer-lattice operations and bounded fresh Schubert images. Shared full-source theorems, A5 intermediary polynomial and previously audited D5+A1/D4+A1 base coefficients; no new source proof or novelty claim.',
                  seconds=round(perf_counter()-start, 3),
                  source_hash_semantics='SHA256 with CRLF normalized to LF',
                  source_sha256={p.relative_to(ROOT).as_posix(): digest(p) for p in paths})
    (ROOT/'relative-image-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'}, indent=2))


if __name__ == '__main__':
    main()
