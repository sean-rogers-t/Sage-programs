"""Record the full regression run for the character-source API release."""
import contextlib
from hashlib import sha256
import io
import json
from pathlib import Path
from time import perf_counter
import sys
import unittest
from sage.env import SAGE_VERSION
import sage_borel

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    started = perf_counter()
    output = io.StringIO()
    suite = unittest.defaultTestLoader.discover(str(ROOT/'tests'))
    with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
        result = unittest.TextTestRunner(stream=output, verbosity=1).run(suite)
    files = sorted((ROOT/'src/sage_borel').glob('*.py'))
    files += sorted((ROOT/'tests').glob('*.py'))
    files += sorted((ROOT/'src/sage_borel/data').glob('*.json'))
    files += [ROOT/'tests/data/exceptional_expected.json',
              ROOT/'tests/data/e7_relative_audited_bases.json', ROOT/'pyproject.toml',
              ROOT/'tests/data/e6_relative_audited_base.json',
              Path(__file__).resolve()]
    report = dict(result='passed' if result.wasSuccessful() else 'failed',
                  version=sage_borel.__version__, sage_version=SAGE_VERSION,
                  test_methods_run=result.testsRun, failures=len(result.failures),
                  errors=len(result.errors), skipped=len(result.skipped),
                  command='PYTHONPATH=src sage -python tools/verify_character_release.py',
                  seconds=round(perf_counter()-started, 3),
                  source_hash_semantics='SHA256 with CRLF normalized to LF.',
                  source_sha256={p.relative_to(ROOT).as_posix():
                     sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest() for p in files})
    (ROOT/'character-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    if not result.wasSuccessful():
        print(output.getvalue())
        raise SystemExit(1)
    print(json.dumps({k: v for k, v in report.items() if k != 'source_sha256'}, indent=2))


if __name__ == '__main__':
    main()
