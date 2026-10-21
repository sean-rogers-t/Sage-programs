"""Run geometry regressions from an extracted wheel, outside the source tree."""
from pathlib import Path
from hashlib import sha256
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    version = re.search(r'^version = "([^"]+)"', (ROOT/'pyproject.toml').read_text(), re.MULTILINE).group(1)
    wheel = ROOT/('dist/sage_borel-%s-py3-none-any.whl' % version)
    with tempfile.TemporaryDirectory(prefix='borel-geometry-wheel-') as temporary:
        target = Path(temporary)
        with zipfile.ZipFile(wheel) as archive:
            archive.extractall(target)
        tests = target/'wheel_tests'
        tests.mkdir()
        shutil.copyfile(ROOT/'tests/test_geometry.py', tests/'test_geometry.py')
        shutil.copyfile(ROOT/'tests/test_e6_geometry.py', tests/'test_e6_geometry.py')
        shutil.copyfile(ROOT/'tests/test_e7_geometry.py', tests/'test_e7_geometry.py')
        shutil.copyfile(ROOT/'tests/test_exceptional_degrees.py', tests/'test_exceptional_degrees.py')
        environment = dict(os.environ, PYTHONPATH=str(target))
        result = subprocess.run([sys.executable, '-m', 'unittest', 'discover',
                                 '-s', str(tests), '-v'], cwd=target,
                                env=environment, text=True, capture_output=True)
        print(result.stdout+result.stderr)
        if result.returncode:
            raise SystemExit(result.returncode)
        imported = subprocess.run([sys.executable, '-c',
            'import sage_borel; print(sage_borel.__file__); print(sage_borel.__version__)'],
            cwd=target, env=environment, text=True, capture_output=True, check=True)
        assert str(target) in imported.stdout and version in imported.stdout
    count = int(re.search(r'Ran (\d+) tests', result.stderr).group(1))
    report = dict(result='passed', version=version, tests=count,
                  command='sage -python tools/verify_geometry_wheel.py',
                  isolated_import_verified=True, wheel=wheel.name,
                  wheel_sha256=sha256(wheel.read_bytes()).hexdigest())
    (ROOT/'geometry-wheel-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
