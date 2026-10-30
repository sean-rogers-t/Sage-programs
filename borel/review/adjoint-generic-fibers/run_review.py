"""Check package integrity, then replay selected checks in an isolated copy."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT/'review/adjoint-generic-fibers/package-manifest.json'


def checked_inputs():
    manifest = json.loads(MANIFEST.read_text())
    seen = set()
    for entry in manifest['files']:
        name = entry['path']
        path = (ROOT/name).resolve()
        if not path.is_relative_to(ROOT.resolve()) or name in seen:
            raise ValueError('Unsafe or duplicate manifest path: '+name)
        seen.add(name)
        data = path.read_bytes()
        if entry['hash_mode'] == 'text-LF':
            data = data.replace(b'\r\n', b'\n')
        if len(data) != entry['bytes'] or sha256(data).hexdigest() != entry['sha256']:
            raise ValueError('Manifest mismatch: '+name)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('integrity', 'native', 'full'), default='integrity')
    parser.add_argument('--engine', help='Built pinned e8-last-f4 executable (full mode)')
    parser.add_argument('--sage', default='sage', help='Sage executable (full mode)')
    parser.add_argument('--output', help='New output directory; must not already exist')
    args = parser.parse_args()
    manifest = checked_inputs()
    print('Manifest passed:', len(manifest['files']), 'files', flush=True)
    if args.mode == 'integrity':
        return
    if args.mode == 'full' and not args.engine:
        parser.error('--engine is required for full mode')
    engine = str(Path(args.engine).resolve()) if args.engine else None
    if args.output:
        out = Path(args.output).resolve()
        out.mkdir(parents=True, exist_ok=False)
    else:
        out = Path(tempfile.mkdtemp(prefix='adjoint-review-'))
    work = out/'borel'
    for entry in manifest['files']:
        target = work/entry['path']
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/entry['path'], target)
    print('Isolated reports and logs:', out, flush=True)
    steps = []
    if args.mode == 'full':
        for prime in (32003, 32009):
            steps.append((f'e8-membership-{prime}', [args.sage, '-python',
                'tools/verify_e8_last_case.py', '--prime', str(prime), '--engine', engine]))
    steps += [
        ('e8-native', [sys.executable, 'tools/audit_e8_last_case.py']),
        ('coverage', [sys.executable, 'tools/audit_adjoint_theorem.py']),
    ]
    completed = []
    for name, command in steps:
        print('Running:', name, flush=True)
        log = out/(name+'.log')
        with log.open('w') as stream:
            result = subprocess.run(command, cwd=work, stdout=stream, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f'{name} failed ({result.returncode}); inspect {log}')
        completed.append(name)
        print('Passed:', name, flush=True)
    reports = ['e8-last-case-audit.json', 'adjoint-theorem-audit.json']
    if args.mode == 'full':
        reports += [f'e8-last-{p}-verification.json' for p in (32003, 32009)]
    for name in reports:
        if json.loads((work/name).read_text())['result'] != 'passed':
            raise RuntimeError('Report did not pass: '+name)
    summary = dict(result='passed', mode=args.mode, proof_commit=manifest['proof_commit'],
                   input_files=len(manifest['files']), steps=completed,
                   geometric_proof_mechanically_verified=False,
                   full_historical_case_replay=False,
                   original_package_modified=False)
    (out/'review-run.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
