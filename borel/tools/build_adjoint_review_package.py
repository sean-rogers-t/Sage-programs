"""Build a compact, hash-inventoried review archive from the pushed proof snapshot."""
import argparse
from hashlib import sha256
import io
import json
from pathlib import Path
import subprocess
import tarfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
REVIEW = 'review/adjoint-generic-fibers'
BASE = '4c21f18ea0c5df4350d86131176739b60f35ac6e'
ADDED = [REVIEW+'/README.md', REVIEW+'/CERTIFICATE_GUIDE.md',
         REVIEW+'/run_review.py', 'tools/build_adjoint_review_package.py']


def payload(name):
    data = (ROOT/name).read_bytes()
    return data if name.endswith('.gz') else data.replace(b'\r\n', b'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, help='New ZIP filename')
    args = parser.parse_args()
    output = Path(args.output).resolve()
    if output.exists():
        raise FileExistsError(output)
    pending = ['e8-last-case-audit.json', 'adjoint-theorem-audit.json']
    proof = set()
    while pending:
        name = pending.pop()
        if name in proof:
            continue
        path = (ROOT/name).resolve()
        if not path.is_relative_to(ROOT.resolve()):
            raise ValueError('Out-of-package input: '+name)
        proof.add(name)
        if path.suffix == '.json':
            data = json.loads(path.read_text())
            if isinstance(data, dict):
                for key in ('source_sha256', 'raw_sha256'):
                    pending.extend(data.get(key, {}))
    # Bind the supplied mathematics to the already pushed commit. Newly
    # prepared review instructions are listed separately in the manifest.
    archive = subprocess.check_output(['git', 'archive', '--format=tar', BASE,
        *['borel/'+name for name in sorted(proof)]], cwd=ROOT.parent)
    with tarfile.open(fileobj=io.BytesIO(archive)) as source:
        for name in sorted(proof):
            data = source.extractfile('borel/'+name).read()
            if not name.endswith('.gz'):
                data = data.replace(b'\r\n', b'\n')
            if data != payload(name):
                raise ValueError('Proof input differs from pinned commit: '+name)
    names = sorted(proof | set(ADDED))
    manifest = dict(format_version=1, prepared='2026-10-27', proof_commit=BASE,
        scope='Consolidated proof dependencies and final E8 replay inputs; not the full historical repository.',
        new_review_files=ADDED,
        files=[dict(path=n, bytes=len(payload(n)), sha256=sha256(payload(n)).hexdigest(),
                    hash_mode='raw' if n.endswith('.gz') else 'text-LF') for n in names])
    manifest_name = REVIEW+'/package-manifest.json'
    (ROOT/manifest_name).write_text(json.dumps(manifest, indent=2)+'\n')
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, 'x', compression=zipfile.ZIP_DEFLATED) as bundle:
        for name in names+[manifest_name]:
            info = zipfile.ZipInfo('borel/'+name, date_time=(2026, 10, 27, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            bundle.writestr(info, payload(name))
    # Verify the archived bytes, not only the source files.
    with zipfile.ZipFile(output) as bundle:
        assert bundle.testzip() is None
        for entry in manifest['files']:
            data = bundle.read('borel/'+entry['path'])
            assert sha256(data).hexdigest() == entry['sha256']
    report = dict(result='passed', proof_commit=BASE, archive_filename=output.name,
                  archive_sha256=sha256(output.read_bytes()).hexdigest(), archive_bytes=output.stat().st_size,
                  proof_files_matched_to_commit=len(proof), input_files=len(names),
                  archived_bytes_checked=True, replay_from_extracted_archive=False,
                  review_not_external=True)
    (ROOT/REVIEW/'package-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
