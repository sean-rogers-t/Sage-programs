"""Compress completed proof certificates and retain readable JSON summaries."""
from hashlib import sha256
from pathlib import Path
import gzip,json

ROOT=Path(__file__).resolve().parents[1]


def main():
    path=ROOT/'experiments/e6-remaining-a2a2a1.json'
    data=json.loads(path.read_text())
    if 'certificate_archive' in data:
        print('Already packaged',data['certificate_archive']);return
    assert data['result']=='passed' and len(data['degrees'])==30
    packed=json.dumps(data,separators=(',',':')).encode()
    archive=path.with_name('e6-remaining-a2a2a1-certificates.json.gz')
    archive.write_bytes(gzip.compress(packed,compresslevel=9,mtime=0))
    summary=dict(data)
    summary['sources']=[dict(degree=r['degree'],rank=len(r['exponents']),coefficient_rows=r['coefficient_rows'],
                             certificate='Full basis and primitive-column certificate in the archive.') for r in data['sources']]
    summary['degrees']=[{k:v for k,v in r.items() if k not in ('source_labels','root_image_columns')} for r in data['degrees']]
    summary['certificate_archive']=archive.relative_to(ROOT).as_posix()
    summary['certificate_archive_sha256']=sha256(archive.read_bytes()).hexdigest()
    summary['archive_hash_semantics']='Raw bytes; gzip mtime fixed at zero.'
    path.write_text(json.dumps(summary,indent=2)+'\n')
    print('Packaged',len(packed),'to',archive.stat().st_size,'bytes',flush=True)


if __name__=='__main__':main()
