"""Remove redundant expanded higher polynomials from the source certificate.

F and C determine A=F C uniquely. The native verifier reconstructs every
coefficient of A, checks integrality and its full primitive mod-two rank.
The seed polynomials are retained for direct reading and the first witness.
"""
from pathlib import Path
from hashlib import sha256
import gzip,json
ROOT=Path(__file__).resolve().parents[1]

def main():
 path=ROOT/'experiments/e7-d6-sources.json';report=json.loads(path.read_text());archive=ROOT/report['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==report['certificate_archive_sha256']
 with gzip.open(archive,'rt') as stream:data=json.load(stream)
 for source in data.values():
  for row in source['degrees']:
   if row['degree']>20:row.pop('basis_polynomials',None)
 with gzip.GzipFile(filename=str(archive),mode='wb',mtime=0) as stream:stream.write(json.dumps(data,separators=(',',':')).encode())
 report['certificate_archive_sha256']=sha256(archive.read_bytes()).hexdigest();report['certificate_encoding']='Expanded seed through degree 20; higher basis polynomials are uniquely specified by rational coordinates C in the classical basis F. Reconstruct A=F C and check integrality and primitive rank.';report['source_sha256'][Path(__file__).relative_to(ROOT).as_posix()]=sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest();path.write_text(json.dumps(report,indent=2)+'\n');print('packed bytes',archive.stat().st_size)

if __name__=='__main__':main()
