"""Generate readable all-degree tables only from fully passed proof audits."""
from hashlib import sha256
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]


def read(name):
    data=json.loads((ROOT/'experiments'/name).read_text())
    assert data['result']=='passed',name
    return data


def main():
    aa=read('e6-remaining-a5-audit.json');bb=read('e6-remaining-a2a2a1-audit.json')
    A=read('e6-remaining-a5.json');B=read('e6-remaining-a2a2a1.json')
    bgg=read('e6-remaining-bgg-check.json');structure=read('e6-remaining-structure.json')
    for audit,report,top in ((aa,A,21),(bb,B,29)):
        assert len(audit['degree_rows'])==len(audit['source_rows'])==len(report['degrees'])==top+1
        for path,expected in audit['source_sha256'].items():
            assert sha256((ROOT/path).read_bytes().replace(b'\r\n',b'\n')).hexdigest()==expected,path
        for path,expected in audit.get('raw_sha256',{}).items():
            assert sha256((ROOT/path).read_bytes()).hexdigest()==expected,path
        for primary,independent in zip(report['degrees'],audit['degree_rows']):
            assert primary['rank']==independent['rank']
            assert [primary['smith_exponents'].count(a) for a in (1,2,3)]==independent['multiplicities']
    for path,expected in bgg['source_sha256'].items():
        assert sha256((ROOT/path).read_bytes().replace(b'\r\n',b'\n')).hexdigest()==expected,path
    for path,expected in bgg.get('raw_sha256',{}).items():
        assert sha256((ROOT/path).read_bytes()).hexdigest()==expected,path
    lines=['# Complete tables for the three exceptional E6 parabolics','',
      'Generated from the passed all-degree certificates and independent audits.',
      'See the [theorem and proof](E6_REMAINING_CASES.md) and',
      '[proof/source audit](E6_REMAINING_AUDIT.md). Polynomial degree d corresponds',
      'to H^(2d). Each tuple (m1,m2,m3) denotes',
      '(Z/3)^m1 plus (Z/9)^m2 plus (Z/27)^m3 in N_I^d/M_I^d.',
      'The rank is the degree-d Schubert rank; its difference from m1+m2+m3',
      'is the number of unit inclusion factors. Values outside the tables are zero.','']
    tables=[('I_A={1,3,4,5,6}, derived Levi A5',aa['degree_rows']),
            ('I_B={1,2,3,5,6}, derived Levi A2+A2+A1',bb['degree_rows'])]
    base={r['degree']:r for r in bb['degree_rows']}
    extension=[]
    for d in range(31):
        now=base.get(d,dict(rank=0,multiplicities=[0]*3))
        prev=base.get(d-1,dict(rank=0,multiplicities=[0]*3))
        extension.append(dict(degree=d,rank=now['rank']+prev['rank'],
                        multiplicities=[a+b for a,b in zip(now['multiplicities'],prev['multiplicities'])]))
    assert sum(r['rank'] for r in extension)==1440
    tables.append(('I_C={1,3,5,6}, derived Levi A2+A2',extension))
    for kind,(_,entries) in zip(('A5','A2A2A1','A2A2'),tables):
        assert [r['rank'] for r in entries]==structure['schubert_ranks'][kind]
    for title,entries in tables:
        lines+=['## '+title,'','| d | Schubert rank | m1 | m2 | m3 |','|---:|---:|---:|---:|---:|']
        for row in entries:lines.append('| '+' | '.join(map(str,[row['degree'],row['rank']]+row['multiplicities']))+' |')
        totals=[sum(r['multiplicities'][i] for r in entries) for i in range(3)]
        lines+=['',f'Total multiplicities: {tuple(totals)}. Total Schubert rank: {sum(r["rank"] for r in entries)}.','']
    lines+=['The third table is deduced from the proved common (1,z) basis:',
            'each row is the sum of the two adjacent rows of the second table.',
            'It is not a separate all-degree BGG computation.','']
    (ROOT/'E6_REMAINING_TABLES.md').write_text('\n'.join(lines))
    print('Generated three complete tables from 52 audited base degrees and the structural P1 extension.')


if __name__=='__main__':main()
