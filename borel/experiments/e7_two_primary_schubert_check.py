"""Independent Chevalley basis and direct BGG bridge for the E7/P7 proof."""
from collections import defaultdict
from hashlib import sha256
from pathlib import Path
from time import perf_counter
import json,sys
from sage.all import QQ,ZZ,matrix,vector

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from tests.test_exceptional import RootChevalleyOracle
from sage_borel import exceptional_borel


def local(A):return A.denominator()%2!=0


def main():
    start=perf_counter();b=exceptional_borel('E7',7);oracle=RootChevalleyOracle('E7')
    paths=[ROOT/'experiments/e7-character-form-results.json',ROOT/'experiments/e7-adjoint-comparison-results.json']
    oldroot,oldweight=[json.loads(p.read_text()) for p in paths]
    anchors={0:(),5:(2,4,5,6,7),9:(1,5,4,2,3,4,5,6,7)}
    towers={}
    for birth,word in anchors.items():
        expansion={oracle.W.from_reduced_word(word):QQ.one()}
        for d in range(birth,28):
            towers[birth,d]={tuple(w.reduced_word()):c for w,c in expansion.items()}
            # Weyl words need not be canonical under independent parents.
            towers[birth,d]={b.W.from_reduced_word(w):c for w,c in towers[birth,d].items()}
            out=defaultdict(QQ)
            for w,c in expansion.items():
                for v,weights in oracle.covers(w):
                    if weights[6]:out[v]+=c*weights[6]
            expansion={w:c for w,c in out.items() if c}
    # Literal published generator images in these three Chevalley towers.
    formulas=[('omega7',1,[(0,1)]),('psi2',2,[(0,-144)]),
              ('psi5',5,[(0,403200),(5,-933120)]),
              ('psi6',6,[(0,4720896),(5,-11197440)]),
              ('psi8',8,[(0,-622420992),(5,1463132160)]),
              ('psi9',9,[(0,-1289723904),(5,-188116992),(9,15237476352)]),
              ('psi12',12,[(0,3066946633728),(5,-23932781199360),(9,78538306682880)])]
    for name,d,terms in formulas:
        # Expand the actual polynomial directly by divided differences.
        # Runtime catalogue data contains input definitions, not image coefficients.
        direct=b.expand(b.generator_polynomial(name))
        expected=defaultdict(QQ)
        for birth,c in terms:
            for w,a in towers[birth,d].items():expected[w]+=c*a
        assert {w:c for w,c in direct.items() if c}=={w:c for w,c in expected.items() if c},name
        print('Fresh BGG generator',name,round(perf_counter()-start,2),flush=True)
    checks=[]
    for d in range(28):
        labels=list(b.basis(d));basis=[(d,0)]
        for birth in (5,9):
            if birth<=d<birth+14:basis.append((d-birth,birth))
        columns=[vector(QQ,[(1 if birth==0 else 2)*towers[birth,d].get(w,0) for w in labels]) for i,birth in basis]
        N=matrix(QQ,columns).transpose();assert N.is_square() and N.det()!=0
        M=matrix(QQ,[c*(2 if i%2 else 1) for c,(i,birth) in zip(columns,basis)]).transpose()
        for data,key,factor,proposed in ((oldroot,'image_basis',12,M),(oldweight,'image_basis',36,N)):
            row=data['degrees'][d];words=[b.W.from_reduced_word(w) for w in row['schubert_words']]
            old=matrix(QQ,[row[key][words.index(w)] for w in labels])/(factor**d)
            change=proposed.inverse()*old
            assert local(change) and change.det().valuation(2)==0
        checks.append(dict(degree=d,rank=len(labels),basis_change_2_unit=True))
    files=[Path(__file__).resolve(),ROOT/'tests/test_exceptional.py']+paths
    report=dict(result='passed',fresh_bgg_generators=7,chevalley_basis_degrees=28,
                total_basis_classes=56,degrees=checks,seconds=round(perf_counter()-start,3),
                scope='Direct divided differences for all seven catalogue generators; positive-root Chevalley multiplication from published h,s,t Schubert anchors for every proposed image basis. Historical image lattices are compared, not taken as a new source proof.',
                source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in files})
    (ROOT/'experiments/e7-two-primary-schubert-check.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASSED',report['seconds'],flush=True)


if __name__=='__main__':main()
