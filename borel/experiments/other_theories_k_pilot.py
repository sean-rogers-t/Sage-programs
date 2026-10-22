"""Research pilot: literal Springer characters under the K characteristic map.

K coordinates are Euler characteristics on Schubert varieties, obtained by
isobaric Demazure operators. The Euler pairing is integral and unimodular.
Ordinary coordinates use the separate BGG engine. This is an experiment,
not a new public K-theory backend.
"""
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from time import perf_counter

from sage.all import QQ, ZZ, RootSystem, matrix, prod
from sage_borel import BorelMap, SpringerCharacterMap

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def multiply(f,g):
    out={}
    for a,c in f.items():
        for b,d in g.items():
            p=tuple(x+y for x,y in zip(a,b));out[p]=out.get(p,0)+c*d
    return {p:c for p,c in out.items() if c}


def demazure(f,i,alpha):
    out={}
    for weight,c in f.items():
        k=weight[i]
        terms=((j,1) for j in range(k+1)) if k>=0 else ((-j,-1) for j in range(1,-k))
        for j,sign in terms:
            p=tuple(x-j*y for x,y in zip(weight,alpha));out[p]=out.get(p,0)+sign*c
    return {p:c for p,c in out.items() if c}


def euler(f,word,roots):
    for node in reversed(word): f=demazure(f,node-1,roots[node-1])
    return sum(f.values())


def run(kind,highest):
    start=perf_counter();character=SpringerCharacterMap(kind,highest)
    borel=BorelMap(kind);r=character.rank
    words=sorted((tuple(map(int,w.reduced_word())) for w in borel.W),key=lambda w:(len(w),w))
    elements=[borel.W.from_reduced_word(w) for w in words];top=max(map(len,words))
    lattice=RootSystem(kind).weight_lattice();rs=lattice.simple_roots()
    roots=[tuple(int(rs[i+1].monomial_coefficients().get(j+1,0)) for j in range(r)) for i in range(r)]
    pulls=[{tuple(map(int,p)):int(c) for p,c in f.dict().items()} for f in character.linear_pullbacks]
    exponents=sorted((p for p in product(range(top+1),repeat=r) if sum(p)<=top),key=lambda p:(sum(p),p))
    Kcols=[];Hcols=[]
    for powers in exponents:
        f={(0,)*r:1};h=borel.ring.one()
        for j,power in enumerate(powers):
            for _ in range(power):f=multiply(f,pulls[j])
            ell=sum(character.linear_basis_matrix[i,j]*borel.weights[i] for i in range(r))
            h*=ell**power
        Kcols.append([euler(f,w,roots) for w in words])
        image=borel.expand(h)
        Hcols.append([int(image.get(w,0)) for w in elements])
    K=matrix(ZZ,Kcols).transpose();H=matrix(ZZ,Hcols).transpose()
    assert K.change_ring(QQ).row_space()==H.change_ring(QQ).row_space()
    ks=[abs(int(x)) for x in K.elementary_divisors() if x]
    hs=[abs(int(x)) for x in H.elementary_divisors() if x]
    assert len(ks)==len(words)==len(hs) and prod(ks)==prod(hs)
    report=dict(cartan_type=kind,highest_weights=[list(w) for w in highest],representation_dimension=character.dimension,
                linear_lattice_rows=[list(map(int,row)) for row in character.linear_basis_matrix],
                positive_roots=top,ambient_rank=len(words),source_monomials=len(exponents),
                same_rational_source_kernel=True,K_smith_factors=ks,ordinary_total_smith_factors=hs,
                same_total_index=True,total_index=int(prod(ks)),seconds=round(perf_counter()-start,3))
    cert=dict(cartan_type=kind,roots=[list(a) for a in roots],words=[list(w) for w in words],
              exponents=[list(p) for p in exponents],
              pullbacks=[[dict(weight=list(p),coefficient=c) for p,c in sorted(f.items())] for f in pulls],
              K_rows=[list(map(int,row)) for row in K],H_rows=[list(map(int,row)) for row in H])
    print(kind,highest,'K',ks,'H',hs,flush=True)
    return report,cert


def same_lattice_pair():
    results=[];b=BorelMap('A2',omitted_nodes=1)
    for last in ((2,0),(0,2)):
        ch=SpringerCharacterMap('A2',((1,0),)*4+(last,))
        x,y=ch.ring.gens();f1=27*x;f2=81*(2*x-y)*(x+y)
        roots=[(2,-1),(-1,2)];words=[(),(1,),(2,1)]
        fs=[ch.torus_ring.one(),ch.pullback(f1),ch.pullback(f2)]
        columns=[]
        for f in fs:
            d={tuple(map(int,p)):int(c) for p,c in f.dict().items()}
            columns.append([euler(d,w,roots) for w in words])
        K=matrix(ZZ,columns).transpose();ks=[abs(int(c)) for c in K.elementary_divisors()]
        assert list(b.expand(b.ring(f1)).values())==[27]
        assert list(b.expand(b.ring(f2)).values())==[243]
        # Equality with 9Q is checked as embedded column lattices.
        expected=matrix(ZZ,[[18,-9],[-9,18]])
        assert ch.linear_basis_matrix.change_ring(ZZ).transpose().row_module()==expected.transpose().row_module()
        results.append(dict(last_highest_weight=list(last),dimension=ch.dimension,
                            lattice='9Q',euler_matrix=[list(map(int,row)) for row in K],
                            K_smith_factors=ks,ordinary_smith_factors=[1,27,243],
                            pullbacks=[[dict(weight=list(p),coefficient=int(c)) for p,c in f.dict().items()] for f in fs]))
    assert results[0]['K_smith_factors']==[1,3,2187]
    assert results[1]['K_smith_factors']==[1,9,729]
    return results


def main():
    cases=[];certificates=[]
    for kind,high in [('A2',((1,0),)),('A2',((1,1),)),('C2',((1,0),)),
                      ('G2',((1,0),)),('G2',((0,1),)),('A3',((1,0,0),))]:
        a,c=run(kind,high);cases.append(a);certificates.append(c)
    pair=same_lattice_pair()
    archive=ROOT/'experiments/other-theories-k-certificates.json'
    archive.write_text(json.dumps(dict(cases=certificates),indent=2)+'\n')
    report=dict(result='passed',date='2026-10-22',cases=cases,same_lattice_different_K_pair=pair,
                scope='Finite full-flag checks and a complete projective-plane pair; corroborates a separately written general kernel/filtered-image argument. No novelty claim or gamma-filtration identification.',
                source_sha256={'experiments/other_theories_k_pilot.py':digest(Path(__file__))},
                raw_sha256={archive.relative_to(ROOT).as_posix():digest(archive)})
    (ROOT/'experiments/other-theories-k-pilot.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
