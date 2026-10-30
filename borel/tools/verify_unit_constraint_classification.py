"""Exhaustive E6/E7/E8 cubic all-A branching and unit-constraint census."""
from collections import Counter, deque
from hashlib import sha256
from pathlib import Path
import json
import gzip
from sage.all import QQ, LaurentPolynomialRing, PolynomialRing, prod, CartanType, RootSystem, vector
from verify_affine_branching_classification import levi_data, subsets
from verify_branching_norm_fibers import charge_polynomial
from linear_character_elimination import eliminate, as_laurent
from unit_constraint_systems import classify

ROOT=Path(__file__).resolve().parents[1]


def e6_transfers(rows):
    A=CartanType('E6').cartan_matrix();Ai=A.inverse()
    ambient=RootSystem('E6').ambient_space();cr=ambient.simple_coroots()
    roots=[tuple(int(r.scalar(cr[j+1])) for j in range(6)) for r in ambient.roots()]
    def subsystem(nodes):
        outside=set(range(6))-{i-1 for i in nodes}
        return tuple(sorted(r for r in roots if all((Ai*vector(r))[i]==0 for i in outside)))
    seeds={(1,1,1):(2,3,5),(1,1,2):(2,3,5,6),(1,2,2):(1,2,3,5,6)}
    answer=[]
    for shape,seed in seeds.items():
        targets={subsystem(row['character']['nodes']):row['character']['nodes'] for row in rows
                 if row['character']['type']=='E6' and
                 tuple(sorted(len(o) for o in row['character']['factor_orders']))==shape}
        start=subsystem(seed);queue=deque([(start,())]);seen={start}
        while targets and queue:
            current,word=queue.popleft()
            if current in targets:
                answer.append(dict(seed=list(seed),target=targets.pop(current),word=list(word)))
            for i in range(6):
                nxt=tuple(sorted(tuple(r[j]-r[i]*int(A[j,i]) for j in range(6)) for r in current))
                if nxt not in seen:seen.add(nxt);queue.append((nxt,word+(i+1,)))
        assert not targets
    assert len(answer)==11
    return sorted(answer,key=lambda r:r['target'])


def residual(row,choice):
    n=len(row['nodes']);s=len(row['central_cocharacters'])
    R=LaurentPolynomialRing(QQ,['z'+str(i) for i in range(s)]);zs=R.gens()
    P=PolynomialRing(R.fraction_field(),['x'+str(i) for i in range(n)]);xs=P.gens()
    phi=sum(t['coefficient']*prod(z**e for z,e in zip(zs,t['charge']))*
            prod(x**e for x,e in zip(xs,t['exponent'])) for t in row['terms'])
    if choice:
        cert=eliminate(phi,[xs[i] for i in choice['multipliers']],
                       [xs[i] for i in choice['paired']],R)
        reduced=cert['reduced'];graph=cert['graph']
        remain=[i for i in range(n) if xs[i] not in graph]
        assert all(phi.derivative(xs[i]).subs(graph)==reduced.derivative(xs[i]) for i in remain)
        for z in zs:
            original=P({e:z*v.derivative(z) for e,v in phi.dict().items()})
            after=P({e:z*v.derivative(z) for e,v in reduced.dict().items()})
            assert original.subs(graph)==after
    else:
        # All A1: keep the canonical first square, not a claimed multiplier.
        assert all(f=='A1' for f in row['factors'])
        x=xs[0];coefficient=phi.derivative(x).derivative(x)/2
        assert coefficient==1
        value=-phi.derivative(x).subs({x:0})/2
        graph={x:value};reduced=phi.subs(graph);remain=list(range(1,n))
        assert phi.derivative(x).subs(graph)==0
        assert all(phi.derivative(xs[i]).subs(graph)==reduced.derivative(xs[i]) for i in remain)
        for z in zs:
            original=P({e:z*v.derivative(z) for e,v in phi.dict().items()})
            after=P({e:z*v.derivative(z) for e,v in reduced.dict().items()})
            assert original.subs(graph)==after
    records=[]
    for e,coef in sorted(reduced.dict().items()):
        for charge,c in sorted(as_laurent(R,coef).dict().items()):
            # Explicit records also retain affine-independent terms.
            if s==1:charge=(charge,)
            records.append(dict(exponent=[int(e[i]) for i in remain],
                                charge=list(map(int,charge)),coefficient=str(c)))
    return dict(kind='simultaneous_multipliers' if choice else 'one_A1_square',
                remaining_nodes=[row['nodes'][i] for i in remain],
                degree=int(reduced.total_degree()),terms=records,
                original_target_derivatives_verified=True)


def main():
    rows=[];summary={}
    for kind in ('E6','E7','E8'):
        counts=Counter();degrees=Counter();ranks=Counter()
        for I in subsets(int(kind[1:])):
            factors,_=levi_data(kind,I)
            if any(not name.startswith('A') for name,order in factors):
                counts['non_A']+=1;continue
            counts['all_A']+=1
            row=charge_polynomial(kind,I)
            shape=tuple(sorted(len(order) for name,order in factors))
            if kind=='E6' and shape in ((1,1,1),(1,1,2),(1,2,2)):
                assert row['max_degree']==3
            if row['max_degree']<=2:counts['quadratic']+=1;continue
            assert row['max_degree']==3
            counts['cubic']+=1
            result=classify(row)
            assert bool(result['chosen'])==any(len(order)>=2 for name,order in factors)
            assert (result['maximum_pairs']==2)==(shape in ((1,2,3),(1,2,4)))
            assert result['maximum_pairs']<=2
            red=residual(row,result['chosen']);degrees[red['degree']]+=1;ranks[result['maximum_pairs']]+=1
            assert red['degree']==4 and len(red['remaining_nodes'])<=4
            rows.append(dict(character=row,constraints=result,residual=red))
        summary[kind]=dict(counts=counts,maximum_pair_counts=dict(ranks),residual_degrees=dict(degrees))
        print(kind,json.dumps(summary[kind]),flush=True)
    archive='unit-constraint-cases.json.gz'
    payload=gzip.compress((json.dumps(rows,separators=(',',':'))+'\n').encode(),mtime=0)
    (ROOT/archive).write_bytes(payload)
    report=dict(result='passed',summary=summary,case_archive=archive,
                case_archive_sha256=sha256(payload).hexdigest(),
                E6_Weyl_transfers=e6_transfers(rows),
                all_E6_proper_Levi_generic_regular_fibers_smooth_or_empty=True,
                cases=[dict(type=r['character']['type'],nodes=r['character']['nodes'],
                            factors=r['character']['factors'],
                            maximum_pairs=r['constraints']['maximum_pairs'],
                            remaining_nodes=r['residual']['remaining_nodes'],
                            residual_degree=r['residual']['degree']) for r in rows],
                original_coordinate_exhaustive=True,recursive_reduction_classified=False,
                all_adjoint_generic_smoothness_proved=False,literature_priority_claim=False)
    names=['EXCEPTIONAL_QUARTIC_REDUCTION.md','tools/unit_constraint_systems.py',
           'tools/verify_unit_constraint_classification.py',
           'tools/linear_character_elimination.py','tools/verify_branching_norm_fibers.py',
           'QUADRATIC_BRANCHING_FIBERS.md','EXCEPTIONAL_ADJOINT_CONTACT_FIBERS.md',
           'CUBIC_LEVI_REFINEMENTS.md','BRANCHING_NORM_FIBERS.md']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'unit-constraint-classification-verification.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
