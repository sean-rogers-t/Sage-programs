"""Check the equivariant K determinant via two independent coordinate methods.

Fixed-point evaluations use Weyl orbits; Euler coordinates use Laurent
Demazure strings. Reflection exponents use small weight-orbit BFS rather
than enumerating the exceptional Weyl groups.
"""
from hashlib import sha256
import json
from pathlib import Path
from sage.all import QQ,RootSystem,matrix,prod,vector
from sage_borel import BorelMap,SpringerCharacterMap
from other_theories_k_pilot import demazure

ROOT=Path(__file__).resolve().parents[1]


def digest(p):return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def root_data(kind):
    lattice=RootSystem(kind).weight_lattice();r=lattice.dimension()
    roots=[vector(QQ,[lattice.simple_roots()[i+1].monomial_coefficients().get(j+1,0) for j in range(r)]) for i in range(r)]
    return lattice,roots


def inverse_word_action(v,word,roots):
    v=vector(QQ,v)
    for node in word:v-=v[node-1]*roots[node-1]
    return v


def character_value(weight,point):return prod(t**int(a) for t,a in zip(point,weight))


def run_case(primary,cert,adjoint):
    kind=primary['cartan_type'];b=BorelMap(kind);r=b.rank
    lattice,roots=root_data(kind);words=[tuple(w) for w in cert['words']]
    elements=[b.W.from_reduced_word(w) for w in words];simple=b.W.simple_reflections()
    H=matrix(QQ,cert['H_rows']);pivots=H.pivots();Hbasis=H.matrix_from_columns(pivots)
    exponents=[cert['exponents'][i] for i in pivots]
    B=matrix(QQ,primary['linear_lattice_rows'])
    inputs=[prod(sum(B[i,j]*b.weights[i] for i in range(r))**p for j,p in enumerate(powers)) for powers in exponents]
    steinberg=[]
    for w,word in zip(elements,words):
        weight=vector(QQ,[int((simple[i+1]*w).length()<w.length()) for i in range(r)])
        steinberg.append(-inverse_word_action(weight,word,roots))
    positive=[]
    for a in lattice.positive_roots():
        positive.append(vector(QQ,[a.monomial_coefficients().get(i+1,0) for i in range(r)]))
    maps=[SpringerCharacterMap(kind,tuple(tuple(x) for x in primary['highest_weights'])),
          SpringerCharacterMap(kind,(adjoint,))]
    reports=[]
    for point0 in ((2,3),(3,5),(2,5)):
        point=tuple(map(QQ,point0));delta=prod(point)*prod(1-character_value(-a,point) for a in positive)
        restrictions=matrix(QQ,[[character_value(inverse_word_action(a,w,roots),point) for a in steinberg] for w in words])
        assert abs(restrictions.det()/delta**(len(words)//2))==1
        values=[];dets=[];euler_dets=[]
        for ch in maps:
            theta=vector(QQ,[ch.pullback(x)(*point) for x in ch.ring.gens()])
            jac=prod(sum(a[i]*theta[i] for i in range(r)) for a in positive)/delta
            evaluation=matrix(QQ,[[prod(sum(inverse_word_action(B.column(j),w,roots)[i]*theta[i] for i in range(r))**p for j,p in enumerate(powers)) for powers in exponents] for w in words])
            assert abs(evaluation.det()/restrictions.det()/Hbasis.det()/jac**(len(words)//2))==1
            columns=[]
            for f in inputs:
                pull=ch.pullback(ch.ring(f));terms={tuple(map(int,p)):c for p,c in pull.dict().items()}
                column=[]
                for word in words:
                    value=terms
                    for node in reversed(word):value=demazure(value,node-1,roots[node-1])
                    column.append(sum(c*character_value(p,point) for p,c in value.items()))
                columns.append(column)
            euler_det=matrix(QQ,columns).transpose().det()
            values.append(jac);dets.append(evaluation.det()/restrictions.det());euler_dets.append(euler_det)
        ratio=(values[1]/values[0])**(len(words)//2)
        assert dets[1]/dets[0]==euler_dets[1]/euler_dets[0]==ratio
        reports.append(dict(point=list(map(str,point)),minimal_root_quotient=str(values[0]),
                            adjoint_root_quotient=str(values[1]),
                            fixed_point_vs_Steinberg_determinant_verified=True,
                            independent_Demazure_determinant_ratio_verified=True))
    return dict(cartan_type=kind,flag_rank=len(words),points=reports)


def orbit_exponent(kind,node):
    lattice,roots=root_data(kind);r=len(roots);start=tuple(int(i==node-1) for i in range(r))
    seen={start};todo=[start]
    while todo:
        a=todo.pop()
        for i in range(r):
            v=vector(QQ,a)-a[i]*roots[i];t=tuple(map(int,v))
            if t not in seen:seen.add(t);todo.append(t)
    ambient=RootSystem(kind).ambient_space();simple=ambient.simple_roots()
    lengths={i:simple[i+1].inner_product(simple[i+1]) for i in range(r)}
    result={}
    for length in set(lengths.values()):
        nodes=[i for i in range(r) if lengths[i]==length]
        counts={sum(a[i]!=0 for a in seen)//2 for i in nodes}
        assert len(counts)==1
        result[str(length)]=counts.pop()
    return dict(cartan_type=kind,omitted_node=node,cosets=len(seen),moved_pair_exponents=result)


def main():
    primary=json.loads((ROOT/'experiments/other-theories-k-pilot.json').read_text())
    cert=json.loads((ROOT/'experiments/other-theories-k-certificates.json').read_text())
    cases=[]
    for kind,adj in [('A2',(1,1)),('C2',(2,0)),('G2',(0,1))]:
        i=next(i for i,c in enumerate(primary['cases']) if c['cartan_type']==kind)
        cases.append(run_case(primary['cases'][i],cert['cases'][i],adj))
        print(kind,'equivariant determinant independently checked',flush=True)
    orbits=[orbit_exponent(kind,node) for kind,node in [('G2',1),('G2',2),('F4',1),('F4',4),('E6',1),('E7',7),('E8',8)]]
    assert orbits[4]['cosets']==27 and set(orbits[4]['moved_pair_exponents'].values())=={6}
    assert orbits[5]['cosets']==56 and set(orbits[5]['moved_pair_exponents'].values())=={12}
    assert orbits[6]['cosets']==240 and set(orbits[6]['moved_pair_exponents'].values())=={57}
    report=dict(result='passed',date='2026-10-22',rank_two_determinants=cases,independent_coset_orbit_counts=orbits,
                scope='Finite exact characteristic-zero determinant corroboration and independent exceptional orbit counts. The all-parabolic theorem is the separate discriminant/localization proof; no integral exceptional Smith classification follows.',
                source_sha256={'experiments/other_theories_equivariant_k_pilot.py':digest(Path(__file__)),
                               'experiments/other_theories_k_pilot.py':digest(ROOT/'experiments/other_theories_k_pilot.py')},
                raw_sha256={name:digest(ROOT/name) for name in ('experiments/other-theories-k-pilot.json','experiments/other-theories-k-certificates.json')})
    (ROOT/'experiments/other-theories-equivariant-k-pilot.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
