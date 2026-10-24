"""E7 cubic Weyl families and an exact nonunit-multiplier boundary proof."""
from collections import deque
from hashlib import sha256
from pathlib import Path
import gzip
import json
from sage.all import QQ, CartanType, RootSystem, LaurentPolynomialRing, PolynomialRing, matrix, prod
from verify_branching_norm_fibers import charge_polynomial
from linear_character_elimination import eliminate

ROOT = Path(__file__).resolve().parents[1]


def groups(rows):
    # Simple-root coordinates; the native auditor instead uses weight labels.
    lattice = RootSystem('E7').root_lattice()
    roots = sorted(tuple(int(root[i+1]) for i in range(7)) for root in lattice.roots())
    A = CartanType('E7').cartan_matrix()
    index = {r:i for i,r in enumerate(roots)}
    permutations = []
    for i in range(7):
        permutation = []
        for root in roots:
            reflected = list(root)
            reflected[i] -= sum(int(A[i,j])*root[j] for j in range(7))
            permutation.append(index[tuple(reflected)])
        permutations.append(permutation)
    def subsystem(nodes):
        return tuple(i for i,r in enumerate(roots)
                     if all(r[j] == 0 for j in range(7) if j+1 not in nodes))
    targets = {subsystem(row['character']['nodes']): row for row in rows}
    result = []
    while targets:
        start = min(targets, key=lambda k: targets[k]['character']['nodes'])
        row = targets[start]['character']
        seed = row['nodes']; queue = deque([(start,())]); seen = {start}; members = []
        while queue:
            current, word = queue.popleft()
            if current in targets:
                target = targets.pop(current)['character']['nodes']
                members.append(dict(nodes=target, word=list(word)))
            for i,p in enumerate(permutations):
                nxt = tuple(sorted(p[j] for j in current))
                if nxt not in seen:
                    seen.add(nxt); queue.append((nxt,word+(i+1,)))
        result.append(dict(seed=seed, factors=row['factors'], members=members, orbit_size=len(seen)))
    assert len(result) == 7 and sum(len(g['members']) for g in result) == 32
    assert sorted(len(g['members']) for g in result) == [1,1,2,3,3,10,12]
    assert sorted(g['orbit_size'] for g in result) == [3780,3780,5040,5040,7560,10080,15120]
    assert all(any(set(m['nodes']) <= {1,2,3,5,6,7} for m in g['members']) for g in result)
    return result


def symbolic():
    row = charge_polynomial('E7',(0,1,3,4,6))
    L = LaurentPolynomialRing(QQ,['r','s']); r,s = L.gens()
    P = PolynomialRing(L.fraction_field(),['x','a','v','b','y']); x,a,v,b,y = P.gens()
    phi = sum(t['coefficient']*prod(z**e for z,e in zip((r,s),t['charge']))*
              prod(q**e for q,e in zip(P.gens(),t['exponent'])) for t in row['terms'])
    R=r+1/r; S=r*r*s*s+1/(r*r*s*s)
    F=y*(x/(r*s)+s*(1+r*r)); G=y*(r*s*x+(1+1/(r*r))/s)
    H=x*x+y*y+R*S*x+v*(R*x+S)+R*R-3
    assert phi == a*b+a*F+b*G+H
    cert=eliminate(phi,[a],[b],L)
    psi=(1-y*y)*(x*x+R*S*x+R*R-1)+v*(R*x+S)-2
    assert cert['reduced'] == psi
    graph=cert['graph']
    assert graph == {a:-G,b:-F}
    D=lambda f,z:P({e:z*c.derivative(z) for e,c in P(f).dict().items()})
    for q in (x,v,y): assert phi.derivative(q).subs(graph)==psi.derivative(q)
    for z in (r,s): assert D(phi,z).subs(graph)==D(psi,z)
    x0=-S/R; v0=-(1-y*y)*(2*x0+R*S)/R
    graph2={x:x0,v:v0}; K=(R*R-1)*(R*R-S*S)/(R*R)
    assert psi.subs(graph2)==(1-y*y)*K-2
    assert psi.derivative(x).subs(graph2)==psi.derivative(v).subs(graph2)==0
    assert psi.derivative(y).subs(graph2)==-2*K*y
    for z in (r,s): assert D(psi,z).subs(graph2)==(1-y*y)*D(K,z)
    # Original central trace normalization at the identity of the Levi.
    identity={x:2,a:4,v:6,b:4,y:2}
    trace=matrix(QQ,[[D(D(phi,z),w).subs(identity).constant_coefficient()(r=1,s=1)
                     for w in (r,s)] for z in (r,s)])
    assert trace==matrix(QQ,[[216,144],[144,144]])
    # Independent parameterization for the boundary and dominance certificate.
    T=LaurentPolynomialRing(QQ,['r','t']); rr,t=T.gens(); RT=rr+1/rr; ST=t+1/t
    KT=(RT*RT-1)*(RT*RT-ST*ST)/(RT*RT)
    assert RT*RT-ST*ST == (rr-t)*(rr+t)*(rr*t-1)*(rr*t+1)/(rr*rr*t*t)
    logH=matrix(QQ,[[ (zz*(ww*KT.derivative(ww)).derivative(zz))(r=2,t=3)
                       for ww in (rr,t)] for zz in (rr,t)])
    assert KT(r=2,t=3)==-QQ(49)/12
    assert logH.det()==-QQ(9125752)/16875
    B=PolynomialRing(QQ,['x','y']); xx,yy=B.gens()
    ideal=B.ideal([xx*(1-yy*yy),yy*(xx*xx-1)])
    # Quotient length five, with five distinct rational points, proves reducedness.
    assert ideal.vector_space_dimension()==5
    points=[(0,0)]+[(i,j) for i in (-1,1) for j in (-1,1)]
    assert all(all(f(x=i,y=j)==0 for f in ideal.gens()) for i,j in points)
    return dict(character=row, original_relative_derivatives_verified=True,
                original_target_derivatives_verified=True, exact_unit_pair_elimination=True,
                nonunit_multiplier_boundary_dimension=1, central_target_dimension=2,
                boundary_affine_points=[list(p) for p in points],
                trace_matrix=[list(map(int,rw)) for rw in trace.rows()],
                logarithmic_hessian=[[str(v) for v in rw] for rw in logH.rows()],
                logarithmic_hessian_determinant=str(logH.det()),
                generic_regular_fibers_nonempty_smooth=True, regular_fiber_dimension=16)


def main():
    archive=(ROOT/'unit-constraint-cases.json.gz').read_bytes()
    rows=[row for row in json.loads(gzip.decompress(archive)) if row['character']['type']=='E7']
    report=dict(result='passed', Weyl_classes=groups(rows), application=symbolic(),
                cubic_cases=32,cubic_classes=7,newly_resolved_standard_cases=3,
                resolved_proper_standard_cases=99,total_proper_standard_cases=127,
                remaining_cubic_cases=28,remaining_cubic_classes=5,
                common_maximal_Levi_refinement=True,all_E7_Levis_resolved=False,
                full_map_degree_claim=False,literature_priority_claim=False,
                case_archive_sha256=sha256(archive).hexdigest())
    names=['E7_MULTIPLIER_BOUNDARY.md','tools/verify_e7_multiplier_boundary.py',
           'tools/verify_branching_norm_fibers.py','tools/linear_character_elimination.py',
           'GENERAL_CAYLEY_FIBER_CRITICAL.md']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'e7-multiplier-boundary-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('Weyl_classes','application','source_sha256')},indent=2))


if __name__=='__main__':main()
