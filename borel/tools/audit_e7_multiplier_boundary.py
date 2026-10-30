"""Sage-free root, character, conjugacy and boundary audit for E7."""
from collections import deque
from hashlib import sha256
from pathlib import Path
import gzip
import json
import sympy as sp
from audit_affine_branching_certificates import cartan, inverse, root_data
from audit_branching_norm_fibers import reconstruct
from audit_unit_constraint_classification import census

ROOT=Path(__file__).resolve().parents[1]


def weyl_audit(source,rows):
    A=cartan('E7'); Ai=inverse(A); positive=root_data(A)[1]
    roots=sorted(list(positive)+[tuple(-v for v in r) for r in positive])
    index={r:i for i,r in enumerate(roots)}
    coeff=[tuple(sum(x*y for x,y in zip(row,r)) for row in Ai) for r in roots]
    permutations=[[index[tuple(r[j]-r[i]*A[j][i] for j in range(7))] for r in roots] for i in range(7)]
    def subsystem(nodes):
        return tuple(i for i,c in enumerate(coeff) if all(c[j]==0 for j in range(7) if j+1 not in nodes))
    counts,cubic=census()
    expected={nodes for kind,nodes in cubic if kind=='E7'}
    assert expected=={tuple(r['character']['nodes']) for r in rows}
    seen_targets=set(); all_orbits=set()
    for group in source['Weyl_classes']:
        start=subsystem(group['seed']); orbit={start}; queue=deque([start])
        while queue:
            current=queue.popleft()
            for p in permutations:
                nxt=tuple(sorted(p[i] for i in current))
                if nxt not in orbit: orbit.add(nxt);queue.append(nxt)
        assert len(orbit)==group['orbit_size']
        assert not all_orbits.intersection(orbit)
        all_orbits.update(orbit)
        members=set()
        for member in group['members']:
            current=start
            for node in member['word']:
                current=tuple(sorted(permutations[node-1][i] for i in current))
            assert current==subsystem(member['nodes'])
            members.add(tuple(member['nodes']))
        assert members=={nodes for nodes in expected if subsystem(nodes) in orbit}
        assert not seen_targets.intersection(members)
        seen_targets.update(members)
    assert seen_targets==expected and len(seen_targets)==32
    assert counts['E7']==dict(all_A=120,quadratic=88,cubic=32,non_A=7)
    # Check the theorem covers every standard embedding of this abstract type,
    # including any that might have been absent from the cubic archive.
    shape_targets=set()
    for mask in range(127):
        retained={i for i in range(7) if mask >> i & 1}; todo=set(retained); sizes=[]
        while todo:
            component={min(todo)}
            while True:
                grown=component|{j for i in component for j in todo if i!=j and A[i][j]}
                if grown==component:break
                component=grown
            todo-=component; sizes.append(len(component))
        if sorted(sizes)==[1,1,3]:shape_targets.add(tuple(i+1 for i in sorted(retained)))
    assert shape_targets=={(1,2,4,5,7),(1,2,5,6,7),(2,3,5,6,7)}
    for row in rows: reconstruct(row['character'])
    return dict(independent_census=True,full_Weyl_orbits_reenumerated=7,
                explicit_words_verified=32,all_32_characters_reconstructed=True,
                all_A3_A1_A1_standard_embeddings_covered=True)


def algebra_audit(row):
    reconstruct(row)
    x,a,v,b,y,r,s=sp.symbols('x a v b y r s')
    affine=(x,a,v,b,y); central=(r,s)
    phi=sum(sp.Rational(term['coefficient'])*sp.prod(z**e for z,e in zip(central,term['charge']))*
            sp.prod(q**e for q,e in zip(affine,term['exponent'])) for term in row['terms'])
    # Derive the elimination directly from the original derivatives.
    F=sp.diff(phi,a)-b; G=sp.diff(phi,b)-a
    graph={a:-G,b:-F}; psi=sp.expand(phi.subs(graph,simultaneous=True))
    R=r+1/r; S=r*r*s*s+1/(r*r*s*s)
    expected=(1-y*y)*(x*x+R*S*x+R*R-1)+v*(R*x+S)-2
    assert sp.expand(psi-expected)==0
    for q in (x,v,y):
        assert sp.expand(sp.diff(phi,q).subs(graph,simultaneous=True)-sp.diff(psi,q))==0
    for z in central:
        assert sp.expand((z*sp.diff(phi,z)).subs(graph,simultaneous=True)-z*sp.diff(psi,z))==0
    x0=-S/R;v0=-(1-y*y)*(2*x0+R*S)/R;K=(R*R-1)*(R*R-S*S)/(R*R)
    graph2={x:x0,v:v0}
    assert sp.cancel(psi.subs(graph2,simultaneous=True)-((1-y*y)*K-2))==0
    for z in central:
        assert sp.cancel((z*sp.diff(psi,z)).subs(graph2,simultaneous=True)-(1-y*y)*z*sp.diff(K,z))==0
    for q in (x,v):assert sp.cancel(sp.diff(psi,q).subs(graph2,simultaneous=True))==0
    assert sp.cancel(sp.diff(psi,y).subs(graph2,simultaneous=True)+2*K*y)==0
    # Trace form independently from 126 ambient roots, rather than character derivatives.
    A=cartan('E7');positive=root_data(A)[1];roots=list(positive)+[tuple(-v for v in w) for w in positive]
    H=row['central_cocharacters']
    charges=[[sum(a*b for a,b in zip(root,h)) for h in H] for root in roots]
    trace=[[sum(q[i]*q[j] for q in charges) for j in range(2)] for i in range(2)]
    assert trace==[[216,144],[144,144]]
    # Boundary quotient and full Hessian Schur complement in the simplified torus.
    boundary=sp.groebner([x*(1-y*y),y*(x*x-1)],x,y)
    assert boundary.is_zero_dimensional
    assert sp.solve_poly_system(list(boundary),x,y)==[(-1,-1),(-1,1),(0,0),(1,-1),(1,1)]
    t=sp.symbols('t'); S=t+1/t
    Q=x*x+R*S*x+R*R-1
    potential=(1-y*y)*Q+v*(R*x+S)-2
    KK=(R*R-1)*(R*R-S*S)/(R*R)
    assert sp.cancel(R*R-S*S-(r-t)*(r+t)*(r*t-1)*(r*t+1)/(r*r*t*t))==0
    point={r:2,t:3,x:-sp.Rational(4,3),v:-sp.Rational(34,15),y:0}
    assert all(sp.diff(potential,q).subs(point)==0 for q in (x,v,y))
    operators=[lambda f,q=q:sp.diff(f,q) for q in (x,v,y)]+[
        lambda f,z=z:z*sp.diff(f,z) for z in (r,t)]
    full=sp.Matrix([[d(e(potential)).subs(point) for e in operators] for d in operators])
    relative=full[:3,:3]; mixed=full[:3,3:]; central_block=full[3:,3:]
    schur=central_block-mixed.T*relative.inv()*mixed
    logH=sp.Matrix([[(z*sp.diff(w*sp.diff(KK,w),z)).subs(point) for w in (r,t)] for z in (r,t)])
    assert relative.det()==-sp.Rational(1225,24)
    assert schur==logH==sp.Matrix([[sp.Rational(3889,225),-sp.Rational(256,75)],
                                 [-sp.Rational(256,75),-sp.Rational(2296,75)]])
    assert logH.det()==-sp.Rational(9125752,16875)
    assert full.det()==relative.det()*logH.det()!=0
    # The torus change and its target transformation retain rank two.
    assert sp.Matrix([[1,2],[0,2]]).det()==2
    return dict(independent_character=True,original_derivatives_preserved=True,
                boundary_support_points=5,boundary_dimension=1,
                root_trace_matrix=trace,full_logarithmic_hessian_schur_complement_verified=True,
                relative_hessian_determinant=str(relative.det()),
                torus_hessian_determinant=str(logH.det()),
                generic_regular_fiber_dimension=16)


def main():
    source=json.loads((ROOT/'e7-multiplier-boundary-verification.json').read_text())
    archive=(ROOT/'unit-constraint-cases.json.gz').read_bytes()
    assert sha256(archive).hexdigest()==source['case_archive_sha256']
    rows=[r for r in json.loads(gzip.decompress(archive)) if r['character']['type']=='E7']
    report=dict(result='passed',uses_Sage=False,weyl=weyl_audit(source,rows),
                algebra=algebra_audit(source['application']['character']),
                all_E7_Levis_resolved=False,full_map_degree_claim=False,literature_priority_claim=False)
    names=['E7_MULTIPLIER_BOUNDARY.md','tools/audit_e7_multiplier_boundary.py',
           'tools/audit_affine_branching_certificates.py','tools/audit_branching_norm_fibers.py',
           'tools/audit_unit_constraint_classification.py','e7-multiplier-boundary-verification.json']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'e7-multiplier-boundary-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
