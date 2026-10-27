"""Independent E8 audit from the 112 integral and 128 half-integral roots.

No Sage or generator imports. Enumerate full signed-root subsystem orbits,
not the positive-root-line representation used by the certificate generator.
"""
from collections import Counter, deque
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def main():
    source = json.loads((ROOT/'e8-weyl-transfers-verification.json').read_text())
    roots = set()
    for i,j in combinations(range(8),2):
        for a,b in product((-2,2),repeat=2):
            v = [0]*8; v[i] = a; v[j] = b; roots.add(tuple(v))
    for v in product((-1,1),repeat=8):
        if sum(x < 0 for x in v) % 2 == 0:
            roots.add(v)
    roots = sorted(roots)
    assert len(roots) == 240
    simple = [(1,-1,-1,-1,-1,-1,-1,1),(2,2,0,0,0,0,0,0),(-2,2,0,0,0,0,0,0)]
    for i in range(1,6):
        v = [0]*8; v[i] = -2; v[i+1] = 2; simple.append(tuple(v))
    dot = lambda a,b:sum(x*y for x,y in zip(a,b))
    A = [[dot(a,b)//4 for b in simple] for a in simple]
    edges = {(1,3),(2,4),(3,4),(4,5),(5,6),(6,7),(7,8)}
    assert all(A[i][j] == (2 if i == j else -1 if tuple(sorted((i+1,j+1))) in edges else 0)
               for i in range(8) for j in range(8))
    B = s.Matrix.hstack(*(s.Matrix(a) for a in simple)); Bi = B.inv()
    assert abs(B.det()) == 256
    coeff = [tuple(Bi*s.Matrix(r)) for r in roots]
    assert all(all(c.q == 1 for c in row) and (min(row)>=0 or max(row)<=0) for row in coeff)
    index = {r:i for i,r in enumerate(roots)}
    perms = []
    for a in simple:
        p = []
        for r in roots:
            assert dot(r,a) % 4 == 0
            v = tuple(x-dot(r,a)//4*y for x,y in zip(r,a))
            p.append(index[v])
        assert sorted(p) == list(range(240)) and all(p[p[i]]==i for i in range(240))
        perms.append(bytes(p+list(range(240,256))))
    subsets = {}; cubic = set(); counts = Counter()
    for mask in range(255):
        nodes = tuple(i+1 for i in range(8) if mask>>i&1)
        I = {i-1 for i in nodes}
        sub = bytes(k for k,c in enumerate(coeff) if all(c[j]==0 for j in range(8) if j not in I))
        subsets[sub] = nodes
        if any(sum(A[i][j] == -1 for j in I) > 2 for i in I):
            counts['non_A'] += 1
        else:
            labels = [tuple(dot(r,a)//4 for a in simple) for r in roots]
            degree = max(sum(r[i] for i in I) for r in labels if all(r[i]>=0 for i in I))
            assert degree <= 3
            counts['cubic' if degree==3 else 'quadratic'] += 1
            if degree == 3:
                cubic.add(nodes)
    assert dict(counts) == source['counts']
    inverse = {v:k for k,v in subsets.items()}
    covered = set(); all_members = set(); orbit_sizes = []; transfer_words = 0
    proof_seeds = {(1,2,3,4,6,8),(1,2,3,5,6,8),(1,2,3,5,6,7,8)}
    for group in source['Weyl_classes']:
        seed = tuple(group['seed']); start = inverse[seed]
        seen = {start}; queue = deque([start]); actual = set()
        while queue:
            cur = queue.popleft()
            if cur in subsets:
                actual.add(subsets[cur])
            for p in reversed(perms):
                nxt = bytes(sorted(cur.translate(p)))
                if nxt not in seen:
                    seen.add(nxt); queue.append(nxt)
        assert len(seen) == group['orbit_size']
        assert actual == {tuple(m['nodes']) for m in group['members']}
        assert actual <= cubic and not actual & all_members
        all_members |= actual
        assert group['resolved'] == (seed in proof_seeds)
        if group['resolved']:
            assert group['proof']
            covered |= actual
        for member in group['members']:
            cur = start
            # Also propagate the entire simple basis, independently of root IDs.
            basis = list(simple)
            for node in member['word']:
                p = perms[node-1]; cur = bytes(sorted(cur.translate(p)))
                a = simple[node-1]
                basis = [tuple(x-dot(v,a)//4*y for x,y in zip(v,a)) for v in basis]
            assert cur == inverse[tuple(member['nodes'])]
            W = s.Matrix.hstack(*(s.Matrix(v) for v in basis))*Bi
            assert W.T*W == s.eye(8) and W.det() in (-1,1)
            integral = Bi*W*B
            assert all(x.q == 1 for x in integral)
            target_I = {n-1 for n in member['nodes']}
            assert all(integral[j,i-1] == 0 for i in seed for j in range(8) if j not in target_I)
            # Orthogonal complement = Lie center; integral isometry already
            # proves its exact transfer, without choosing an isogeny cover.
            transfer_words += 1
        orbit_sizes.append(len(seen))
        print(json.dumps(dict(shape=group['derived_type'],full_signed_root_orbit=len(seen),
                              standard_members=len(actual),all_words_checked=True)),flush=True)
    assert all_members == cubic and len(all_members)==93
    assert len(covered)==source['covered_cubic_standard_subsets']==7
    assert source['unresolved_standard_subsets']==86
    assert source['cubic_Weyl_classes']==10 and source['unresolved_cubic_classes']==7
    assert source['covered_proper_standard_subsets']==169
    # Source proof scope checked separately in the written geometric argument.
    result = dict(result='passed',independent_model='Euclidean E8: 112 integral plus 128 half-integral roots, scaled by 2; full signed-root orbit BFS',
                  independent_census=dict(counts),full_orbits_reenumerated=len(orbit_sizes),
                  full_orbit_sizes=orbit_sizes,explicit_words_checked=transfer_words,
                  root_lattice_and_Killing_isometries=True,Lie_center_transfer_checked=True,
                  resolved_classes=3,unresolved_classes=7,covered_standard_subsets=169,
                  unresolved_standard_subsets=86,all_simple_adjoint_conjecture_proved=False,
                  scope='Independent conjugacy/count/isometry audit; does not rerun the three seed polynomial proofs or establish smoothness for the seven remaining classes.')
    names=['tools/audit_e8_weyl_transfers.py','E8_WEYL_TRANSFERS.md',
           'e8-weyl-transfers-verification.json']
    result['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                             for n in names}
    (ROOT/'e8-weyl-transfers-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))

if __name__=='__main__':
    main()
