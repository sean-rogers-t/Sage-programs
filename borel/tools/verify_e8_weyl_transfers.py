"""Exhaustive E8 root-line orbits and explicit proper-Levi proof transfers."""
from array import array
from collections import Counter
from hashlib import sha256
from pathlib import Path
import gzip
import json
from sage.all import CartanType, RootSystem

ROOT = Path(__file__).resolve().parents[1]
PROVED = {
    (1,2,3,5,6,7,8): 'MULTIPLIER_BRANCHING_FIBERS.md: Theorem 3 (node 4 omitted)',
    (1,2,3,4,6,8): 'UNIFORM_ADJOINT_PROOF_TARGET.md: equation (46), entire critical surface smooth',
    (1,2,3,5,6,8): 'UNIFORM_ADJOINT_PROOF_TARGET.md: quasifiniteness theorem and six boundary exclusions',
}


def main():
    A = CartanType('E8').cartan_matrix()
    lattice = RootSystem('E8').root_lattice()
    roots = sorted(tuple(int(r[i+1]) for i in range(8)) for r in lattice.roots())
    positive = sorted(r for r in roots if min(r) >= 0)
    assert len(roots) == 240 and len(positive) == 120
    def canonical(r):
        return r if min(r) >= 0 else tuple(-v for v in r)
    index = {r:i for i,r in enumerate(positive)}
    perms = []
    for i in range(8):
        p = []
        for r in positive:
            v = list(r)
            v[i] -= sum(int(A[i,j])*r[j] for j in range(8))
            p.append(index[canonical(tuple(v))])
        assert sorted(p) == list(range(120))
        perms.append(bytes(p + list(range(120,256))))
    subsets = {}
    shapes = {}
    cubic = set()
    counts = Counter()
    for mask in range(255):
        nodes = tuple(i+1 for i in range(8) if mask >> i & 1)
        I = {i-1 for i in nodes}
        subsystem = bytes(k for k,r in enumerate(positive)
                          if all(r[j] == 0 for j in range(8) if j not in I))
        subsets[subsystem] = nodes
        todo = set(I); components = []
        while todo:
            component = {min(todo)}
            while True:
                grown = component | {j for i in component for j in todo if i != j and A[i,j]}
                if grown == component:
                    break
                component = grown
            todo -= component
            components.append(component)
        non_A = any(sum(i != j and A[i,j] != 0 for j in c) > 2 for c in components for i in c)
        if non_A:
            counts['non_A'] += 1
            continue
        shapes[nodes] = '+'.join('A'+str(len(c)) for c in sorted(components,key=len)) or 'T'
        labels = [tuple(sum(int(A[i,j])*r[j] for j in range(8)) for i in range(8)) for r in roots]
        degree = max(sum(r[i] for i in I) for r in labels if all(r[i] >= 0 for i in I))
        assert degree <= 3
        counts['cubic' if degree == 3 else 'quadratic'] += 1
        if degree == 3:
            cubic.add(nodes)
    assert dict(counts) == dict(quadratic=147,cubic=93,non_A=15)
    archive = json.loads(gzip.decompress((ROOT/'unit-constraint-cases.json.gz').read_bytes()))
    assert cubic == {tuple(r['character']['nodes']) for r in archive if r['character']['type'] == 'E8'}
    node_to_subsystem = {nodes:s for s,nodes in subsets.items()}
    pending = set(cubic); groups = []
    while pending:
        seed = min(pending & PROVED.keys()) if pending & PROVED.keys() else min(pending)
        start = node_to_subsystem[seed]
        states = [start]; seen = {start:0}; parents = array('i',[-1]); steps = bytearray([0])
        members = []
        for k,current in enumerate(states):
            if current in subsets:
                nodes = subsets[current]
                assert nodes in cubic, ('conjugate outside cubic census',seed,nodes)
                word = []; cursor = k
                while parents[cursor] != -1:
                    word.append(steps[cursor]); cursor = parents[cursor]
                word.reverse()
                members.append(dict(nodes=nodes,word=word))
            for i,p in enumerate(perms):
                nxt = bytes(sorted(current.translate(p)))
                if nxt not in seen:
                    seen[nxt] = len(states); states.append(nxt); parents.append(k); steps.append(i+1)
        member_nodes = {tuple(m['nodes']) for m in members}
        assert member_nodes <= pending
        assert {shapes[n] for n in member_nodes} == {shapes[seed]}
        pending -= member_nodes
        groups.append(dict(seed=seed,derived_type=shapes[seed],orbit_size=len(states),
                           members=sorted(members,key=lambda m:m['nodes']),
                           proof=PROVED.get(seed),resolved=seed in PROVED))
        print(json.dumps(dict(seed=seed,shape=shapes[seed],orbit=len(states),
                              standard_members=len(members),resolved=seed in PROVED)),flush=True)
    covered = sum(len(g['members']) for g in groups if g['resolved'])
    result = dict(result='passed',scope='Complex simply connected E8, adjoint Springer map; Weyl conjugacy and transfer of existing generic regular-source proofs only.',
                  counts=dict(counts),cubic_Weyl_classes=len(groups),
                  resolved_cubic_classes=sum(g['resolved'] for g in groups),
                  unresolved_cubic_classes=sum(not g['resolved'] for g in groups),
                  covered_cubic_standard_subsets=covered,
                  unresolved_standard_subsets=93-covered,
                  covered_proper_standard_subsets=162+covered,
                  Weyl_classes=sorted(groups,key=lambda g:(len(g['seed']),g['seed'])),
                  word_convention='Apply listed simple reflections in order, left to right; the represented group element is s_last ... s_first.',
                  exhaustive_orbits=True,standard_subsets_checked=255,
                  all_simple_adjoint_conjecture_proved=False)
    names=['tools/verify_e8_weyl_transfers.py','E8_WEYL_TRANSFERS.md',
           'MULTIPLIER_BRANCHING_FIBERS.md','UNIFORM_ADJOINT_PROOF_TARGET.md']
    result['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                             for n in names}
    result['archive_sha256']=sha256((ROOT/'unit-constraint-cases.json.gz').read_bytes()).hexdigest()
    (ROOT/'e8-weyl-transfers-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('Weyl_classes','source_sha256')},indent=2))

if __name__ == '__main__':
    main()
