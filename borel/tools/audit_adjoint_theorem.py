"""Consolidation audit: fresh integer roots, exact transfers, scheme controls.

No Sage, no imported root-data implementation, and no claim to rerun the
historical ideal-membership calculations. Their current bindings are checked.
The geometric proof remains in ADJOINT_GENERIC_FIBER_THEOREM.md.
"""
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
PROOFS = {
    6: {
        '235': 'CUBIC_LEVI_REFINEMENTS.md',
        '2356': 'CUBIC_LEVI_REFINEMENTS.md',
        '12356': 'BRANCHING_NORM_FIBERS.md',
    },
    7: {
        '125': 'BINARY_BRANCHING_FIBERS.md',
        '1235': 'BINARY_BRANCHING_FIBERS.md',
        '1257': 'E7_SYMMETRIC_PENCIL.md',
        '12457': 'E7_MULTIPLIER_BOUNDARY.md',
        '123567': 'ADJOINT_PENCIL_SCOPE.md',
        '12356': 'E7_FINAL_FAMILIES.md',
        '12357': 'E7_FINAL_FAMILIES.md',
    },
    8: {
        '125': 'E8_A1_TRIPLE.md',
        '1235': 'E8_A2_A1_PAIR.md',
        '1257': 'E8_FINAL_CASES.md',
        '12356': 'E8_FINAL_CASES.md',
        '12357': 'E8_LAST_ADJOINT_FIBER.md',
        '12457': 'E8_FINAL_CASES.md',
        '123468': 'UNIFORM_ADJOINT_PROOF_TARGET.md',
        '123567': 'E8_FINAL_CASES.md',
        '123568': 'UNIFORM_ADJOINT_PROOF_TARGET.md',
        '1235678': 'MULTIPLIER_BRANCHING_FIBERS.md',
    },
}
REPORTS = (
    'e7-companion-membership-verification.json',
    'e7-companion-geometry-audit.json',
    'e8-final-cases-audit.json',
    'e8-last-case-audit.json',
    'e8-last-32003-verification.json',
    'e8-last-32009-verification.json',
    'f4-adjoint-fibers-audit.json',
    'g2-adjoint-fibers-verification.json',
)


def digest(name, raw=False):
    data = (ROOT/name).read_bytes()
    return sha256(data if raw else data.replace(b'\r\n', b'\n')).hexdigest()


def read(name):
    return json.loads((ROOT/name).read_text())


def root_data(rank):
    edges = [(1, 3), (2, 4), (3, 4)] + [(i, i+1) for i in range(4, rank)]
    cartan = [[2*int(i == j) for j in range(rank)] for i in range(rank)]
    for i, j in edges:
        cartan[i-1][j-1] = cartan[j-1][i-1] = -1

    def reflect(root, node):
        i = node-1
        v = list(root)
        v[i] -= sum(cartan[i][j]*root[j] for j in range(rank))
        return tuple(v)

    roots = {tuple(int(j == i) for j in range(rank)) for i in range(rank)}
    pending = list(roots)
    while pending:
        root = pending.pop()
        for node in range(1, rank+1):
            image = reflect(root, node)
            if image not in roots:
                roots.add(image)
                pending.append(image)
    assert len(roots) == {6: 72, 7: 126, 8: 240}[rank]
    assert all(sum(v[i]*cartan[i][j]*v[j] for i in range(rank)
                   for j in range(rank)) == 2 for v in roots)
    return cartan, roots, reflect


def coverage(rank, classification):
    cartan, roots, reflect = root_data(rank)
    subsets = {}
    for mask in range(2**rank-1):
        nodes = tuple(i+1 for i in range(rank) if mask >> i & 1)
        # Every proper E-type induced graph is a disjoint union of paths
        # (all A) unless it retains a trivalent vertex (a D/E factor).
        if any(sum(cartan[i-1][j-1] == -1 for j in nodes) > 2 for i in nodes):
            kind = 'non_A'
        else:
            dominant = []
            for root in roots:
                labels = [sum(cartan[i-1][j]*root[j] for j in range(rank))
                          for i in nodes]
                if all(a >= 0 for a in labels):
                    dominant.append(sum(labels))
            degree = max(dominant, default=0)
            assert degree <= 3
            kind = 'cubic' if degree == 3 else 'quadratic'
        subsets[nodes] = kind
    counts = Counter(subsets.values())
    expected = classification['summary']['E'+str(rank)]['counts']
    assert all(counts[k] == expected[k] for k in ('quadratic', 'cubic', 'non_A'))
    cubic = {nodes for nodes, kind in subsets.items() if kind == 'cubic'}
    archive_cubic = {tuple(c['nodes']) for c in classification['cases']
                     if c['type'] == 'E'+str(rank)}
    assert cubic == archive_cubic

    if rank == 6:
        transfers = classification['E6_Weyl_transfers']
    else:
        filename = ('e7-multiplier-boundary-verification.json' if rank == 7
                    else 'e8-weyl-transfers-verification.json')
        transfers = [dict(seed=c['seed'], target=m['nodes'], word=m['word'])
                     for c in read(filename)['Weyl_classes'] for m in c['members']]

    def subsystem(nodes):
        omitted = set(range(rank)) - {i-1 for i in nodes}
        return {root for root in roots if all(root[j] == 0 for j in omitted)}

    seen = set()
    seeds = Counter()
    for record in transfers:
        seed, target, word = record['seed'], tuple(record['target']), record['word']
        seed_name = ''.join(map(str, seed))
        assert seed_name in PROOFS[rank]
        assert (ROOT/PROOFS[rank][seed_name]).is_file()
        assert target in cubic and target not in seen
        images = subsystem(seed)
        for node in word:
            images = {reflect(root, node) for root in images}
        assert images == subsystem(target), (rank, seed, target)
        seen.add(target)
        seeds[seed_name] += 1
    assert seen == cubic and set(seeds) == set(PROOFS[rank])
    assert sum(counts.values()) == 2**rank-1
    return dict(roots=len(roots), proper_subsets=len(subsets), counts=dict(counts),
                cubic_word_replays=len(transfers), seed_multiplicities=dict(seeds),
                proof_dependencies=PROOFS[rank], complete_standard_subset_coverage=True)


def scheme_controls():
    x, y, z, c = sp.symbols('x y z c')
    field = sp.QQ.frac_field(c)
    # An embedded nilpotent survives a generic target even though the
    # unique minimal component is generically reduced (invert y).
    ideal = sp.groebner([x*x, x*y, z-c], x, y, z, domain=field)
    assert ideal.reduce(x)[1] == x
    assert ideal.reduce(x*x)[1] == ideal.reduce(x*y)[1] == 0
    # An incidence nilpotent supported at target zero is legitimately
    # excluded by generic target localization, with the actual ideal.
    ideal = sp.groebner([x*x, x*y, x*z, z-c], x, y, z, domain=field)
    assert ideal.reduce(x)[1] == 0
    # Smooth excess: rank-two minors of (x,0) vanish identically although
    # the quotient is a smooth polynomial ring in two variables.
    jac = sp.Matrix([[1, 0, 0], [0, 0, 0]])
    assert jac.rank() == 1 and all(jac[:, cols].det() == 0
                                  for cols in ([0, 1], [0, 2], [1, 2]))
    # Affine modular dimension is not a characteristic-zero certificate.
    # x*(p*x-1), y*(p*x-1) has a line over Q but only a point mod p.
    p = 101
    affine = [x*(p*x-1), y*(p*x-1)]
    assert all(f.subs(x, sp.Rational(1, p)) == 0 for f in affine)
    assert sp.groebner(affine, x, y, modulus=p).is_zero_dimensional
    return dict(embedded_nilpotent_detected=True,
                actual_target_localization_eliminates_supported_nilpotent=True,
                smooth_excess_not_rejected_by_redundant_minors=True,
                affine_modular_dimension_shortcut_rejected=True)


def main():
    classification = read('unit-constraint-classification-verification.json')
    cases = {f'E{n}': coverage(n, classification) for n in (6, 7, 8)}
    bindings = {}
    for name in REPORTS:
        report = read(name)
        assert report['result'] == 'passed', name
        for key in ('source_sha256', 'raw_sha256'):
            for path, wanted in report.get(key, {}).items():
                assert digest(path, raw=(key == 'raw_sha256')) == wanted, (name, path)
                bindings[(key, path)] = wanted
    for prime in (32003, 32009):
        r = read(f'e8-last-{prime}-verification.json')
        assert r['original_ideal_memberships_replayed']
        assert r['homogeneous_dimension_upper_bound'] == 4
        assert r['homogeneous_degree_upper_bound'] == r['accounted_infinity_degree'] == 21
        assert r['characteristic_zero_affine_singular_dimension_bound'] == 2
        assert not r['nonemptiness_claim']
    names = set(REPORTS) | {'ADJOINT_GENERIC_FIBER_THEOREM.md',
        'tools/audit_adjoint_theorem.py', 'unit-constraint-classification-verification.json',
        'e7-multiplier-boundary-verification.json', 'e8-weyl-transfers-verification.json',
        'GENERAL_CAYLEY_FIBER_CRITICAL.md', 'QUADRATIC_BRANCHING_FIBERS.md',
        'CLASSICAL_ADJOINT_FIBERS.md', 'EXCEPTIONAL_QUARTIC_REDUCTION.md',
        'EXCEPTIONAL_ADJOINT_CONTACT_FIBERS.md', 'G2_ADJOINT_FIBERS.md',
        'F4_ADJOINT_FIBERS.md'}
    names.update(path for proofs in PROOFS.values() for path in proofs.values())
    report = dict(result='passed', date='2026-10-26', uses_Sage=False,
                  exceptional_E_coverage=cases,
                  proper_E_subsets_independently_enumerated=445,
                  cubic_standard_subset_transfers_replayed=136,
                  reports_bound=len(REPORTS), unique_dependency_bindings_checked=len(bindings),
                  exact_scheme_controls=scheme_controls(),
                  historical_heavy_membership_calculations_rerun=False,
                  full_Weyl_orbits_reenumerated=False, external_peer_review=False,
                  geometric_proof_audit='ADJOINT_GENERIC_FIBER_THEOREM.md',
                  scope='Fresh coverage and full-subsystem word replay; exact control ideals; saved certificate bindings. Written scheme/geometric audit is a separate mathematical argument.',
                  source_sha256={n: digest(n) for n in sorted(names)})
    (ROOT/'adjoint-theorem-audit.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items()
                      if k not in ('source_sha256', 'exceptional_E_coverage')}, indent=2))


if __name__ == '__main__':
    main()
