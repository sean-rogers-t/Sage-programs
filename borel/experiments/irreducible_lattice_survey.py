"""Closed lattice classification versus full Laurent coefficient lattices."""
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from time import perf_counter

from sage.all import CartanMatrix, QQ
from sage_borel import irreducible_character_lattice, springer_character_map, SpringerCharacterMap

ROOT = Path(__file__).resolve().parents[1]


def rows(A):
    return [[str(x) for x in r] for r in A.rows()]


def cases():
    for kind,n,total in [("G2",2,6),("F4",4,4),("E6",6,3),("E7",7,3),("E8",8,3)]:
        for highest in product(range(total+1),repeat=n):
            if 0 < sum(highest) <= total:
                yield kind,highest
    for kind,weights in [("A1",[(i,) for i in range(1,9)]),
                         ("A2",[(1,0),(2,0),(1,1),(2,1),(3,3)]),
                         ("A3",[(1,0,0),(0,1,0),(0,0,2),(1,1,0),(0,2,0)]),
                         ("B2",[(1,0),(0,1),(2,0),(0,2),(1,1)])]:
        for highest in weights:
            yield kind,highest
    for kind,n in [("G2",2),("F4",4),("E6",6),("E7",7),("E8",8)]:
        yield kind,(100,)+(0,)*(n-1)


def main():
    started = perf_counter()
    data = dict(result="running",cases=[],full_weight_dimension_limit=10000)
    for kind,highest in cases():
        fast = irreducible_character_lattice(kind,highest)
        entry = dict(cartan_type=kind,highest_weight=highest,dimension=str(fast.dimension),
            trace_index=str(fast.trace_index),central_character_order=fast.central_character_order,
            lattice_type=fast.lattice_type,trace_gram=rows(fast.trace_gram),
            moment_basis=rows(fast.moment_lattice_basis),linear_basis=rows(fast.linear_basis_matrix),
            multiplicity_one_paths=fast.multiplicity_one_paths,full_weight_check=False)
        if fast.dimension <= data['full_weight_dimension_limit']:
            full = springer_character_map(kind,highest)
            assert full.dimension == fast.dimension and full.trace_gram == fast.trace_gram
            T = fast.linear_basis_matrix.change_ring(QQ).inverse()*full.linear_basis_matrix
            assert T.denominator() == 1 and abs(T.det()) == 1
            weights = dict(full.weight_multiplicities)
            C = CartanMatrix(kind)
            for path in fast.multiplicity_one_paths:
                mu = tuple(highest)
                assert weights[mu] == 1
                for j in path:
                    mu = tuple(mu[i]-C[i,j-1] for i in range(len(highest)))
                    assert weights[mu] == 1
            entry.update(full_weight_check=True,full_linear_basis=rows(full.linear_basis_matrix),
                         distinct_weights=len(weights))
        data['cases'].append(entry)
        if len(data['cases']) % 50 == 0:
            print('Checked',len(data['cases']),'elapsed',round(perf_counter()-started,2),flush=True)
    minimal = springer_character_map("F4",(0,0,0,1))
    triple = SpringerCharacterMap("F4",[(0,0,0,1)]*3)
    adjoint = springer_character_map("F4",(1,0,0,0))
    assert triple.trace_gram == adjoint.trace_gram
    T = triple.linear_basis_matrix.inverse()*minimal.linear_basis_matrix
    assert T.denominator() == 1 and abs(T.det()) == 1
    T = (3*triple.linear_basis_matrix).inverse()*adjoint.linear_basis_matrix
    assert T.denominator() == 1 and abs(T.det()) == 1
    data['reducible_control'] = dict(type="F4",representation="3*V(omega4)",
        trace_index=18,linear_scale=6,adjoint_trace_index=18,adjoint_linear_scale=18,
        same_normalized_map_as_minimal=True)
    files = [Path(__file__).resolve(),ROOT/'src/sage_borel/character_lattices.py',
        ROOT/'src/sage_borel/characters.py',ROOT/'src/sage_borel/__init__.py',
        ROOT/'IRREDUCIBLE_CHARACTER_LATTICES.md',ROOT/'tests/test_character_lattices.py']
    data.update(result="passed",total_cases=len(data['cases']),
        full_weight_cases=sum(row['full_weight_check'] for row in data['cases']),
        seconds=round(perf_counter()-started,3),
        source_hash_semantics="SHA256 with CRLF normalized to LF.",
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in files})
    (ROOT/'experiments/irreducible-lattice-survey.json').write_text(json.dumps(data,indent=2)+'\n')
    print({k:v for k,v in data.items() if k not in ('cases','source_sha256')},flush=True)


if __name__ == '__main__':
    main()
