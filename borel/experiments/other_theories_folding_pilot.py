"""Finite folding tests: trace-map compatibility and exact integral source defect."""
from hashlib import sha256
import json
from pathlib import Path
from sage.all import QQ,ZZ,RootSystem,matrix,vector
from sage_borel import SpringerCharacterMap

ROOT=Path(__file__).resolve().parents[1]


def simple_matrix(kind):
    W=RootSystem(kind).weight_lattice();r=W.dimension()
    return matrix(QQ,[[W.simple_roots()[j+1].monomial_coefficients().get(i+1,0) for j in range(r)] for i in range(r)])


def run(G,H,root_nodes,high_G,high_H,copies,trivials,expected):
    g=SpringerCharacterMap(G,high_G);h=SpringerCharacterMap(H,(high_H,))
    A=matrix(QQ,h.rank,g.rank,lambda i,j:int(root_nodes[j]==i+1))
    R=simple_matrix(H)*A*simple_matrix(G).inverse()
    assert all(x.denominator()==1 for x in R.list())
    def restrict(f):
        terms={}
        for p,c in f.dict().items():
            weight=tuple(map(int,R*vector(QQ,p)));terms[weight]=terms.get(weight,0)+c
        return h.torus_ring(terms)
    for i in range(g.rank):assert restrict(g.theta[i])==sum(R[j,i]*h.theta[j] for j in range(h.rank))
    # Reconstruct the restriction of the character itself from the actual
    # weight dictionaries, rather than assuming the branching formula.
    restricted={}
    for p,m in g.weight_multiplicities:
        weight=tuple(map(int,R*vector(QQ,p)));restricted[weight]=restricted.get(weight,0)+m
    expected_weights={tuple(p):copies*m for p,m in h.weight_multiplicities}
    zero=(0,)*h.rank
    if trivials:expected_weights[zero]=expected_weights.get(zero,0)+trivials
    assert restricted==expected_weights
    coordinates=h.linear_basis_matrix.inverse()*R*g.linear_basis_matrix
    assert all(x.denominator()==1 for x in coordinates.list())
    factors=[int(x) for x in coordinates.change_ring(ZZ).elementary_divisors()]
    assert factors==expected
    return dict(group=G,fixed_group=H,root_node_map=root_nodes,
                representation_highest_weights=[list(w) for w in high_G],
                restricted_defining_copies=copies,restricted_trivial_copies=trivials,
                restriction_weight_matrix=[list(map(int,row)) for row in R],
                all_torus_Springer_coordinates_agree=True,
                quotient_of_full_fixed_source_by_restricted_source=factors)


def main():
    cases=[run('A3','C2',[1,2,1],((1,0,0),(0,0,1)),(1,0),2,0,[1,2]),
           run('E6','F4',[4,1,3,2,3,4],((1,0,0,0,0,0),(0,0,0,0,0,1)),(0,0,0,1),2,2,[1,1,2,2]),
           run('D4','G2',[1,2,1,1],((1,0,0,0),(0,0,1,0),(0,0,0,1)),(1,0),3,3,[1,3])]
    report=dict(result='passed',date='2026-10-22',cases=cases,
                scope='Exact finite root restriction, actual branching-weight dictionaries, Laurent Springer coordinates and integral source quotients; the uniform fixed-subgroup compatibility is proved separately, not inferred from these three examples.',
                source_sha256={'experiments/other_theories_folding_pilot.py':sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
    (ROOT/'experiments/other-theories-folding-pilot.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASSED three exact folding/source comparisons')


if __name__=='__main__':main()
