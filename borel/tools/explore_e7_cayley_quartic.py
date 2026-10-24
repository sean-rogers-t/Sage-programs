"""Reproduce the companion dimension calculation; not an independent certificate.

The rank-one dimension theorem in E7_SYMMETRIC_PENCIL.md does not depend on
this experiment. To promote its companion, verify reverse ideal inclusion
for the selected members, e.g. by supplying explicit multiplier identities.
"""
from hashlib import sha256
from pathlib import Path
import gzip
import json
from sage.all import QQ, GF, PolynomialRing, matrix
from sage.interfaces.singular import singular

ROOT=Path(__file__).resolve().parents[1]


def main():
    P=PolynomialRing(QQ,names='x,y,w,A,B,C',order='degrevlex')
    x,y,w,A,B,C=P.gens()
    h=A*A+B*B+C*C-A*B*C-4
    psi=x*x+y*y+w*w-x*x*y*y-x*x*w*w-y*y*w*w-x*y*w*(A*x+B*y+C*w)-A*B*x*y-A*C*x*w-B*C*y*w-A*B*C-3
    equations=[psi.derivative(v) for v in (x,y,w)]+[h]
    jacobian=matrix(P,[[f.derivative(v) for v in P.gens()] for f in equations])
    generators=equations+jacobian.minors(4)
    singular.lib('modstd.lib')
    # The documented exactness option is used, but the repository's stronger
    # independent inclusion audit is still pending. Do not promote this flag.
    basis_singular=singular.modStd(singular(P.ideal(generators)),1)
    basis=list(basis_singular.sage().gens())
    dimension=int(singular.dim(basis_singular))
    assert dimension==1
    modular=P.change_ring(GF(101))
    modular_dimension=modular.ideal([modular(f) for f in generators]).dimension()
    assert modular_dimension==1
    support=lambda f:frozenset(i for i,e in enumerate(f.exponents()[0]) if e)
    supports={support(f) for f in basis}
    minimal={s for s in supports if not any(t<s for t in supports)}
    chosen=[min((f for f in basis if support(f)==s),key=lambda f:(len(str(f)),f.total_degree()))
            for s in sorted(minimal,key=lambda s:tuple(s))]
    assert minimal=={frozenset(s) for s in ((0,),(1,),(2,4),(2,5),(3,),(4,5))}
    encode=lambda f:[[list(map(int,e)),str(c)] for e,c in f.dict().items()]
    archive=dict(variables=list(P.variable_names()),term_order='degrevlex',
                 original_generators=[encode(f) for f in generators],
                 computed_basis=[encode(f) for f in basis],
                 selected_members=[encode(f) for f in chosen],
                 reverse_ideal_inclusion_independently_certified=False)
    payload=gzip.compress((json.dumps(archive,separators=(',',':'))+'\n').encode(),mtime=0)
    name='e7-quartic-boundary-candidate.json.gz';(ROOT/name).write_bytes(payload)
    report=dict(status='exploratory',computed_singular_locus_dimension=dimension,
                basis_size=len(basis),selected_leading_monomials=[str(f.lm()) for f in chosen],
                source_algorithm='Singular modStd, exactness=1, over QQ',
                prime101_dimension_corroboration=int(modular_dimension),
                independent_reverse_ideal_inclusion_complete=False,
                companion_generic_smoothness_promoted=False,
                archive=name,archive_sha256=sha256(payload).hexdigest(),
                next_gate='Verify six polynomial memberships in the original Jacobian-minor ideal; their leading supports then prove dimension at most one without another Groebner computation.',
                source_sha256={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for n in ('E7_SYMMETRIC_PENCIL.md','tools/explore_e7_cayley_quartic.py')})
    (ROOT/'e7-quartic-boundary-exploration.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
