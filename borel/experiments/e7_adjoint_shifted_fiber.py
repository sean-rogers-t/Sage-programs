"""Direct exact Groebner fiber after a polynomial square-completion change."""
from hashlib import sha256
from pathlib import Path
from time import perf_counter
import json
import os
import re
import sys
from sage.all import PolynomialRing
from exceptional_degree_fibers import prepare_fiber,msolve_leading


def square_completed(R,equations,index=2):
    B=PolynomialRing(R.base_ring(),R.ngens()-1,names=['y%d'%(i+1) for i in range(R.ngens()-1)])
    U=PolynomialRing(B,'v')
    v=U.gen()
    images=list(B.gens())
    images.insert(index,v)
    source=[R.hom(images,U)(f) for f in equations]
    lc=source[0].leading_coefficient()
    assert source[0].degree()==2 and lc.is_constant() and lc
    Q=U(source[0]/lc)
    shift=Q[1]/2
    d=(Q[1]**2-4*Q[0])/4
    q=v*v-d
    assert Q(v-shift)==q
    completed=[q]
    for f in source[1:]:
        moved=f(v-shift)
        quotient,remainder=moved.quo_rem(q)
        assert moved==quotient*q+remainder and remainder.degree()<=1
        completed.append(remainder)
    from_base=B.hom([R.gen(i) for i in range(R.ngens()) if i!=index],R)
    def lift(f):
        return sum(from_base(c)*R.gen(index)**i for i,c in enumerate(f.list()))
    answer=[lift(f) for f in completed]
    return answer,dict(variable_index=index,polynomial_shift=str(shift),
        scalar_normalization=str(lc),exact_remainder_identities=len(source),
        original_total_degrees=[int(f.total_degree()) for f in equations],
        transformed_total_degrees=[int(f.total_degree()) for f in answer])


def main():
    prime=int(sys.argv[1]) if len(sys.argv)>1 else 65521
    seed=int(sys.argv[2]) if len(sys.argv)>2 else 20261021
    R,equations,targets,eliminated,chart=prepare_fiber('E7','adjoint',prime,seed,mode='eliminate-center')
    equations,change=square_completed(R,equations)
    print('Completed-square degrees:',change['transformed_total_degrees'],flush=True)
    start=perf_counter()
    leading,digest,version=msolve_leading(R,equations,prime,os.environ['MSOLVE_BINARY'])
    ideal=R.ideal([R({e:1}) for e in leading])
    dimension=int(ideal.dimension())
    assert dimension==0
    length=int(ideal.vector_space_dimension())
    report=dict(status='specialized-fiber-only',cartan_type='E7',representation='adjoint',
        characteristic=prime,seed=seed,targets=list(map(str,targets)),
        fiber_dimension=dimension,chart_fiber_length=length,fiber_length=chart['cover_degree']*length,
        engine=version+'; exact linear algebra option '+os.environ.get('MSOLVE_LINEAR_ALGEBRA','2'),
        solver_threads=int(os.environ.get('MSOLVE_THREADS','2')),
        solver_max_pairs=int(os.environ.get('MSOLVE_MAX_PAIRS','0')),
        eliminated_coordinate=eliminated,center_chart=chart,square_completion=change,
        monomial_order=str(R.term_order()),leading_exponents=leading,basis_sha256=digest,
        seconds=perf_counter()-start)
    path=Path(__file__).with_name('e7-adjoint-fiber-msolve-shifted-%d-%d.json'%(prime,seed))
    encoded=json.dumps(report,indent=2)
    encoded=re.sub(r'\[\n\s+(-?\d+(?:,\n\s+-?\d+)*)\n\s+\]',
                   lambda m:'['+','.join(re.findall(r'-?\d+',m.group(1)))+']',encoded)
    path.write_text(encoded+'\n')
    print({k:v for k,v in report.items() if k!='leading_exponents'},flush=True)


if __name__=='__main__':
    main()
