"""Direct BGG checks of the transfer formula on four different flag varieties.

Neither check obtains its image matrix from the Cayley-plane tensor model.
Full flag: degrees 0..4; Levi D4 (omitted nodes 1,6): degrees 0..6.
Retained A1+A1 and maximal A4+A1: degrees 0..4. The latter constructs
the full integer invariant kernels instead of assuming generators.
"""
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from time import perf_counter
import json
from sage.all import QQ, ZZ, PolynomialRing, matrix, prod
from sage_borel import BorelMap
from e6_root_weight_comparison import powers

ROOT=Path(__file__).resolve().parents[1]


def main():
    start=perf_counter()
    baseline=json.loads((ROOT/'experiments/e6-root-weight-comparison.json').read_text())
    output=[]
    for omitted,bound in ((None,4),((1,6),6),((2,3,4,5),4),((3,),4)):
        b=BorelMap('E6',omitted)
        w=b.weights
        if omitted is None:
            inputs=list(w); degrees=(1,)*6
            roots=[b.simple_roots[i] for i in range(1,7)]
            fiber=baseline['d5_refinements'][0]
        elif omitted==(1,6):
            ys=[-w[0]-2*w[1]+2*w[2],-w[0]+2*w[1]+2*w[2]-2*w[3],
                -w[0]+2*w[3]-2*w[4],-w[0]+2*w[4]-2*w[5],-w[0]+2*w[5]]
            def elementary(k):return sum(prod(a*a for a in subset) for subset in combinations(ys[:4],k))
            inputs=[w[0],ys[4],elementary(1),elementary(2),prod(ys[:4]),elementary(3)]
            roots=[3*inputs[0]]+inputs[1:]
            degrees=(1,1,2,4,4,6)
            fiber=next(row for row in baseline['d5_refinements'] if set(row['e6_retained'])=={2,3,4,5})
        elif omitted==(2,3,4,5):
            alpha=b.simple_roots
            inputs=[w[i-1] for i in (2,3,4,5)]+[alpha[1]**2,alpha[6]**2]
            roots=[alpha[2],alpha[3]+alpha[1]/2,alpha[4],alpha[5]+alpha[6]/2]+inputs[4:]
            degrees=(1,1,1,1,2,2)
            fiber=next(row for row in baseline['proper_parabolics'] if row['e6_retained']==[1,6])
        else:
            fiber=next(row for row in baseline['proper_parabolics'] if row['e6_retained']==[1,2,4,5,6])
        for d in range(bound+1):
            labels=b.basis(d)
            if omitted!=(3,):
                es=list(powers(degrees,d))
            def image(gs):
                columns=[]
                for e in es:
                    expanded=b.expand(prod(f**a for f,a in zip(gs,e)))
                    columns.append([expanded.get(v,0) for v in labels])
                return matrix(QQ,columns).transpose()
            if omitted!=(3,):
                N=image(inputs); M=image(roots)
                source_columns=len(es)
            else:
                # Full integral invariant kernels in root and weight bases;
                # no chosen invariant-generator presentation is assumed.
                R=PolynomialRing(QQ,6,'x'); xs=R.gens()
                es=list(powers((1,)*6,d))
                mon=[prod(x**a for x,a in zip(xs,e)) for e in es]
                C=matrix(ZZ,6,6,lambda i,j:b.simple_roots[i+1].monomial_coefficient(w[j]))
                def invariant_image(root_basis):
                    equations=[]
                    for i in b.retained_nodes:
                        j=i-1
                        if root_basis:
                            reflected=[xs[k]-C[k,j]*xs[j] for k in range(6)]
                        else:
                            alpha=sum(C[j,k]*xs[k] for k in range(6))
                            reflected=[xs[k]-(alpha if k==j else 0) for k in range(6)]
                        diffs=[f(*reflected)-f for f in mon]
                        equations+= [[f.monomial_coefficient(prod(x**a for x,a in zip(xs,e))) for f in diffs] for e in es]
                    K=matrix(ZZ,equations).right_kernel().basis_matrix()
                    originals=list(b.simple_roots.values()) if root_basis else list(w)
                    inputs0=[prod(f**a for f,a in zip(originals,e)) for e in es]
                    polys=[sum(c*f for c,f in zip(row,inputs0)) for row in K]
                    columns=[]
                    for f in polys:
                        expansion=b.expand(f)
                        columns.append([expansion.get(v,0) for v in labels])
                    return matrix(QQ,columns).transpose()
                N=invariant_image(False);M=invariant_image(True)
                source_columns=M.ncols()
            denominator=N.denominator()
            H=matrix(ZZ,denominator*N).transpose().hermite_form(include_zero_rows=False).transpose()/denominator
            assert H.is_square() and H.rank()==len(labels)
            transition=H.change_ring(QQ).inverse()*M
            assert transition.denominator()%3
            smith=matrix(ZZ,transition.denominator()*transition).elementary_divisors()
            exponents=[int(a.valuation(3)) for a in smith]
            nonunit={str(a):exponents.count(a) for a in (1,2,3)}
            expected={str(a):(fiber['multiplicities'][str(a)][d] if d<len(fiber['multiplicities'][str(a)]) else 0) for a in (1,2,3)}
            assert nonunit==expected,(omitted,d,nonunit,expected)
            assert max(exponents,default=0)<=3
            output.append(dict(omitted_nodes=list(b.omitted_nodes),degree=d,rank=len(labels),
                               source_columns=source_columns,multiplicities=nonunit))
            print('direct',b.omitted_nodes,d,nonunit,round(perf_counter()-start,2),flush=True)
    paths=[Path(__file__).resolve(),ROOT/'experiments/e6-root-weight-comparison.json',
           ROOT/'experiments/e6_root_weight_comparison.py']
    result=dict(result='passed',cases=output,seconds=round(perf_counter()-start,3),
                scope='22 direct BGG image comparisons on full flags, the D4 refinement, retained A1+A1 and maximal A4+A1. The last uses full integral invariant kernels; not all 60 parabolics are directly recomputed.',
                source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'experiments/e6-root-weight-refinements.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASSED',result['seconds'],flush=True)


if __name__=='__main__':main()
