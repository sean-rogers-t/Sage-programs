"""Exact fiber length using a quadratic algebra as a free rank-two module.

No projection, saturation or radical is taken: nilpotent multiplicities are
retained. Run with Sage; optional arguments are prime, seed, dp/weighted.
"""
from hashlib import sha256
from pathlib import Path
from time import perf_counter
import json
import os
import re
import subprocess
import sys
import tempfile

from sage.all import GF, PolynomialRing
from exceptional_degree_fibers import prepare_fiber


def quadratic_module(R,equations,index):
    """Return B and relations for B[v]/(v^2-d, r_i), free over B before r_i."""
    B=PolynomialRing(R.base_ring(),R.ngens()-1,names=['x%d'%(i+1) for i in range(R.ngens()-1)])
    U=PolynomialRing(B,'v')
    v=U.gen()
    images=list(B.gens())
    images.insert(index,v)
    hom=R.hom(images,U)
    source=[hom(f) for f in equations]
    lc=source[0].leading_coefficient()
    assert source[0].degree()==2 and lc.is_constant() and lc
    Q=source[0]/lc
    shift=Q[1]/2
    d=(Q[1]**2-4*Q[0])/4
    q=v*v-d
    assert Q(v-shift)==q
    columns=[]
    identities=1
    for f in source[1:]:
        moved=f(v-shift)
        quotient,r=moved.quo_rem(q)
        assert moved==quotient*q+r and r.degree()<=1
        a,b=r[1],r[0]
        assert (r*v).quo_rem(q)[1]==d*a+b*v
        columns.extend([(b,a),(d*a,b)])
        identities+=2
    return B,columns,dict(quadratic_variable=index,quadratic_rank=2,
        square_completion_shift=str(shift),square_rhs_terms=len(d.dict()),
        square_rhs_degree=int(d.total_degree()),exact_reduction_identities=identities)


def module_length(B,columns,*,weights=None,label='module'):
    rank=len(columns[0])
    assert all(len(col)==rank for col in columns)
    root=Path(os.environ.get('MSOLVE_LOG_DIR',tempfile.gettempdir()))
    root.mkdir(parents=True,exist_ok=True)
    work=Path(tempfile.mkdtemp(prefix=label+'-',dir=str(root)))
    source=work/'module.sing'
    log=work/'progress.log'
    order='(dp,C)' if weights is None else '(am('+','.join(map(str,weights))+'),dp,C)'
    text='ring r=%d,(%s),%s;\n'%(B.base_ring().characteristic(),','.join(B.variable_names()),order)
    text+='option(prot);\nmodule rel='+','.join('['+','.join(map(str,col))+']' for col in columns)+';\n'
    text+='module gb=slimgb(rel);\nprint("");\n'
    text+='print("RESULT_DIM="+string(dim(gb)));\nprint("RESULT_LENGTH="+string(vdim(gb)));\n'
    text+='int j,k; vector v;\nfor(j=1;j<=size(gb);j++){v=lead(gb[j]);for(k=1;k<=%d;k++){if(v[k]!=0){print("LM="+string(k)+","+string(leadexp(v[k])));}}}\nquit;\n'%rank
    source.write_text(text)
    print('Singular module progress:',log,flush=True)
    start=perf_counter()
    with log.open('w') as stream:
        run=subprocess.run(['Singular','-q',str(source)],stdout=stream,stderr=subprocess.STDOUT)
    run.check_returncode()
    output=log.read_text()
    dim=int(re.search(r'^RESULT_DIM=(-?\d+)$',output,re.M)[1])
    length=int(re.search(r'^RESULT_LENGTH=(-?\d+)$',output,re.M)[1])
    leading=[]
    ideals=[[] for _ in range(rank)]
    for match in re.finditer(r'^LM=([0-9]+),([0-9, ]+)$',output,re.M):
        component=int(match[1])
        exponents=tuple(map(int,match[2].split(',')))
        assert len(exponents)==B.ngens()
        leading.append([component]+list(exponents))
        ideals[component-1].append(B({exponents:1}))
    assert leading and dim==0 and length>0
    component_lengths=[int(B.ideal(gens).vector_space_dimension()) for gens in ideals]
    assert sum(component_lengths)==length
    return dict(fiber_dimension=dim,chart_fiber_length=length,
        module_leading_terms=leading,module_component_lengths=component_lengths,
        base_rank=B.ngens(),module_rank=rank,monomial_order=order,
        engine=subprocess.check_output(['Singular','--version'],text=True).splitlines()[0]+'; slimgb module',
        input_sha256=sha256(source.read_bytes()).hexdigest(),
        basis_sha256=sha256(json.dumps(leading).encode()).hexdigest(),
        seconds=perf_counter()-start)


def preflight():
    B=PolynomialRing(GF(32003),2,names=('x','y'))
    x,y=B.gens()
    kill_y=[(y,B.zero()),(B.zero(),y)]
    # B[v]/(v^2-x,v-1) has length one.
    first=module_length(B,[(-B.one(),B.one()),(x,-B.one())]+kill_y,label='module-smoke')
    assert first['chart_fiber_length']==1
    # B[v]/(v^2-x,x^2) has length four, including its nilpotents.
    second=module_length(B,[(x*x,B.zero()),(B.zero(),x*x)]+kill_y,label='module-nilpotent')
    assert second['chart_fiber_length']==4


def main():
    prime=int(sys.argv[1]) if len(sys.argv)>1 else 32003
    seed=int(sys.argv[2]) if len(sys.argv)>2 else 20261020
    order=sys.argv[3] if len(sys.argv)>3 else 'weighted'
    assert order in ('dp','weighted')
    preflight()
    R,equations,targets,eliminated,chart=prepare_fiber('E7','adjoint',prime,seed,mode='eliminate-center')
    B,columns,reduction=quadratic_module(R,equations,2)
    weights=(2,1,2,2,2,0,5) if order=='weighted' else None
    report=module_length(B,columns,weights=weights,label='e7-adjoint-module')
    report.update(status='specialized-fiber-only',cartan_type='E7',representation='adjoint',
        characteristic=prime,seed=seed,targets=list(map(str,targets)),
        fiber_length=chart['cover_degree']*report['chart_fiber_length'],
        eliminated_coordinate=eliminated,center_chart=chart,quadratic_module_reduction=reduction,
        preflight_lengths=[1,4])
    path=Path(__file__).with_name('e7-adjoint-fiber-module-%d-%d-%s.json'%(prime,seed,order))
    encoded=json.dumps(report,indent=2)
    encoded=re.sub(r'\[\n\s+(-?\d+(?:,\n\s+-?\d+)*)\n\s+\]',
                   lambda m:'['+','.join(re.findall(r'-?\d+',m.group(1)))+']',encoded)
    path.write_text(encoded+'\n')
    print({k:v for k,v in report.items() if k!='module_leading_terms'},flush=True)


if __name__=='__main__':
    main()
