"""Exploratory exact fibers; specialized lengths are not generic degrees.

Run with Sage: kind E6/E7, representation minimal/adjoint, prime, seed.
The target point is recorded and exact finite-field Groebner arithmetic is
used. Lifting a specialized answer to characteristic zero needs a separate
flatness/standard-basis certificate and is not asserted by this script.
"""
from pathlib import Path
from time import perf_counter
import json
import random
import sys
import os
import re
import subprocess
import tempfile
from hashlib import sha256
from contextlib import nullcontext
from sage.all import GF, QQ, PolynomialRing, TermOrder
from sage_borel import e6_adjoint_quotient,e7_adjoint_quotient


def msolve_leading(R,equations,prime,binary):
    version=subprocess.check_output([binary,'-V'],text=True).strip()
    threads=int(os.environ.get('MSOLVE_THREADS','2'))
    linear=int(os.environ.get('MSOLVE_LINEAR_ALGEBRA','2'))
    pairs=int(os.environ.get('MSOLVE_MAX_PAIRS','0'))
    assert threads>0 and linear in (1,2), 'use positive threads and exact linear algebra 1 or 2'
    assert pairs>=0
    log_root=os.environ.get('MSOLVE_LOG_DIR')
    if log_root:
        Path(log_root).mkdir(parents=True,exist_ok=True)
        context=nullcontext(tempfile.mkdtemp(prefix='fiber-',dir=log_root))
    else:
        context=tempfile.TemporaryDirectory(prefix='borel-degree-fiber-')
    with context as tmp:
        source=Path(tmp)/'input.ms'
        output=Path(tmp)/'basis.ms'
        log=Path(tmp)/'progress.log'
        source.write_text(','.join(R.variable_names())+'\n'+str(prime)+'\n'+',\n'.join(map(str,equations))+'\n')
        command=[binary,'-f',str(source),'-o',str(output),'-g','1',
                 '-t',str(threads),'-l',str(linear),'-m',str(pairs),'-v','2']
        if log_root:
            print('msolve progress:',log,flush=True)
        with log.open('w') as stream:
            run=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
        (Path(tmp)/'run.json').write_text(json.dumps(dict(command=command,returncode=run.returncode,
            input_sha256=sha256(source.read_bytes()).hexdigest()),indent=2)+'\n')
        run.check_returncode()
        raw=output.read_bytes()
        body=''.join(line for line in raw.decode().splitlines() if not line.startswith('#'))
        body=body.strip('[]: \n')
        leading=[]
        for monomial in body.split(','):
            if not monomial:
                continue
            e=[0]*R.ngens()
            for i,a in re.findall(r'x([0-9]+)\^([0-9]+)',monomial):
                e[int(i)-1]=int(a)
            if not any(e):
                assert monomial.strip()=='1',monomial
            leading.append(tuple(e))
    return leading,sha256(raw).hexdigest(),version


def check_unit(R,equations,prime):
    binary=os.environ.get('MSOLVE_BINARY')
    if binary:
        leading,_,_=msolve_leading(R,equations,prime,binary)
        assert leading==[(0,)*R.ngens()]
    else:
        assert list(R.ideal(equations).groebner_basis(algorithm='libsingular:slimgb'))==[R.one()]


def prepare_fiber(kind,name,prime,seed,*,engine='drl',mode=None):
    """Exact quotient equations and optional globally checked chart changes."""
    if kind not in ('E6','E7') or name not in ('minimal','adjoint'):
        raise ValueError('choose E6/E7 and minimal/adjoint')
    field=GF(prime) if prime else QQ
    q=(e6_adjoint_quotient if kind=='E6' else e7_adjoint_quotient)(representation=name)
    marks=(1,2,2,3,2,1) if kind=='E6' else (2,2,3,4,3,2,1)
    R=PolynomialRing(field,len(marks),names=['x%d'%(i+1) for i in range(len(marks))],
                     order=(TermOrder('wdegrevlex',marks) if engine=='weighted'
                            else 'degrevlex'))
    rng=random.Random(seed)
    targets=[field(rng.randrange(1,1000)) for d in q.degrees]
    equations=[R({tuple(e):field(c) for e,c in f.dict().items()})-a
               for f,a in zip(q.character_polynomials,targets)]
    eliminated=None
    if mode and mode.startswith('eliminate'):
        j={('E6','minimal'):2,('E6','adjoint'):3,
           ('E7','minimal'):5,('E7','adjoint'):2}[kind,name]
        x=R.gen(j)
        c=equations[0].derivative(x)
        assert c.is_constant() and c and equations[0].degree(x)==1
        replacement=-(equations[0]-c*x)/field(c)
        S=PolynomialRing(field,len(marks)-1,names=['x%d'%(i+1) for i in range(len(marks)-1)])
        substitute=S({tuple(e[:j])+tuple(e[j+1:]):a for e,a in replacement.dict().items()})
        images=list(S.gens())
        images.insert(j,substitute)
        hom=R.hom(images,S)
        assert hom(equations[0])==0
        eliminated=dict(variable=str(x),constant_coefficient=str(c),replacement=str(replacement))
        equations=[hom(f) for f in equations[1:]]
        R=S
    center_chart=None
    if mode=='eliminate-center':
        assert (kind,name)==('E7','adjoint')
        # On chi7 != 0 use u=chi2/chi7, v=chi5/chi7, z=chi7^2.
        # Every source polynomial has even central parity. This is the
        # quotient by the free central involution on this chart.
        check_unit(R,equations+[R.gen(5)],prime)
        print('original chi7=0 divisor empty',flush=True)
        mapped=[]
        for f in equations:
            coefficients={}
            for e,c in f.dict().items():
                degree=e[1]+e[3]+e[5]
                assert degree%2==0
                exponent=tuple(e[:5])+(degree//2,)
                coefficients[exponent]=coefficients.get(exponent,field.zero())+c
            mapped.append(R(coefficients))
        equations=mapped
        check_unit(R,equations+[R.gen(5)],prime)
        print('quotient z=0 divisor empty',flush=True)
        center_chart=dict(coordinates='chi1,chi2/chi7,chi4,chi5/chi7,chi6,chi7^2',
                          excluded_original_divisor_empty=True,
                          excluded_quotient_divisor_empty=True,cover_degree=2)
    return R,equations,targets,eliminated,center_chart


def main():
    kind,name=sys.argv[1:3]
    prime=int(sys.argv[3]) if len(sys.argv)>3 else 32003
    seed=int(sys.argv[4]) if len(sys.argv)>4 else 20261020
    engine=sys.argv[5] if len(sys.argv)>5 else 'drl'
    mode=sys.argv[6] if len(sys.argv)>6 else None
    R,equations,targets,eliminated,center_chart=prepare_fiber(kind,name,prime,seed,engine=engine,mode=mode)
    I=R.ideal(equations)
    print(kind,name,'field',R.base_ring(),'targets',list(map(str,targets)),flush=True)
    start=perf_counter()
    if engine=='msolve':
        binary=os.environ['MSOLVE_BINARY']
        leading,basis_hash,version=msolve_leading(R,equations,prime,binary)
        J=R.ideal([R({e:1}) for e in leading])
        engine_name=version+'; exact linear algebra option '+os.environ.get('MSOLVE_LINEAR_ALGEBRA','2')
    else:
        G=I.groebner_basis(algorithm='libsingular:slimgb')
        print('basis',len(G),'seconds',perf_counter()-start,flush=True)
        leading=[tuple(f.lm().exponents()[0]) for f in G]
        basis_hash=sha256('\n'.join(map(str,G)).encode()).hexdigest()
        J=R.ideal([f.lm() for f in G])
        engine_name='Singular libsingular:slimgb'
    dimension=int(J.dimension())
    length=int(J.vector_space_dimension()) if dimension==0 else None
    chart_length=length
    if center_chart is not None and length is not None:
        length*=center_chart['cover_degree']
    print('fiber dimension',dimension,'length',length,'seconds',perf_counter()-start,flush=True)
    report=dict(status='specialized-fiber-only',cartan_type=kind,representation=name,
                characteristic=prime,seed=seed,targets=list(map(str,targets)),
                fiber_dimension=dimension,fiber_length=length,
                engine=engine_name,monomial_order=str(R.term_order()),
                eliminated_coordinate=eliminated,
                center_chart=center_chart,chart_fiber_length=chart_length,
                basis_sha256=basis_hash,leading_exponents=leading,
                seconds=perf_counter()-start)
    suffix='-center' if center_chart else ''
    path=Path(__file__).with_name('%s-%s-fiber-%s-%s-%s%s.json'%(kind.lower(),name,engine,prime,seed,suffix))
    encoded=json.dumps(report,indent=2)
    encoded=re.sub(r'\[\n\s+(-?\d+(?:,\n\s+-?\d+)*)\n\s+\]',
                   lambda m:'['+','.join(re.findall(r'-?\d+',m.group(1)))+']',encoded)
    path.write_text(encoded+'\n')


if __name__=='__main__':
    main()
