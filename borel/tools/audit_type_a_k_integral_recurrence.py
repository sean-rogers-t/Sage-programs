"""Independent integer replay of difference-trace annihilation witnesses.

No Sage, HNF, ideal solver, SymPy, or builder import is used for replay.
The all-rank theorem is a written proof, not a conclusion from this sample.
"""
from pathlib import Path
from hashlib import sha256
from math import comb, gcd
import gzip
import json

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    data = path.read_bytes()
    return sha256(data if path.suffix == '.gz' else data.replace(b'\r\n',b'\n')).hexdigest()


def add(out, poly, scalar=1):
    for e, value in poly.items():
        out[e] = out.get(e,0) + scalar*value
        if out[e] == 0:
            del out[e]


def multiply(left, right):
    result = {}
    for a,x in left.items():
        for b,y in right.items():
            e = tuple(i+j for i,j in zip(a,b))
            result[e] = result.get(e,0)+x*y
    return {e:v for e,v in result.items() if v}


def powers(poly, bound, one):
    result=[one]
    for _ in range(bound):
        result.append(multiply(result[-1],poly))
    return result


def evaluate(terms, values, one):
    maximum=[max((e[i] for e,v in terms),default=0) for i in range(len(values))]
    cache=[powers(value,n,one) for value,n in zip(values,maximum)]
    out={}
    for e,v in terms:
        assert type(v) is int and v != 0
        term=one
        for i,k in enumerate(e):
            assert type(k) is int and k>=0
            if k:
                term=multiply(term,cache[i][k])
        add(out,term,v)
    return out


def coefficient_model(m):
    zero=(0,)*m;one={zero:1}
    L={(1,)+(0,)*(m-1):1}
    lp=powers(L,2*m+4,one)
    c=[one]
    for k in range(1,m):
        e=list(zero);e[k]=1;c.append({tuple(e):1})
    cm={}
    for k in range(m):
        add(cm,multiply(c[k],lp[m-k]),(-1)**(m+1+k))
    c.append(cm)
    f=[]
    # Coefficients of F(L+t), in a polynomial ambient ring with c_m
    # eliminated by the actual monic root equation.
    for j in range(1,m):
        value={}
        for k in range(j+1):
            add(value,multiply(c[k],lp[j-k]),(-1)**k*comb(m-k,j-k))
        f.append(value)
    return c[1:]+f,one


def root_model(m):
    zero=(0,)*m;one={zero:1};roots=[]
    for i in range(m):
        e=list(zero);e[i]=1;roots.append({tuple(e):1})
    def elementary(values):
        out=[one]
        for x in values:
            out.append({})
            for j in range(len(out)-1,0,-1):
                add(out[j],multiply(x,out[j-1]))
        return out[1:]
    differences=[]
    for x in roots[1:]:
        d=dict(roots[0]);add(d,x,-1);differences.append(d)
    return elementary(roots)+elementary(differences),one


def main():
    path=ROOT/'experiments/type-a-k-integral-recurrence.json.gz'
    data=json.loads(gzip.decompress(path.read_bytes()))
    assert data['result']=='passed'
    for name,value in data['source_sha256'].items():
        assert digest(ROOT/name)==value
    cases=[];root_checks=0
    for case in data['cases']:
        m=case['m'];d=case['degree'];weights=list(range(1,m+1))+list(range(1,m))
        for key in ('relation','witness'):
            assert len({tuple(e) for e,v in case[key]})==len(case[key])
            assert all(len(e)==2*m-1 and sum(x*w for x,w in zip(e,weights))==d
                       for e,v in case[key])
        assert all(sum(e[:m])>0 for e,v in case['witness'])
        values,one=coefficient_model(m)
        residual=evaluate(case['relation'],values,one)
        residual={e:m*v for e,v in residual.items()}
        add(residual,evaluate(case['witness'],values,one),-1)
        assert not residual,(m,case['indices'])
        if m<=4:
            values,one=root_model(m)
            residual=evaluate(case['relation'],values,one)
            residual={e:m*v for e,v in residual.items()}
            add(residual,evaluate(case['witness'],values,one),-1)
            assert not residual
            root_checks+=1
        cases.append({k:case[k] for k in ('m','kind','indices','degree','term_count')})

    content_checks=0
    for m in range(2,201):
        content=[1]
        for d in range(1,3*m+1):
            value=0
            for j in range(1,min(m,d)+1):
                value=gcd(value,comb(m,j)*content[d-j])
            content.append(value)
            assert value==m//gcd(m,d)
            if d>=m:
                assert value==content[d-m]
            content_checks+=1
    translation_checks=0
    for m in range(2,41):
        for j in range(m+1):
            for k in range(j+1):
                for r in range(j-k+1):
                    assert comb(m-k,j-k)*comb(j-k,r)==comb(m-j+r,r)*comb(m-k,j-k-r)
                    translation_checks+=1

    names=['TYPE_A_K_INTEGRAL_RECURRENCE.md',
           'experiments/type_a_k_integral_recurrence.py',
           'experiments/type-a-k-integral-recurrence.json.gz',
           'tools/audit_type_a_k_integral_recurrence.py']
    report=dict(result='passed',scope='Finite exact implementation checks of the written all-rank trace recurrence.',
        universal_annihilator_status='written proof; not inferred from finite checks',
        source_image_preserved=True, integer_coefficient_ideal_witnesses=len(cases),
        independent_root_variable_replays=root_checks,
        content_recurrence_checks=content_checks,hasse_binomial_checks=translation_checks,
        cases=cases,source_sha256={name:digest(ROOT/name) for name in names})
    (ROOT/'type-a-k-integral-recurrence-audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('cases','source_sha256')},indent=2))


if __name__=='__main__':
    main()
