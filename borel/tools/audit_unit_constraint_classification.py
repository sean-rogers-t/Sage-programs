"""Sage-free exhaustive census, determinants, and sparse chain-rule audit."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import gzip
import json
import sympy as s
from audit_affine_branching_certificates import cartan, root_data, inverse
from audit_branching_norm_fibers import reconstruct

ROOT=Path(__file__).resolve().parents[1]


def plus(*polys):
    result={}
    for p in polys:
        for e,c in p.items():result[e]=result.get(e,F(0))+c
    return {e:c for e,c in result.items() if c}


def scale(p,c):return {e:v*c for e,v in p.items() if v*c}


def times(p,q):
    result={}
    for a,c in p.items():
        for b,d in q.items():
            e=tuple(x+y for x,y in zip(a,b));result[e]=result.get(e,F(0))+c*d
    return {e:c for e,c in result.items() if c}


def derivative(p,i,log=False):
    result={}
    for e,c in p.items():
        if not e[i]:continue
        v=list(e)
        if not log:v[i]-=1
        result[tuple(v)]=c*e[i]
    return result


def at_zero(p,indices):return {e:c for e,c in p.items() if not any(e[i] for i in indices)}


def substitute(p,graph,length):
    result={};unit={(0,)*length:F(1)};powers={}
    for i,value in graph.items():
        powers[i,0]=unit
        for k in range(1,max((e[i] for e in p),default=0)+1):
            powers[i,k]=times(powers[i,k-1],value)
    for e,c in p.items():
        base=list(e)
        for i in graph:base[i]=0
        term={tuple(base):c}
        for i in graph:term=times(term,powers[i,e[i]])
        result=plus(result,term)
    return result


def from_expr(expr,z,n):
    answer={}
    for term in s.Add.make_args(s.expand(expr)):
        powers=term.as_powers_dict();charge=tuple(int(powers.get(v,0)) for v in z)
        coefficient=s.cancel(term/s.prod(v**e for v,e in zip(z,charge)))
        assert coefficient.is_Rational
        if coefficient:answer[(0,)*n+charge]=F(coefficient)
    return answer


def audit_case(record):
    row=record['character'];reconstruct(row)
    n=len(row['nodes']);d=len(row['central_cocharacters']);length=n+d
    phi={tuple(t['exponent']+t['charge']):F(t['coefficient']) for t in row['terms']}
    z=s.symbols('z0:'+str(d));partials=[derivative(phi,i) for i in range(n)]
    hessian=[[derivative(p,j) for j in range(n)] for p in partials]
    def constant_poly(p):return all(not any(e[:n]) for e in p)
    def expr(p):return sum(s.Rational(c.numerator,c.denominator)*s.prod(v**k for v,k in zip(z,e[n:])) for e,c in p.items())
    accepted=[]
    for k in range(1,n//2+1):
        for tau in combinations(range(n),k):
            if any(hessian[i][j] for i,j in product(tau,repeat=2)):continue
            for paired in combinations([i for i in range(n) if i not in tau],k):
                if any(not constant_poly(hessian[i][j]) for i in tau for j in paired):continue
                C=s.Matrix([[expr(hessian[i][j]) for j in paired] for i in tau])
                det=from_expr(s.expand(C.det()),z,n)
                if len(det)!=1:continue
                e,c=next(iter(det.items()));assert c.denominator==1
                accepted.append(dict(multipliers=list(tau),paired=list(paired),
                                     determinant_charge=list(e[n:]),determinant_coefficient=int(c)))
    maximum=max((len(r['paired']) for r in accepted),default=0)
    best=[r for r in accepted if len(r['paired'])==maximum]
    cert=record['constraints']
    assert len(accepted)==cert['accepted_systems'] and maximum==cert['maximum_pairs']
    assert best==cert['maximum_systems']
    chosen=min(best,key=lambda r:(r['multipliers'],r['paired'])) if best else None
    assert chosen==cert['chosen']
    if chosen:
        tau=chosen['multipliers'];paired=chosen['paired']
        C=s.Matrix([[expr(hessian[i][j]) for j in paired] for i in tau])
        Ci=[[from_expr(v,z,n) for v in row] for row in C.inv().tolist()]
        rhs=[scale(at_zero(partials[i],paired),-1) for i in tau]
        alpha={j:plus(*(times(Ci[k][l],rhs[l]) for l in range(maximum))) for k,j in enumerate(paired)}
        G=at_zero(phi,tau);Ga=[substitute(derivative(G,j),alpha,length) for j in paired]
        graph=dict(alpha)
        for k,i in enumerate(tau):graph[i]=scale(plus(*(times(Ci[l][k],Ga[l]) for l in range(maximum))),-1)
        reduced=substitute(G,alpha,length)
    else:
        assert all(f=='A1' for f in row['factors'])
        graph={0:scale(at_zero(partials[0],[0]),F(-1,2))}
        reduced=substitute(phi,graph,length)
    remaining=[i for i in range(n) if i not in graph]
    expected={}
    for term in record['residual']['terms']:
        e=[0]*n
        for i,v in zip(remaining,term['exponent']):e[i]=v
        expected[tuple(e+term['charge'])]=F(term['coefficient'])
    assert reduced==expected
    assert max(sum(e[:n]) for e in reduced)==record['residual']['degree']==4
    for i in range(length):
        original=substitute(derivative(phi,i,log=i>=n),graph,length)
        target={} if i in graph else derivative(reduced,i,log=i>=n)
        assert original==target,(row['type'],row['nodes'],i)
    return dict(type=row['type'],nodes=row['nodes'],maximum_pairs=maximum,
                remaining_variables=len(remaining),degree=4,
                every_original_derivative_checked=True)


def census():
    results={};cubic=set()
    for kind in ('E6','E7','E8'):
        A=cartan(kind);n=len(A);pos=root_data(A)[1]
        roots=list(pos)+[tuple(-x for x in r) for r in pos];counts=Counter()
        for mask in range(2**n-1):
            I=[i for i in range(n) if mask>>i&1];todo=set(I);parts=[]
            while todo:
                part={min(todo)}
                while True:
                    grown=part|{j for i in part for j in todo if i!=j and A[i][j]}
                    if grown==part:break
                    part=grown
                todo-=part;parts.append(part)
            if any(sum(i!=j and bool(A[i][j]) for j in part)>2 for part in parts for i in part):
                counts['non_A']+=1;continue
            counts['all_A']+=1
            degree=max((sum(r[i] for i in I) for r in roots if all(r[i]>=0 for i in I)),default=0)
            if degree==3:
                counts['cubic']+=1;cubic.add((kind,tuple(i+1 for i in I)))
            else:assert degree<=2;counts['quadratic']+=1
        results[kind]=dict(counts)
    return results,cubic


def main():
    source=json.loads((ROOT/'unit-constraint-classification-verification.json').read_text())
    archive=(ROOT/source['case_archive']).read_bytes()
    assert sha256(archive).hexdigest()==source['case_archive_sha256']
    rows=json.loads(gzip.decompress(archive));counts,expected=census()
    assert {(r['character']['type'],tuple(r['character']['nodes'])) for r in rows}==expected
    assert len(rows)==len(expected)==136
    checks=[]
    for row in rows:
        checks.append(audit_case(row))
        if len(checks)%25==0:print('audited',len(checks),'of',len(rows),flush=True)
    for kind in counts:assert counts[kind]==source['summary'][kind]['counts']
    assert all(r['remaining_variables']<=4 for r in checks)
    # Verify explicit words on complete subsystems, not just diagram types.
    A=cartan('E6');Ai=inverse(A);positive=root_data(A)[1]
    roots=set(positive)|{tuple(-v for v in r) for r in positive}
    def subsystem(nodes):
        outside=set(range(6))-{i-1 for i in nodes}
        return {r for r in roots if all(sum(Ai[i][j]*r[j] for j in range(6))==0 for i in outside)}
    targets=set()
    for witness in source['E6_Weyl_transfers']:
        assert witness['seed'] in ([2,3,5],[2,3,5,6],[1,2,3,5,6])
        current=subsystem(witness['seed'])
        for node in witness['word']:
            i=node-1;current={tuple(r[j]-r[i]*A[j][i] for j in range(6)) for r in current}
        assert current==subsystem(witness['target'])
        targets.add(tuple(witness['target']))
    assert targets=={nodes for kind,nodes in expected if kind=='E6'}
    report=dict(result='passed',uses_Sage=False,cases_audited=len(checks),
                all_445_proper_standard_subsets_enumerated=True,counts=counts,
                maximum_remaining_variables=4,all_residuals_quartic=True,
                independent_subset_character_reconstruction=True,
                independent_symbolic_determinants=True,all_original_target_derivatives=True,
                E6_explicit_Weyl_transfers_verified=len(targets),
                all_E6_proper_Levi_coverage_verified=True,
                recursive_reduction_classified=False,all_adjoint_generic_smoothness_proved=False,
                literature_priority_claim=False,cases=checks)
    names=['EXCEPTIONAL_QUARTIC_REDUCTION.md','tools/audit_unit_constraint_classification.py','tools/audit_affine_branching_certificates.py',
           'tools/audit_branching_norm_fibers.py','unit-constraint-classification-verification.json']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'unit-constraint-classification-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('cases','source_sha256')},indent=2))


if __name__=='__main__':main()
