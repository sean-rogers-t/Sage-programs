"""Independent native root expansion and SymPy boundary audit (no Sage imports)."""
from pathlib import Path
from collections import Counter
from itertools import combinations
from hashlib import sha256
import json,time,gzip
import sympy as s
from audit_affine_branching_certificates import cartan,root_data

ROOT=Path(__file__).resolve().parents[1]
def conv(a,b):
    out=Counter()
    for x,c in a.items():
        for y,d in b.items():out[tuple(u+v for u,v in zip(x,y))]+=c*d
    return out
def decode(data,v):return s.Add(*(s.Rational(c)*s.prod(w**i for w,i in zip(v,e)) for e,c in data))
def monomial_dimension(polys,variables):
    supports=[set(i for i,e in enumerate(p.LM(order='grevlex').exponents) if e) for p in polys]
    for k in range(len(variables)+1):
        for chosen in combinations(range(len(variables)),k):
            if all(set(chosen)&support for support in supports):return len(variables)-k
    raise AssertionError('Monomial dimension failed')

def main():
    rows=json.loads((ROOT/'e8-final-character-models.json').read_text())
    report=json.loads((ROOT/'e8-final-cases-verification.json').read_text())
    A=cartan('E8');positive=root_data(A)[1];roots=list(positive)+[tuple(-v for v in r) for r in positive]
    for row in rows:
        columns=row['columns'];H=s.Matrix.hstack(*(s.Matrix(c) for c in columns));assert H.det()!=0
        k=row['derived']+row['central_semisimple']
        actual=Counter({tuple(sum(r[j]*columns[i][j] for j in range(8)) for i in range(8)):1 for r in roots})
        actual[(0,)*8]=8;fundamental={}
        for path in row['components']:
            N=len(path)+1;letters=[]
            for j in range(N):
                w=[0]*8
                for p,i in enumerate(path):w[i]=int(j==p)-int(j==p+1)
                letters.append(w)
            for j,i in enumerate(path):
                fundamental[i]=Counter(tuple(sum(letters[t][v] for t in subset) for v in range(8)) for subset in combinations(range(N),j+1))
        rebuilt=Counter()
        for q in row['terms']:
            weights={(0,)*k+tuple(q['charge']):q['coefficient']}
            for i,e in enumerate(q['exponent']):
                for _ in range(e):weights=conv(weights,fundamental[i])
            rebuilt.update(weights)
        assert {w:c for w,c in rebuilt.items() if c}==dict(actual)
        # Recompute the actual unit eliminations independently, in a flat
        # SymPy Laurent expression rather than a Sage coefficient tower.
        v=s.symbols(' '.join(row['variables']));phi=s.Add(*(s.Rational(q['coefficient'])*s.prod(w**e for w,e in zip(v,q['exponent']+q['charge'])) for q in row['terms']))
        pair={'1257':None,'12356':([0],[2]),'12357':([0],[2]),'12457':([1],[3]),'123567':([3,4],[0,5])}[''.join(map(str,row['seed']))]
        graph={};removed=[]
        if pair:
            taus=[v[i] for i in pair[0]];aa=[v[i] for i in pair[1]]
            fs=s.Matrix([s.diff(phi,w) for w in taus]);C=fs.jacobian(aa);det=s.expand(C.det())
            assert len(s.Add.make_args(det))==1 and not any(det.has(w) for w in v[:k])
            rhs=C*s.Matrix(aa)-fs;Ci=C.inv();sol=Ci*rhs
            graph=dict(zip(aa,map(s.expand,sol)));G=phi.subs(dict.fromkeys(taus,0))
            tau=-Ci.T*s.Matrix([s.diff(G,w).subs(graph,simultaneous=True) for w in aa])
            graph.update(zip(taus,map(s.expand,tau)));removed=pair[0]+pair[1]
            assert all(s.expand(s.diff(phi,w).subs(graph,simultaneous=True))==0 for w in aa+taus)
        rv=[v[i] for i in range(k) if i not in removed]+list(v[k:])
        expected=s.Add(*(s.Rational(q['coefficient'])*s.prod(w**e for w,e in zip(rv,q['exponent']+q['charge'])) for q in row['residual_terms']))
        assert s.expand(phi.subs(graph,simultaneous=True)-expected)==0
        for w in rv:
            assert s.expand(s.diff(phi,w).subs(graph,simultaneous=True)-s.diff(expected,w))==0
    print('Native Cartan/root expansion: all five full characters pass',flush=True)
    x,y,a,t=s.symbols('x y a t')
    aa=x*x-(t+1/t)*a*x+a*a-1
    bb=a*x**3-2*t*a*a*x*x+a**3*x+(t-1/t)*x*x-a*x+(t-1/t)*(a*a-1)
    cc=-t*a**3*x**3+t*a*x**3+a*a*x*x+t*a**3*x-(2*t+1/t)*a*x-3
    p=s.Poly(s.expand(t**5*(4*aa**2*s.diff(cc,x)-2*aa*bb*s.diff(bb,x)+bb**2*s.diff(aa,x))),x)
    assert s.gcd_list(p.all_coeffs())==t*t
    assert s.gcd(s.expand(t*aa),s.expand(t*bb))==1
    sp=s.Poly(p.as_expr().subs({a:3,t:2}),x);assert sp.degree()==7 and s.gcd(sp,sp.diff()).degree()==0
    x0=(t+1/t)*a/2
    f=s.Poly(s.cancel(aa.subs(x,x0)).as_numer_denom()[0],a)
    g=s.Poly(s.cancel(bb.subs(x,x0)).as_numer_denom()[0],a)
    resultant=s.factor(s.resultant(f,g,a))
    assert s.cancel(resultant/(256*(t-1)**2*(t+1)**2*t**4*(t*t+1)**4*(t**4+1)**2)) in (1,-1)
    print('Binary primitivity, squarefree specialization, and resultant pass',flush=True)
    pair=next(c for c in report['cases'] if c['seed']=='12356');checks=[]
    for cert in pair['boundary_certificates']:
        v=s.symbols(' '.join(cert['variables']));original=[decode(f,v) for f in cert['generators']]
        start=time.time();G=s.groebner(original,*v,order='grevlex',domain=s.QQ)
        assert all(G.reduce(decode(f,v))[1]==0 for f in cert['basis'])
        d=monomial_dimension(G.polys,v);assert d==cert['dimension']
        checks.append(dict(chart=cert['chart'],dimension=d,basis_size=len(G.polys)))
        print('Independent rational boundary basis',checks[-1],'seconds',round(time.time()-start,2),flush=True)
    multiplier=next(c for c in report['cases'] if c['seed']=='12457')['certificate']
    v=s.symbols(' '.join(multiplier['variables']));original=[decode(f,v) for f in multiplier['original_generators']]
    start=time.time();G=s.groebner(original,*v,order='grevlex',modulus=32009)
    d=monomial_dimension(G.polys,v);assert d==3
    print('Independent homogeneous multiplier bound passes',len(G.polys),'seconds',round(time.time()-start,2),flush=True)
    four=[]
    for prime in (32003,32009):
        name=f'e8-a1-four-{prime}-certificate.json.gz'
        data=json.loads(gzip.decompress((ROOT/name).read_bytes()))
        supports=[]
        for member in data['selected_members']:
            lead=max((tuple(e) for e,c in member['polynomial']),key=lambda e:(sum(e),tuple(-j for j in reversed(e))))
            assert lead==tuple(member['leading_exponent'])
            supports.append(set(i for i,e in enumerate(lead) if e))
        allowed=[q for k in range(10) for q in combinations(range(9),k) if all(not p<=set(q) for p in supports)]
        assert max(map(len,allowed))==4
        replay=json.loads((ROOT/f'e8-a1-four-{prime}-verification.json').read_text())
        assert replay['result']=='passed' and replay['original_ideal_memberships_replayed']
        assert replay['raw_sha256'][name]==sha256((ROOT/name).read_bytes()).hexdigest()
        four.append(dict(prime=prime,selected_members=len(supports),independent_monomial_dimension=4))
    result=dict(result='passed',native_full_root_characters=5,independent_actual_scheme_eliminations=5,
                binary_identities_over_Q=True,A1_four_independent_support_checks=four,
                independently_recomputed_rational_boundary_bases=checks,
                multiplier_homogeneous_dimension=3,multiplier_independent_prime=32009,
                geometry_audited='E8_FINAL_CASES.md; no reduction of the actual fiber schemes',
                all_E8_Levis_resolved=False,literature_priority_claim=False)
    names=['tools/audit_e8_final_cases.py','e8-final-character-models.json','e8-final-cases-verification.json','E8_FINAL_CASES.md',
           'e8-a1-four-32003-verification.json','e8-a1-four-32009-verification.json']
    result['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'e8-final-cases-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))
if __name__=='__main__':main()




