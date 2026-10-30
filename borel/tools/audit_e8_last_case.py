"""Native/SymPy audit of the last E8 geometry and degree certificate.

No Sage import. F4 membership replay is explicitly supplied by the two
separate verification reports; it is not inferred from a monomial count.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations,product
import gzip,json
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1]
def digest(path):return sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def key(e):return (sum(e)+sum(e[:4]),)+tuple(-v for v in reversed(e))

def degree_from_local_monomials(exponents):
    supports=[{i for i,v in enumerate(e) if v} for e in exponents]
    allowed=[v for n in range(9) for v in combinations(range(8),n) if not any(q<=set(v) for q in supports)]
    dimension=max(map(len,allowed));assert dimension==4
    records=[]
    for free in [v for v in allowed if len(v)==4]:
        normal=[i for i in range(8) if i not in free]
        projected={tuple(e[i] for i in normal) for e in exponents};minimal=[]
        for e in sorted(projected,key=lambda e:(sum(e),e)):
            if not any(all(a<=b for a,b in zip(q,e)) for q in minimal):minimal.append(e)
        bounds=[min(e[i] for e in minimal if e[i] and sum(bool(v) for v in e)==1) for i in range(4)]
        standard=[m for m in product(*(range(b) for b in bounds)) if not any(all(a<=b for a,b in zip(q,m)) for q in minimal)]
        records.append(dict(free=list(free),normal=normal,length=len(standard),standard=[list(m) for m in standard]))
    assert sum(v['length'] for v in records)==21
    return records

def main():
    # Bind to the already completed independent full-root and actual unit-
    # elimination audit; the character data have not changed in this phase.
    previous=json.loads((ROOT/'e8-final-cases-audit.json').read_text())
    assert previous['result']=='passed' and previous['native_full_root_characters']==5
    for name,want in previous['source_sha256'].items():assert digest(ROOT/name)==want,name
    row=next(v for v in json.loads((ROOT/'e8-final-character-models.json').read_text()) if v['seed']==[1,2,3,5,7])
    x,y,w,a,A,B,C,U,V=s.symbols('x y w a A B C U V')
    p=x*y+w*a;q=x*w+y*a;r=x*a+y*w
    f=x*x+y*y+w*w+a*a-p*p-q*q-r*r+x*y*w*a-A*q*r-B*p*r-C*p*q-B*C*p-A*C*q-A*B*r-A*B*C-4
    traces={A:U+1/U,B:V+1/V,C:U*V+1/(U*V)}
    actual=0
    for term in row['residual_terms']:
        i,j=term['charge'];assert i%3==0
        actual+=s.Rational(term['coefficient'])*s.prod(v**e for v,e in zip((x,y,w,a),term['exponent']))*U**(-i//3-2*j)*V**(i//3+j)
    assert s.expand(f.subs(traces,simultaneous=True)-actual)==0
    assert s.expand((A*A+B*B+C*C-A*B*C-4).subs(traces,simultaneous=True))==0
    for uu,vv in product((-1,1),repeat=2):
        assert all(s.simplify((z*s.diff(t,z)).subs({U:uu,V:vv}))==0 for t in traces.values() for z in (U,V))
    # Exact height input for the three irreducible degree-six components.
    G=s.groebner([p,q,r],x,y,w,a,order='grevlex',domain=s.QQ)
    expected=[w**3*a-w*a**3,x*x*a-w*w*a,y*y*a-w*w*a,p,q,r]
    assert {s.expand(v.as_expr()) for v in G.polys}==set(expected)
    supports=[{i for i,v in enumerate(z.LM(order='grevlex').exponents) if v} for z in G.polys]
    assert max(len(v) for n in range(5) for v in combinations(range(4),n) if not any(t<=set(v) for t in supports))==1
    # Independent sparse polynomial arithmetic for the exact normal jet.
    R,xx,yy,ww,aa=ring('x,y,w,a',s.QQ)
    pp=xx*yy+ww*aa;qq=xx*ww+yy*aa;rr=xx*aa+yy*ww
    M=[[ -z.diff(v) for z in (qq*rr,pp*rr,pp*qq)] for v in (xx,yy,ww)]
    det=lambda t:t[0][0]*(t[1][1]*t[2][2]-t[1][2]*t[2][1])-t[0][1]*(t[1][0]*t[2][2]-t[1][2]*t[2][0])+t[0][2]*(t[1][0]*t[2][1]-t[1][1]*t[2][0])
    D=det(M);assert D==-2*pp*qq*rr*(-aa**3+aa*(xx*xx+yy*yy+ww*ww)-2*xx*yy*ww)
    adj=[]
    for i in range(3):
        row_adj=[]
        for j in range(3):
            rows=[k for k in range(3) if k!=j];cols=[k for k in range(3) if k!=i]
            row_adj.append((-1)**(i+j)*(M[rows[0]][cols[0]]*M[rows[1]][cols[1]]-M[rows[0]][cols[1]]*M[rows[1]][cols[0]]))
        adj.append(row_adj)
    L=[-(pp*pp+qq*qq+rr*rr-xx*yy*ww*aa).diff(v) for v in (xx,yy,ww)]
    ell=[sum(adj[i][j]*L[j] for j in range(3)) for i in range(3)]
    c=D*sum(v*v for v in ell)+ell[0]*ell[1]*ell[2]-4*D**3
    lead=max(c,key=lambda e:(sum(e),)+tuple(-v for v in reversed(e)))
    assert lead==(10,10,6,1);coefficient=c[lead]
    # Universal divided-difference identity, independent of the large R.
    d,la,lb,lc,ta,va,vb,vc=s.symbols('d la lb lc ta va vb vc')
    sa,sb,sc=-ta*la,-ta*lb,-ta*lc
    H=ta*(va*va+vb*vb+vc*vc)-va*vb*vc-4*ta**3
    qa=d*d*ta*va+d*ta*sa-d*d*vb*vc
    qb=d*d*ta*vb+d*ta*sb-d*sa*vc
    qc=d*d*ta*vc+d*ta*sc-sa*sb
    assert s.expand(d**3*H-qa*(d*va-sa)-qb*(d*vb-sb)-qc*(d*vc-sc)-ta**3*(d*(la*la+lb*lb+lc*lc)+la*lb*lc-4*d**3))==0
    cases=[];bound={}
    for prime in (32003,32009):
        assert int(coefficient)%prime
        name=f'e8-last-{prime}-certificate.json.gz';raw=(ROOT/name).read_bytes();cert=json.loads(gzip.decompress(raw))
        orders=[]
        for z in cert['original_affine_generators']:
            deg=max(sum(e) for e,c in z);orders.append(deg-max(sum(e[:4]) for e,c in z))
        assert sorted(orders)[:4]==[1,1,1,3] and sum(v==1 for v in orders)==3
        exponents=[]
        for record in cert['selected_members']:
            assert len({sum(e) for e,co in record['polynomial']})==1
            assert all(int(co)%prime for e,co in record['polynomial'])
            lm=max((tuple(e) for e,co in record['polynomial']),key=key)
            assert lm==tuple(record['leading_exponent']);exponents.append(lm)
        exponents.append(lead+(0,0,0,3))
        count=degree_from_local_monomials(exponents)
        replay=json.loads((ROOT/f'e8-last-{prime}-verification.json').read_text())
        assert replay['result']=='passed' and replay['original_ideal_memberships_replayed']
        assert replay['raw_sha256'][name]==sha256(raw).hexdigest()
        cases.append(dict(prime=prime,selected_members=len(exponents)-1,independent_dimension=4,independent_degree=21,coordinate_lengths=count))
        bound[name]=sha256(raw).hexdigest()
    result=dict(result='passed',independent_compact_residual_from_audited_root_data=True,
        exact_normal_matrix_identity=True,exact_divided_difference_identity=True,
        independent_relation_leading_exponent=list(lead)+( [0,0,0,3]),
        relation_leading_coefficient=str(coefficient),normal_order_lower_bound=3,
        three_sextic_components_height_checked=True,known_infinity_degree_lower_bound=21,
        cases=cases,characteristic_zero_affine_singular_dimension_bound=2,
        covered_proper_E8_standard_subsets=255,remaining_Weyl_classes=0,
        scope='Independent native residual/jet/degree geometry; F4 ideal memberships separately replayed at two primes. No radical or generic parameter replacement.',
        nonemptiness_claim=False,literature_priority_claim=False)
    names=['E8_LAST_ADJOINT_FIBER.md','tools/audit_e8_last_case.py','e8-final-character-models.json',
           'e8-last-32003-verification.json','e8-last-32009-verification.json','e8-final-cases-audit.json']
    result['source_sha256']={n:digest(ROOT/n) for n in names};result['raw_sha256']=bound
    (ROOT/'e8-last-case-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('cases','source_sha256','raw_sha256')},indent=2))
if __name__=='__main__':main()
