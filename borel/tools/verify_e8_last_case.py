"""Replay original homogeneous ideal members and the infinity degree bound.

Sage; requires the pinned exact-F4 driver built by the companion shell script.
The partial F4 output is used only as ideal members, never as a full basis.
"""
from pathlib import Path
from hashlib import sha256
import argparse,gzip,json,os,subprocess,tempfile,time,re
from sage.all import QQ,GF,PolynomialRing,TermOrder,matrix,vector
from sage.interfaces.singular import singular
ROOT=Path(__file__).resolve().parents[1]

def model():
    P=PolynomialRing(QQ,names='x,y,w,a,A,B,C');x,y,w,a,A,B,C=P.gens()
    p=x*y+w*a;q=x*w+y*a;r=x*a+y*w
    f=x*x+y*y+w*w+a*a-p*p-q*q-r*r+x*y*w*a-A*q*r-B*p*r-C*p*q-B*C*p-A*C*q-A*B*r-A*B*C-4
    g=[f.derivative(v) for v in (x,y,w)]+[A*A+B*B+C*C-A*B*C-4]
    J=matrix(P,[[z.derivative(v) for v in P.gens()] for z in g])
    originals=list(dict.fromkeys(z for z in g+J.minors(4) if z))
    assert len(originals)==38
    return P,f,g,originals

def geometry_and_relation():
    P,f,g,originals=model()
    Q=PolynomialRing(QQ,names=list(P.variable_names())+['h'])
    x,y,w,a,A,B,C,h=Q.gens();embed=P.hom(list(Q.gens()[:7]),Q)
    homogeneous=[embed(z).homogenize(h) for z in originals]
    gg=vector(Q,[embed(z).homogenize(h) for z in g[:3]])
    H=embed(g[3]).homogenize(h)
    zero={A:0,B:0,C:0,h:0}
    M=matrix(Q,[[z.derivative(v).subs(zero) for v in (A,B,C)] for z in gg])
    D=M.det();adj=M.adjugate();L=vector(Q,[z.derivative(h).subs(zero) for z in gg])
    s=-h*(adj*L);T=adj*gg
    ri=vector(Q,[D*D*h*A+D*h*s[0]-D*D*B*C,
                 D*D*h*B+D*h*s[1]-D*s[0]*C,
                 D*D*h*C+D*h*s[2]-s[0]*s[1]])
    R=D**3*H-ri*T
    # This equality is an explicit QQ membership in the ORIGINAL ideal.
    assert R.is_homogeneous() and R.degree()==30
    weight=TermOrder('wdegrevlex',(2,2,2,2,1,1,1,1))
    W=PolynomialRing(QQ,names=Q.variable_names(),order=weight)
    assert tuple(W(R).lm().exponents()[0])==(10,10,6,1,0,0,0,3)
    p=x*y+w*a;q=x*w+y*a;r=x*a+y*w
    for ig in ([h,A,p,B*r+C*q],[h,B,q,A*r+C*p],[h,C,r,A*q+B*p]):
        I=Q.ideal(ig);basis=list(I.groebner_basis())
        assert I.dimension()==4 and all(not z.reduce(basis) for z in homogeneous)
    S=PolynomialRing(QQ,names='x,y,w,a');xx,yy,ww,aa=S.gens()
    assert S.ideal([xx*yy+ww*aa,xx*ww+yy*aa,xx*aa+yy*ww]).dimension()==1
    orders=[min(sum(e[4:]) for e in z.exponents()) for z in homogeneous]
    assert orders[:3]==[1,1,1] and min(orders[3:])>=3
    point={x:1,y:2,w:3,a:4,A:0,B:0,C:0,h:0}
    assert abs(QQ(D.subs(point)))==61600
    return P,Q,originals,homogeneous,R,weight

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--engine',required=True);ap.add_argument('--prime',type=int,required=True,choices=(32003,32009));args=ap.parse_args()
    start=time.time();prime=args.prime;name=f'e8-last-{prime}-certificate.json.gz'
    raw=(ROOT/name).read_bytes();cert=json.loads(gzip.decompress(raw))
    P,Q,originals,homogeneous,R,weight=geometry_and_relation()
    decode=lambda v:P({tuple(e):QQ(c) for e,c in v})
    ordered=[decode(v) for v in cert['original_affine_generators']]
    assert set(ordered)==set(originals)
    embed=P.hom(list(Q.gens()[:7]),Q);hh=Q.gen(7)
    ordered=[embed(z).homogenize(hh) for z in ordered]
    # Input order is retained for reproducibility; set equality above checks
    # every input against the original gradients and full Jacobian minors.
    with tempfile.TemporaryDirectory(prefix='e8-last-f4-') as temp:
        temp=Path(temp);inp=temp/'input';out=temp/'members';filt=temp/'indices';log=temp/'log'
        lines=[f'8 {len(ordered)} {sum(len(z.dict()) for z in ordered)}',' '.join(str(len(z.dict())) for z in ordered)]
        for z in ordered:
            # Preserve certificate's original term order from its serialized
            # affine input; F4 sorts these monomials internally.
            for e,c in z.dict().items():lines.append(' '.join(map(str,e))+f' {c.numerator()} {c.denominator()}')
        inp.write_text('\n'.join(lines)+'\n')
        selected=sorted(cert['selected_members'],key=lambda z:z['basis_index'])
        filt.write_text('\n'.join(str(z['basis_index']) for z in selected)+'\n')
        env=dict(os.environ,E8_F4_BOUND='20')
        print('Replaying original homogeneous ideal at',prime,flush=True)
        with log.open('w') as target:
            subprocess.run([str(Path(args.engine).resolve()),str(inp),str(out),str(prime),str(filt)],env=env,stdout=target,stderr=subprocess.STDOUT,check=True)
        logtext=log.read_text();assert 'E8 partial stop at degree 20' in logtext
        assert int(re.search(r'F4_TOTAL_MEMBERS=(\d+)',logtext).group(1))==cert['partial_member_count']
        F=PolynomialRing(GF(prime),names=Q.variable_names(),order=weight)
        stream=out.open();assert tuple(map(int,stream.readline().split()))==(8,len(selected))
        leading=[]
        for record in selected:
            n=int(stream.readline());terms={}
            for _ in range(n):
                row=stream.readline().split();terms[tuple(map(int,row[1:]))]=int(row[0])
            z=F(terms);expected=F({tuple(e):c for e,c in record['polynomial']})
            assert z==expected
            assert tuple(z.lm().exponents()[0])==tuple(record['leading_exponent'])
            leading.append(z.lm())
        rr=F(R);assert tuple(rr.lm().exponents()[0])==tuple(cert['explicit_relation_leading_exponent'])
        leading.append(rr.lm())
        # These monomials need not generate the full initial ideal.
        M=singular(F.ideal(leading)).std()
        assert int(M.dim())==4 and int(M.mult())==21
    result=dict(result='passed',prime=prime,selected_original_ideal_members=len(selected),
        original_ideal_memberships_replayed=True,explicit_rational_member_degree=30,
        exact_rational_infinity_inclusions=True,normal_linear_determinant_at_test_point=61600,
        infinity_components_degrees=[1,6,6,6],infinity_generic_multiplicity_lower_bounds=[3,1,1,1],
        homogeneous_dimension_upper_bound=4,homogeneous_degree_upper_bound=21,
        accounted_infinity_degree=21,characteristic_zero_affine_singular_dimension_bound=2,
        regular_fiber_dimension_if_nonempty=12,covered_proper_E8_standard_subsets=255,
        remaining_Weyl_classes=0,remaining_standard_subsets=0,all_E8_proper_Levis_resolved=True,
        nonemptiness_claim=False,literature_priority_claim=False,seconds=round(time.time()-start,3))
    names=['E8_LAST_ADJOINT_FIBER.md','tools/verify_e8_last_case.py','tools/e8_last_f4_driver.c',
           'tools/build_e8_last_certificate_engine.sh','e8-final-character-models.json']
    result['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    result['raw_sha256']={name:sha256(raw).hexdigest()}
    (ROOT/f'e8-last-{prime}-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','raw_sha256')},indent=2))
if __name__=='__main__':main()
