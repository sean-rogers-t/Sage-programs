"""Sage-free root, exact-identity and normalization audit of the binary theorem."""
from collections import Counter
from hashlib import sha256
from pathlib import Path
import gzip
import json
import sympy as s
from audit_affine_branching_certificates import cartan,root_data
from audit_branching_norm_fibers import reconstruct

ROOT=Path(__file__).resolve().parents[1]


def polynomial(row,affine,central):
    reconstruct(row)  # Independent Cartan reflections and exterior subset weights.
    return s.expand(sum(s.Rational(t['coefficient'])*s.prod(v**e for v,e in zip(affine,t['exponent']))*
                        s.prod(v**e for v,e in zip(central,t['charge'])) for t in row['terms']))


def root_cover():
    A=s.Matrix(cartan('E7'));positive=root_data(cartan('E7'))[1]
    roots=list(positive)+[tuple(-v for v in r) for r in positive]
    e=s.eye(7)
    H=s.Matrix.hstack(e[:,0],e[:,1],e[:,4],e[:,6],s.Matrix([1,1,2,2,1,0,0]),
                      s.Matrix([1,1,2,2,2,2,1]),s.Matrix([1,2,2,4,3,2,1]))
    assert H.T*A*H==2*s.eye(7) and abs(H.det())==8
    vs=s.symbols('tX tx ty tw U V W')
    X,x,y,w,a,b,c=[v+1/v for v in vs]
    phi=X*X+x*x+y*y+w*w+a*X*x*y+b*X*x*w+c*X*y*w+b*c*x*y+a*c*x*w+a*b*y*w+a*b*c*X+a*a+b*b+c*c-7
    observed=Counter()
    for term in s.expand(phi).as_ordered_terms():
        coef,mono=term.as_coeff_Mul();powers=mono.as_powers_dict()
        observed[tuple(int(powers.get(v,0)) for v in vs)]+=int(coef)
    expected=Counter({tuple(int(v) for v in s.Matrix([r])*H):1 for r in roots})
    expected[(0,)*7]+=7
    assert observed==expected and sum(observed.values())==133
    def trace(cols):
        E=s.Matrix([[sum(r[i]*h[i] for i in range(7)) for h in cols] for r in roots])
        return E.T*E
    trace1=trace([H[:,j] for j in (4,5,6,3)])
    assert trace1==72*s.eye(4)
    hp=s.Matrix([4,6,8,12,9,6,3]);hq=s.Matrix([2,3,4,6,5,4,2])
    H2=s.Matrix.hstack(e[:,0],e[:,1],e[:,2],e[:,4],e[:,6],hp,hq)
    trace2=trace([hp,hq,e[:,6]])
    assert H2.det()!=0 and trace2.det()!=0
    assert trace2==s.Matrix([[432,216,0],[216,144,0],[0,0,72]])
    return dict(seven_A1_root_character=True,cover_indices=[int(abs(H.det())),int(abs(H2.det()))],
                central_trace_matrices=[list(map(list,trace1.tolist())),list(map(list,trace2.tolist()))])


def identities(source):
    X,x,y,w,A,B,C,p,q,r,U,V,W,t=s.symbols('X x y w A B C p q r U V W t')
    full=X*X+x*x+y*y+w*w+A*X*x*y+B*X*x*w+C*X*y*w+B*C*x*y+A*C*x*w+A*B*y*w+A*B*C*X+A*A+B*B+C*C-7
    original=polynomial(source['characters'][0],(X,x,y,w),(p,q,r))
    values={A:p*q+1/(p*q),B:p*q*r+1/(p*q*r),C:p*q*q*r+1/(p*q*q*r)}
    assert s.expand(original-full.subs(values))==0
    L=A*y+B*w;M=C*y*w+A*B*C
    a=1-L*L/4;b=B*C*y+A*C*w-L*M/2;c=y*y+w*w+A*B*y*w+A*A+B*B+C*C-7-M*M/4
    residual=s.expand(a*x*x+b*x+c)
    assert s.expand(full.subs(X,-(L*x+M)/2)-residual)==0
    for v in (y,A,B,C,w):
        assert s.expand(s.diff(full,v).subs(X,-(L*x+M)/2)-s.diff(residual,v))==0
    P=s.expand(4*a*a*s.diff(c,y)-2*a*b*s.diff(b,y)+b*b*s.diff(a,y))
    assert P.coeff(y,5)==A**4/2
    assert P.subs({A:0,B:0,C:0,w:0})==8*y
    assert s.factor(P.subs({A:1,B:1,C:1,w:0}))==y*(y*y-3)*(y*y-5)/2
    assert s.gcd(s.Poly(a.subs({A:1,B:1,C:1,w:0}),y),s.Poly(b.subs({A:1,B:1,C:1,w:0}),y)).degree()==0
    assert s.expand(s.diff(a,y)+A*(A*y+B*w)/2)==0
    assert s.expand(s.diff(b,y).subs(A,0)+C*B*(w*w-2)/2)==0
    # Direct torus differentiation, independently of the trace-Hessian shortcut.
    tau=(U,V,W,t)
    on_torus=full.subs({A:U+1/U,B:V+1/V,C:W+1/W,w:t+1/t})
    point={X:-4,x:-4,y:-4,U:1,V:1,W:1,t:1}
    D=lambda f,v:v*s.diff(f,v)
    hessian1=s.Matrix([[D(D(on_torus,v),u).subs(point) for u in tau] for v in tau])
    assert hessian1==s.diag(-216,72,72,72)
    assert s.hessian(full,(X,x,y)).subs({X:-4,x:-4,y:-4,A:2,B:2,C:2,w:2}).det()==-216
    assert all(D(s.diff(on_torus,v),u).subs(point)==0 for v in (X,x,y) for u in tau)

    aa,bb=s.symbols('aa bb')
    original=polynomial(source['characters'][1],(aa,x,bb,y,w),(p,q))
    F=s.diff(original,aa)-bb;G=s.diff(original,bb)-aa;graph={aa:-G,bb:-F}
    reduced=s.expand(original.subs(graph,simultaneous=True))
    psi=x*x+y*y+w*w-x*x*y*y-x*x*w*w-y*y*w*w-x*y*w*(A*x+B*y+C*w)-A*B*x*y-A*C*x*w-B*C*y*w-A*B*C-3
    original_traces={A:q+1/q,B:p**3*q+1/(p**3*q),C:p**3*q*q+1/(p**3*q*q)}
    assert s.expand(reduced-psi.subs(original_traces))==0
    for v in (p,q):assert s.expand(D(original,v).subs(graph,simultaneous=True)-D(reduced,v))==0
    traces={A:U+1/U,B:V+1/V,C:U*V+1/(U*V)}
    phi=s.expand(psi.subs(traces))
    a=phi.coeff(x,2);b=phi.coeff(x,1);c=phi.coeff(x,0)
    P=s.expand(4*a*a*s.diff(c,y)-2*a*b*s.diff(b,y)+b*b*s.diff(a,y))
    assert s.expand(P.coeff(y,5)-2*((V-1/V)**2*w*w+4))==0
    assert s.Poly(s.expand(U**5*V**2*P.coeff(y,0)),U).degree()==10
    assert s.expand(U**5*V**2*P.coeff(y,0)).coeff(U,10)==-V**4*w**3
    spec={U:1,V:1,w:0}
    assert s.factor(P.subs(spec))==8*y*(y*y-3)*(y*y+1)
    assert a.subs(spec)==1-y*y and b.subs(spec)==-4*y
    y0=-(U+1/U)*w/2;E=(U*U-1)**2*w*w+4*U*U
    assert s.cancel(a.subs(y,y0)-E/(4*U*U))==0
    expected=-w*(U*U+1)*(U**4+1)*(V*V-1)/(2*U*U*V*(U*U-1))
    difference=s.fraction(s.cancel(b.subs(y,y0)-expected))[0]
    assert s.rem(difference,E,w)==0
    on_torus=phi.subs(w,t+1/t);tau=(U,V,t)
    point={x:-4,y:-4,U:1,V:1,t:1}
    assert all(s.diff(on_torus,v).subs(point)==0 for v in (x,y))
    hessian2=s.Matrix([[D(D(on_torus,v),u).subs(point) for u in tau] for v in tau])
    assert hessian2==s.Matrix([[144,-72,0],[-72,144,0],[0,0,72]])
    assert s.hessian(on_torus,(x,y)).subs(point).det()==-108
    assert all(D(s.diff(on_torus,v),u).subs(point)==0 for v in (x,y) for u in tau)
    assert s.expand(D(on_torus,t)-(t-1/t)*s.diff(phi,w).subs(w,t+1/t))==0
    assert s.expand(D(phi.subs({U:q,V:p**3*q},simultaneous=True),p)-3*D(phi,V).subs({U:q,V:p**3*q},simultaneous=True))==0
    assert s.expand(D(phi.subs({U:q,V:p**3*q},simultaneous=True),q)-(D(phi,U)+D(phi,V)).subs({U:q,V:p**3*q},simultaneous=True))==0
    return dict(original_characters_reconstructed=2,original_target_derivatives=True,
                primitive_and_squarefree_certificates=True,boundary_identities=True,
                direct_logarithmic_hessians=True,ramification_target_hyperplane=True)


def weyl(source):
    A=s.Matrix(cartan('E7'));inv=A.inv();positive=root_data(cartan('E7'))[1]
    roots=list(positive)+[tuple(-v for v in r) for r in positive]
    def subsystem(nodes):return {r for r in roots if all((inv*s.Matrix(r))[j]==0 for j in range(7) if j+1 not in nodes)}
    targets=[]
    for group in source['Weyl_classes']:
        for member in group['members']:
            current=subsystem(group['seed'])
            for node in member['word']:
                i=node-1;current={tuple(r[j]-r[i]*A[j,i] for j in range(7)) for r in current}
            assert current==subsystem(member['nodes']);targets.append(tuple(member['nodes']))
    rows=json.loads(gzip.decompress((ROOT/'unit-constraint-cases.json.gz').read_bytes()))
    expected={tuple(row['character']['nodes']) for row in rows if row['character']['type']=='E7'
              and sorted(row['character']['factors']) in (['A1','A1','A1'],['A1','A1','A2'])}
    assert set(targets)==expected and len(targets)==22
    return dict(Weyl_words_checked=22,entire_subsystems_compared=True,coverage_complete_within_two_cubic_classes=True)


def main():
    source=json.loads((ROOT/'binary-branching-verification.json').read_text())
    report=dict(result='passed',uses_Sage=False,roots=root_cover(),symbolic=identities(source),Weyl=weyl(source),
                generic_regular_fiber_dimensions=[6,10],new_standard_embeddings=22,
                resolved_proper_standard_Levis=123,remaining_embeddings=4,
                companion_memberships_verified=False,all_simple_adjoint_conjecture_proved=False,
                full_map_degree_claim=False,literature_priority_claim=False)
    names=['BINARY_BRANCHING_FIBERS.md','tools/audit_binary_branching_fibers.py',
           'tools/audit_affine_branching_certificates.py','tools/audit_branching_norm_fibers.py',
           'binary-branching-verification.json']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'binary-branching-audit.json').write_text(json.dumps(report,indent=2,default=int)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2,default=int))


if __name__=='__main__':main()
