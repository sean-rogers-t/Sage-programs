"""Native exact checks for the carrier theorem and smooth-relative-fiber route.

The geometric proofs are written in UNIFORM_ADJOINT_PROOF_TARGET.md.
These computations corroborate identities; they do not prove the remaining
adjoint bad-image bound or replace the arbitrary-rank root argument.
"""
from pathlib import Path
from hashlib import sha256
from collections import Counter
import gzip
import json
import sympy as s
from audit_affine_branching_certificates import cartan, root_data

ROOT = Path(__file__).resolve().parents[1]


def carriers():
    rows = json.loads(gzip.decompress((ROOT/'unit-constraint-cases.json.gz').read_bytes()))
    counts = Counter(); checks = []
    for record in rows:
        row = record['character']; C = s.Matrix(cartan(row['type']))
        positive = root_data(tuple(map(tuple,C.tolist())))[1]
        roots = set(positive) | {tuple(-v for v in a) for a in positive}
        I = [i-1 for i in row['nodes']]
        simples = [tuple(C[:,i]) for i in I]
        inverse = C.inv()
        for alpha in sorted(roots):
            if any(tuple(a+b for a,b in zip(alpha,beta)) in roots for beta in simples): continue
            labels = [alpha[i] for i in I]
            if sum(labels) != 3: continue
            assert all(v in (0,1) for v in labels)
            factors = [factor for factor in row['factor_orders'] if any(labels[j] for j in factor)]
            assert len(factors)==3
            sizes = sorted(map(len,factors))
            assert sizes[:2] == [1,1] or (sizes[:2]==[1,2] and 2<=sizes[2]<=4)
            for factor in factors:
                active = [j for j in factor if labels[j]]
                assert len(active)==1 and active[0] in (factor[0],factor[-1])
            # The joining root is -alpha; its inner products with the three
            # selected simple systems are exactly the three negative labels.
            indices = [I[j] for factor in factors for j in factor]
            columns = [C[:,i] for i in indices] + [-s.Matrix(alpha)]
            M = s.Matrix.hstack(*columns)
            gram = M.T*inverse*M
            assert M.rank()==len(indices)+1
            assert all(gram[i,i]==2 for i in range(gram.rows))
            assert all(gram[i,j] in (0,-1) for i in range(gram.rows) for j in range(i))
            if sizes[:2] == [1,1]:
                name='D'+str(sizes[2]+3); determinant=4
            else:
                name='E'+str(sizes[2]+4); determinant=5-sizes[2]
            assert gram.det()==determinant
            # Close under simple reflections independently of the stored
            # character polynomial. All roots of this carrier must be ambient.
            basis=[tuple(v) for v in columns]
            found=set(basis); queue=list(basis)
            dual=[inverse*v for v in columns]
            while queue:
                beta=queue.pop()
                for gamma,coroot in zip(basis,dual):
                    n=sum(a*b for a,b in zip(beta,coroot))
                    reflected=tuple(a-n*b for a,b in zip(beta,gamma))
                    assert reflected in roots
                    if reflected not in found: found.add(reflected); queue.append(reflected)
            rank=gram.rows
            expected=2*rank*(rank-1) if name.startswith('D') else {'E6':72,'E7':126,'E8':240}[name]
            assert len(found)==expected
            counts[name]+=1
            checks.append(dict(type=row['type'],nodes=row['nodes'],root=list(map(int,alpha)),
                               carrier=name,root_count=len(found)))
    r=s.Symbol('r',positive=True,integer=True)
    assert s.cancel(4*(r+1)*(2-(1+r/(r+1)))-4)==0
    assert s.cancel(6*(r+1)*(2-(s.Rational(1,2)+s.Rational(2,3)+r/(r+1)))-(5-r))==0
    return dict(existing_character_records=len(rows),cubic_highest_roots=len(checks),
                counts=dict(sorted(counts.items())),checks=checks,
                arbitrary_rank_determinant_identities=True,
                classification_proof_is_written_not_extrapolated=True)


def controls():
    a,b,x,z,w1,w2,c=s.symbols('a b x z w1 w2 c')
    D=lambda f:s.expand(z*s.diff(f,z))
    phi=a*b+(z+1/z)*x+z*z+z**-2
    field=s.QQ.frac_field(c)
    original=[b,a,z*z+1,s.expand(z*z*(D(phi)-c))]
    gb=s.groebner(original,a,b,x,z,w1,w2,domain=field)
    # z^2+1 implies z is already a unit; the displayed clearing loses none.
    expected=s.groebner([a,b,z*z+1,x+c*z/2],a,b,x,z,w1,w2,domain=field)
    assert gb==expected
    assert s.rem(s.Poly(2*z,z),s.Poly(z*z+1,z)).as_expr()==2*z
    assert s.gcd(2*z,z*z+1)==1
    H=s.hessian(phi,(a,b));L=s.Matrix([s.diff(phi,v).subs({a:0,b:0}) for v in (a,b)])
    determinant=H.det();N=s.expand(phi.subs({a:0,b:0})*determinant-(L.T*H.adjugate()*L)[0]/2)
    assert determinant==-1 and s.expand(N.subs(z,s.I)-(-2)*determinant)==0
    # On each E6 critical source stratum the affine fiber is smooth:
    # x,w1,w2 are free; H_aff has rank 2, B raises rank to 3.
    variables=(a,b,x,w1,w2)
    gradient=s.Matrix([s.diff(phi,v) for v in variables])
    hessian=gradient.jacobian(variables)
    mixed=gradient.applyfunc(D)
    at={a:0,b:0,z:s.I}
    assert hessian.subs(at).rank()==2
    assert hessian.row_join(mixed).subs(at).rank()==3
    # A non-quadratic smooth-excess control with no absent affine variable.
    u,v=s.symbols('u v');norm=(u*v-1)**2+z
    crit=s.groebner([s.diff(norm,u),s.diff(norm,v)],u,v)
    product=s.groebner([u*(u*v-1),v*(u*v-1)],u,v)
    assert crit==product
    assert s.expand((u*v-1)-u*v)==-1  # the two component ideals are comaximal
    assert s.hessian(norm,(u,v)).subs({u:0,v:0}).det()==-4
    assert s.cancel(s.hessian(norm,(u,v)).det().subs(v,1/u))==0
    assert s.cancel(s.hessian(norm,(u,v))[1,1].subs(v,1/u))==2*u*u
    # Reducedness alone fails even for reduced critical fibers: three axes.
    p,q,r=s.symbols('p q r');axes=p*q*r+z
    assert [s.diff(axes,t) for t in (p,q,r)]==[q*r,p*r,p*q]
    assert D(axes)==z
    assert s.hessian(axes,(p,q,r)).subs({p:0,q:0,r:0})==s.zeros(3)
    thick=(z-1)**2*x+z
    thick_ideal=s.groebner([(z-1)**2],z,x)
    assert thick_ideal.reduce(s.expand((D(thick)-1)**2))[1]==0
    assert thick_ideal.reduce(s.expand(D(thick)-1))[1]!=0
    # Persistent nilpotents from the already established altered-Levi local model.
    nil=s.groebner([u*u*v,u*v*v],u,v)
    assert nil.reduce(u*v)[1]==u*v and nil.reduce((u*v)**2)[1]==0
    # Special-target Sp4 bad image, with every fiber in the affine variable smooth.
    A=z*z+z**-2+1; delta=D(A)
    assert s.rem(s.Poly(s.expand(z**4*(delta**2+12)),z),s.Poly(z**4+z*z+1,z)).is_zero
    return dict(E6_original_generic_ideal_recovered=True,E6_excess_dimension=2,
        E6_source_fiber_tangent_dimension=3,E6_mixed_rank_increment=1,
        E6_entire_pencil_zero_on_dominant_source_stratum=True,
        nonlinear_smooth_excess_not_global_affine_space=True,
        nonreduced_incidence_smooth_relative_fiber_target_annihilator=True,
        reduced_singular_axes_control=True,axes_control_is_representation=False,
        altered_Levi_nilpotent_retained=True,Sp4_bad_target_polynomial='c^2+12',
        Sp4_special_failure_not_generic_failure=True)


def invariant_operator_control():
    C=s.Matrix(cartan('E6'));positive=root_data(tuple(map(tuple,C.tolist())))[1]
    roots=list(positive)+[tuple(-v for v in alpha) for alpha in positive]
    M=s.Matrix(roots);pairings=M*C.inv()*M.T
    values={}
    for degree in (2,4,6,8):
        eigenvalues={sum(v**degree for v in pairings.row(i)) for i in range(len(roots))}
        assert len(eigenvalues)==1
        eigenvalue=int(next(iter(eigenvalues)));values[str(degree)]=eigenvalue
        # Every diagonal entry of this operator is the same; the eigen-equation
        # imposes no relation between the 72 root multiplicities.
        assert all(sum(v**degree for v in pairings.row(i))==eigenvalue for i in range(len(roots)))
    return dict(root_moment_eigenvalues=values,all_operator_claim_has_written_proof=True,
                flat_invariant_operators_force_equal_root_multiplicities=False,
                radial_group_Casimir_claim=False)


def overlap_checks():
    """Exact identities, not a proof of the uniform overlap implication."""
    C=s.Matrix(cartan('E6')); inverse=C.inv()
    positive=root_data(tuple(map(tuple,C.tolist())))[1]
    roots=set(positive)|{tuple(-v for v in a) for a in positive}
    pairing=lambda a,b: (s.Matrix(a).T*inverse*s.Matrix(b))[0]
    zero=(0,)*6; derivatives={}
    for alpha in roots:
        Q=Counter({zero:2,alpha:2})
        Q.update(beta for beta in roots if pairing(beta,alpha)==1)
        product=dict(Q)
        for beta,n in Q.items():
            key=tuple(a-b for a,b in zip(beta,alpha))
            product[key]=product.get(key,0)-n
        product={k:v for k,v in product.items() if v}
        expected={beta:pairing(beta,alpha) for beta in roots if pairing(beta,alpha)}
        assert product==expected
        derivatives[alpha]=expected
    pairs=0
    for alpha in roots:
        for beta in roots:
            if pairing(alpha,beta)!=-1: continue
            gamma=tuple(a+b for a,b in zip(alpha,beta)); assert gamma in roots
            for weight in roots:
                assert derivatives[gamma].get(weight,0)==derivatives[alpha].get(weight,0)+derivatives[beta].get(weight,0)
            pairs+=1
    central=(2,3,4,6,4,2); retained=(0,1,2,4,5)
    added={a for a in roots if tuple(a[i] for i in retained) in ((0,1,0,0,0),(0,-1,0,0,0))
           and abs(sum(v*h for v,h in zip(a,central)))==3}
    assert len(added)==4
    witness=next((a,tuple(a[j]-a[3]*C[j,3] for j in range(6))) for a in sorted(added)
                 if tuple(a[j]-a[3]*C[j,3] for j in range(6)) not in added)
    assert witness[1] in roots

    # Actual SL4 adjoint, proper centralizer GL2 x (C*)^2 (determinant one).
    # Root sums give a separate construction from the Sage matrix-trace check.
    t=s.symbols('t0:4',nonzero=True);r=s.Symbol('r')
    mod=s.Poly(10*r*r+12*r+5,r)
    def nf(expr):
        num,den=s.fraction(s.cancel(expr))
        return s.rem(s.Poly(num,r)*s.invert(s.Poly(den,r),mod),mod).as_expr()
    point={t[0]:1,t[1]:2,t[2]:r,t[3]:1/(2*r)}
    def char(indices):
        return len(indices)-1+sum(t[i]/t[j] for i in indices for j in indices if i!=j)
    full=char(range(4)); c123=char((0,1,2));c124=char((0,1,3))
    assert s.expand(full-c123-c124+char((0,1))-char((2,3))+1)==0
    gradient=[nf((q*s.diff(full,q)).subs(point)) for q in t]
    assert gradient[0]==gradient[1]==-s.Rational(27,10)
    assert all(nf(gradient[i]-gradient[j])!=0 for i,j in ((0,2),(0,3),(2,3)))
    assert all(nf((t[i]-t[j]).subs(point))!=0 for i in range(4) for j in range(i))
    directional=lambda f:nf((t[0]*s.diff(f,t[0])-t[1]*s.diff(f,t[1])).subs(point))
    values=[directional(f) for f in (c123,c124,char((0,1)))]
    assert all(values) and nf(values[0]+values[1]-values[2])==0

    # Full original ideals in the two charts used by the arbitrary-m proof.
    charts=[]
    for m in (1,2):
        x,d,l=s.symbols('x d l');u=s.symbols(f'u0:{m}');v=s.symbols(f'v0:{m}')
        w=1-s.Add(*(a*b for a,b in zip(u,v)));P=x*x-d;f=w*P
        J=[s.diff(f,q) for q in (x,)+u+v]; variables=(l,x,d)+u+v
        gb=lambda fs:s.groebner(fs,*variables)
        assert gb(J+[l*w-1])==gb([x]+[d*q for q in u+v]+[l*w-1])
        assert gb(J+[l*u[0]-1])==gb([P,x*w,l*u[0]-1])
        charts.append(m)
    # The discriminant-nonzero hypothesis cannot be dropped.
    x,w,d,d1,d2=s.symbols('x w d d1 d2')
    jac=s.Matrix([[2*x,0,-d1,-d2],[w,x,0,0]])
    minors=[jac[:,[i,j]].det() for i in range(4) for j in range(i+1,4)]
    singular=s.groebner([x*x-d,x*w]+minors,x,w,d,d1,d2)
    assert all(singular.reduce(e)[1]==0 for e in (d,w*d1,w*d2))
    x,w=s.symbols('x w');nil=s.groebner([x*x,x*w],x,w)
    assert nil.reduce(x)[1]!=0 and nil.reduce(x*x)[1]==0

    # Exact coefficient alignment in the existing E6 control, before centering.
    a,b,u,v,x,z,ap,bp,X,k=s.symbols('a b u v x z ap bp X k')
    S=z**3+z**-3;A=z**-2+z*x;B=z*z+x/z
    phi=a*b+u*v+x*x+k*S*x+A*a*v+B*b*u-2
    changed=s.expand(phi.subs({a:ap-B*u,b:bp-A*v},simultaneous=True))
    assert s.expand(changed-(ap*bp+(1-u*v)*x*(x+S)+(k-1)*S*x-2))==0
    assert s.expand(changed.subs({k:1,x:X-S/2})-(ap*bp+(1-u*v)*(X*X-S*S/4)-2))==0
    assert s.cancel((z*s.diff(S*S/4,z))/(S)-s.Rational(3,2)*(z**3-z**-3))==0
    return dict(root_string_identities=72,A2_directional_syzygies=pairs,
        SL4_proper_centralizer_overlap_cancellation=True,
        SL4_unnormalized_target=[str(q) for q in gradient],
        SL4_subsystem_derivatives=[str(q) for q in values],
        aligned_quadratic_original_ideal_chart_checks=charts,
        aligned_quadratic_singular_scheme_target_factorization=True,
        zero_discriminant_nilpotent_control=True,
        E6_full_coefficient_alignment=True,altered_linear_coefficient_breaks_alignment=True,
        E6_reflection_forces_k_equals_one=[list(map(int,a)) for a in witness],
        arbitrary_m_bad_image_theorem_is_written_not_computational=True,
        general_root_overlap_implies_alignment=False,
        general_remaining_bad_image_bound_proved=False)


def main():
    root=carriers(); print('carrier checks:',root['counts'],flush=True)
    checked=controls()
    report=dict(result='passed',uses_Sage=False,uses_modular_specialization=False,
        carriers=root,controls=checked,invariant_operators=invariant_operator_control(),
        overlap=overlap_checks(),
        geometric_argument_verified_by_computer=False,
        uniform_adjoint_conjecture_proved=False,new_Levi_coverage=False,novelty_claim=False,
        root_carriers_force_squarefree_schur_pencils=False,
        smooth_relative_fiber_bad_image_exclusion_proved_in_note=True,
        remaining_locus='Sing(X) intersect {points singular in their actual relative critical fiber}',
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
            for p in (Path(__file__),ROOT/'UNIFORM_ADJOINT_PROOF_TARGET.md',ROOT/'tools/audit_affine_branching_certificates.py')},
        archive_sha256=sha256((ROOT/'unit-constraint-cases.json.gz').read_bytes()).hexdigest())
    (ROOT/'adjoint-proof-paths-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('carriers','source_sha256')},indent=2))


if __name__=='__main__':main()
