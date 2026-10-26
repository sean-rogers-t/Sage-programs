"""Sage root-lattice and ideal check, separate from the native weight audit."""
from pathlib import Path
from collections import Counter
from hashlib import sha256
import gzip
import json
from sage.all import RootSystem, PolynomialRing, QQ, matrix, NumberField, diagonal_matrix

ROOT=Path(__file__).resolve().parents[1]


def check_overlap_independently():
    lattice=RootSystem(['E',6]).root_lattice();roots=list(lattice.roots())
    zero=lattice.zero(); deriv={}
    for alpha in roots:
        Q={zero:QQ(2),alpha:QQ(2)}
        for beta in roots:
            if beta.scalar(alpha.associated_coroot())==1:
                Q[beta]=Q.get(beta,0)+1
        product=dict(Q)
        for beta,n in Q.items():product[beta-alpha]=product.get(beta-alpha,0)-n
        product={k:v for k,v in product.items() if v}
        deriv[alpha]={beta:beta.scalar(alpha.associated_coroot()) for beta in roots
                      if beta.scalar(alpha.associated_coroot())}
        assert product==deriv[alpha]
    pairs=0
    for alpha in roots:
        for beta in roots:
            if alpha.scalar(beta.associated_coroot())!=-1:continue
            assert all(deriv[alpha+beta].get(g,0)==deriv[alpha].get(g,0)+deriv[beta].get(g,0) for g in roots)
            pairs+=1
    simple=lattice.simple_roots();central=(2,3,4,6,4,2);retained=(1,2,3,5,6)
    labels=lambda a:tuple(a.scalar(simple[i].associated_coroot()) for i in retained)
    charge=lambda a:sum(n*a.scalar(simple[i].associated_coroot()) for i,n in enumerate(central,1))
    added={a for a in roots if labels(a) in ((0,1,0,0,0),(0,-1,0,0,0)) and abs(charge(a))==3}
    assert len(added)==4
    crossing=[(a,a-a.scalar(simple[4].associated_coroot())*simple[4]) for a in added]
    assert any(b not in added and b in roots for a,b in crossing)
    R=PolynomialRing(QQ,'r');r=R.gen();F=NumberField(10*r*r+12*r+5,'q');q=F.gen()
    g=diagonal_matrix(F,[1,2,q,1/(2*q)]);h=g*g.inverse().trace()-g.inverse()*g.trace()
    assert g.det()==1 and h[0,0]==h[1,1]==-QQ(27)/10
    assert all(h[i,i]!=h[j,j] for i,j in ((0,2),(0,3),(2,3)))
    assert all(g[i,i]!=g[j,j] for i in range(4) for j in range(i))
    def subderivative(indices):
        block=diagonal_matrix(F,[g[i,i] for i in indices])
        H=diagonal_matrix(F,[1,-1]+[0]*(len(indices)-2))
        return (block*H).trace()*block.inverse().trace()-block.trace()*(H*block.inverse()).trace()
    ds=[subderivative(I) for I in ((0,1,2),(0,1,3),(0,1))]
    assert all(ds) and ds[0]+ds[1]-ds[2]==0
    charts=[]
    for m in (1,2):
        P=PolynomialRing(QQ,names=['x','d','l']+[f'u{i}' for i in range(m)]+[f'v{i}' for i in range(m)])
        x,d,l=P.gens()[:3];u=P.gens()[3:3+m];v=P.gens()[3+m:]
        w=1-sum(a*b for a,b in zip(u,v));p=x*x-d;f=w*p
        J=[f.derivative(a) for a in (x,)+u+v]
        assert P.ideal(J+[l*w-1])==P.ideal([x]+[d*a for a in u+v]+[l*w-1])
        assert P.ideal(J+[l*u[0]-1])==P.ideal([p,x*w,l*u[0]-1])
        charts.append(m)
    P=PolynomialRing(QQ,names='x,w,d,d1,d2');x,w,d,d1,d2=P.gens()
    jac=matrix(P,[[2*x,0,-d1,-d2],[w,x,0,0]])
    singular=P.ideal([x*x-d,x*w]+jac.minors(2))
    assert all(singular.reduce(e)==0 for e in (d,w*d1,w*d2))
    P=PolynomialRing(QQ,names='x,w');x,w=P.gens();J=P.ideal(x*x,x*w)
    assert J.reduce(x)!=0 and J.reduce(x*x)==0
    C=PolynomialRing(QQ,'z').fraction_field();z=C.gen()
    P=PolynomialRing(C,names='a,b,u,v,x,ap,bp,X,k');a,b,u,v,x,ap,bp,X,k=P.gens()
    S=z**3+z**-3;A=z**-2+z*x;B=z*z+x/z
    phi=a*b+u*v+x*x+k*S*x+A*a*v+B*b*u-2
    completed=phi(a=ap-B*u,b=bp-A*v)
    assert completed==ap*bp+(1-u*v)*x*(x+S)+(k-1)*S*x-2
    assert completed(k=1,x=X-S/2)==ap*bp+(1-u*v)*(X*X-S*S/4)-2
    native=json.loads((ROOT/'adjoint-proof-paths-audit.json').read_text())['overlap']
    assert native['root_string_identities']==len(roots) and native['A2_directional_syzygies']==pairs
    assert native['aligned_quadratic_original_ideal_chart_checks']==charts
    return dict(root_string_identities=len(roots),A2_directional_syzygies=pairs,
        SL4_matrix_trace_and_exact_number_field_check=True,
        aligned_quadratic_original_ideal_chart_checks=charts,
        aligned_quadratic_singular_scheme_target_factorization=True,
        E6_full_coefficient_alignment=True,zero_discriminant_nilpotent_control=True,
        E6_reflection_forces_k_equals_one=True,
        general_remaining_bad_image_bound_proved=False)


def main():
    rows=json.loads(gzip.decompress((ROOT/'unit-constraint-cases.json.gz').read_bytes()))
    counts=Counter(); cached={}
    for record in rows:
        row=record['character']; kind=row['type']
        if kind not in cached:
            lattice=RootSystem([kind[0],int(kind[1:])]).root_lattice()
            cached[kind]=(list(lattice.roots()),lattice.simple_roots(),matrix(QQ,lattice.cartan_type().cartan_matrix()))
        roots,simple,C=cached[kind]; rootset=set(roots)
        I=row['nodes']
        for alpha in roots:
            if any(alpha+simple[i] in rootset for i in I): continue
            labels=[alpha.scalar(simple[i].associated_coroot()) for i in I]
            if sum(labels)!=3: continue
            factors=[f for f in row['factor_orders'] if any(labels[j] for j in f)]
            assert len(factors)==3 and all(v in (0,1) for v in labels)
            for factor in factors:
                active=[j for j in factor if labels[j]]
                assert len(active)==1 and active[0] in (factor[0],factor[-1])
            sizes=sorted(map(len,factors))
            selected=[simple[I[j]] for factor in factors for j in factor]+[-alpha]
            gram=matrix(QQ,[[a.scalar(b.associated_coroot()) for b in selected] for a in selected])
            assert gram.is_positive_definite()
            if sizes[:2]==[1,1]:
                name='D'+str(sizes[2]+3); det=4
            else:
                assert sizes[:2]==[1,2] and 2<=sizes[2]<=4
                name='E'+str(sizes[2]+4); det=5-sizes[2]
            assert gram.det()==det
            generated=set(selected); queue=list(selected)
            while queue:
                beta=queue.pop()
                for gamma in selected:
                    reflection=beta-beta.scalar(gamma.associated_coroot())*gamma
                    assert reflection in rootset
                    if reflection not in generated: generated.add(reflection);queue.append(reflection)
            rank=gram.nrows()
            assert len(generated)==(2*rank*(rank-1) if name[0]=='D' else {'E6':72,'E7':126,'E8':240}[name])
            counts[name]+=1
    native=json.loads((ROOT/'adjoint-proof-paths-audit.json').read_text())
    assert native['result']=='passed' and dict(sorted(counts.items()))==native['carriers']['counts']
    roots=cached['E6'][0];moments={}
    for degree in (2,4,6,8):
        values={sum(alpha.scalar(beta.associated_coroot())**degree for beta in roots) for alpha in roots}
        assert len(values)==1
        moments[str(degree)]=int(next(iter(values)))
    assert moments==native['invariant_operators']['root_moment_eigenvalues']
    L=PolynomialRing(QQ,'c').fraction_field();c=L.gen()
    P=PolynomialRing(L,names='a,b,x,z,w1,w2');a,b,x,z,w1,w2=P.gens()
    phi=a*b+(z+1/z)*x+z*z+z**-2
    original=P.ideal([P(phi.derivative(v).numerator()) for v in (a,b,x,w1,w2)]+[P((z*phi.derivative(z)-c).numerator())])
    expected=P.ideal([a,b,z*z+1,2*x+c*z])
    assert original==expected
    assert P.ideal([z*z+1,2*z]).is_one()
    Q=PolynomialRing(QQ,names='u,v,z');u,v,z=Q.gens()
    f=(u*v-1)**2+z
    J=Q.ideal([f.derivative(u),f.derivative(v)])
    assert J==Q.ideal([u,v]).intersection(Q.ideal(u*v-1))
    assert Q.ideal([u,v,u*v-1]).is_one()
    A=PolynomialRing(QQ,names='p,q,r');p,q,r=A.gens()
    axes=A.ideal([p*q,p*r,q*r])
    assert axes==A.ideal([p,q]).intersection(A.ideal([p,r])).intersection(A.ideal([q,r]))
    assert axes.dimension()==1
    nil=Q.ideal([u*u*v,u*v*v])
    assert nil.reduce(u*v)!=0 and nil.reduce((u*v)**2)==0
    thick=(z-1)**2*u+z;mu=z*thick.derivative(z);thick_ideal=Q.ideal((z-1)**2)
    assert thick_ideal.reduce((mu-1)**2)==0 and thick_ideal.reduce(mu-1)!=0
    report=dict(result='passed',engine='Sage root lattice and Singular ideals',
        independent_of_native_weight_coordinates=True,carrier_counts=dict(sorted(counts.items())),
        flat_invariant_operator_control=moments,geometric_proof_is_written=True,
        overlap=check_overlap_independently(),
        E6_actual_generic_ideal_equal=True,nonlinear_smooth_excess_decomposition_verified=True,
        reduced_singular_axes_decomposition_verified=True,altered_Levi_nilpotent_verified=True,
        nonreduced_incidence_target_annihilator_verified=True,
        uniform_adjoint_conjecture_proved=False,new_Levi_coverage=False,
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                       for p in (Path(__file__),ROOT/'UNIFORM_ADJOINT_PROOF_TARGET.md',ROOT/'adjoint-proof-paths-audit.json')})
    (ROOT/'adjoint-proof-paths-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
