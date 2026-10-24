"""Exact structural checks for exceptional adjoint contact-Levi fibers."""
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json

from sage.all import (QQ, ZZ, CartanType, RootSystem, PolynomialRing,
                      CyclotomicField, matrix, vector, prod)
from sage.env import SAGE_VERSION
from verify_affine_branching_classification import levi_data, character_weights

ROOT=Path(__file__).resolve().parents[1]
CONTACT=[('G2',2,'A1',4),('F4',1,'C3',14),('E6',2,'A5',20),
         ('E7',1,'D6',32),('E8',8,'E7',56)]


def serialize(weights):
    return [dict(weight=list(w),multiplicity=int(m)) for w,m in sorted(weights.items())]


def graded_characters():
    rows=[]
    for kind,node,factor,dimV in CONTACT:
        n=int(kind[1:]); I=tuple(i for i in range(n) if i!=node-1)
        components,_=levi_data(kind,I)
        assert len(components)==1 and components[0][0]==factor
        order=tuple(I[j] for j in components[0][1])
        A=CartanType(kind).cartan_matrix()
        H=A.transpose().inverse()*vector(QQ,[int(i==node-1) for i in range(n)])
        assert all(v in ZZ for v in H) and H[node-1]==2
        R=RootSystem(kind).ambient_space(); coroots=R.simple_coroots()
        grades={k:Counter() for k in (-2,-1,0,1,2)}
        roots=[]
        for alpha in R.roots():
            w=tuple(int(alpha.scalar(coroots[j+1])) for j in range(n))
            roots.append(w)
            grade=int(vector(QQ,w).dot_product(H)); assert grade in grades
            grades[grade][tuple(w[j] for j in order)]+=1
        grades[0][(0,)*(n-1)]+=n
        local=RootSystem(factor).ambient_space(); lc=local.simple_coroots()
        adj=Counter(tuple(int(alpha.scalar(lc[j+1])) for j in range(n-1))
                    for alpha in local.roots())
        adj[(0,)*(n-1)]+=n
        assert grades[0]==adj
        assert grades[2]==grades[-2]==Counter({(0,)*(n-1):1})
        assert grades[1]==grades[-1] and sum(grades[1].values())==dimV
        candidates={'A1':[(3,)],'C3':[(0,0,1)],'A5':[(0,0,1,0,0)],
                    'D6':[(0,0,0,0,1,0),(0,0,0,0,0,1)],
                    'E7':[(0,0,0,0,0,0,1)]}[factor]
        matches=[hw for hw in candidates if Counter(character_weights(factor,hw))==grades[1]]
        assert len(matches)==1
        trace=sum(int(vector(QQ,w).dot_product(H))**2 for w in roots)
        assert trace==2*dimV+8
        cover=matrix(ZZ,n,n,lambda i,j:int(H[i]) if j==node-1 else int(i==j))
        assert abs(cover.det())==2
        rows.append(dict(type=kind,omitted_node=node,derived_type=factor,
                         derived_standard_nodes=[j+1 for j in order],
                         grading_cocharacter=list(map(int,H)),V_highest=list(matches[0]),
                         V_dimension=dimV,trace_normalization=trace,central_cover_degree=2,
                         grades={str(k):serialize(v) for k,v in grades.items()}))
        print('Contact grading checked:',kind,factor,flush=True)
    return rows


def exact_ideals():
    F=PolynomialRing(QQ,'c').fraction_field();c=F.gen()
    R=PolynomialRing(F,['x1','x2','x3','x4','x5','z','u'])
    x1,x2,x3,x4,x5,z,u=R.gens()
    chi=x1*x5+z*z+u*u+(z+u)*x3
    eq=[chi.derivative(x) for x in (x1,x2,x3,x4,x5)]
    eq += [z*chi.derivative(z)-u*chi.derivative(u)-c,z*u-1]
    I=R.ideal(eq)
    J=R.ideal([x1,x5,z*z+1,u+z,2*z*x3-c])
    assert I==J and I.dimension()==2
    nonzero=[f for f in eq if f]
    jac=matrix(R,[[f.derivative(v) for v in R.gens()] for f in nonzero])
    assert I+R.ideal(jac.minors(5))==R.ideal(1)
    deck=R.hom([-x1,x2,-x3,x4,-x5,-z,-u],R)
    assert R.ideal([deck(f) for f in I.gens()])==I
    assert deck(x2)==x2 and deck(x4)==x4

    G=PolynomialRing(F,['x','z','u']);x,z,u=G.gens()
    ff=2*x+(z+u)*(3*x*x-2)
    mu=(z-u)*(2*(z+u)+x**3-2*x)
    curve=G.ideal([ff,z*u-1])
    curve_j=matrix(G,[[f.derivative(v) for v in G.gens()] for f in curve.gens()])
    assert curve+G.ideal(curve_j.minors(2))==G.ideal(1)
    target=G.ideal([ff,z*u-1,mu-c])
    assert target.dimension()==0 and target!=G.ideal(1)
    target_j=matrix(G,[[f.derivative(v) for v in G.gens()] for f in target.gens()])
    assert target+G.ideal(target_j.det())==G.ideal(1)
    return dict(E6_exact_ideal=True,E6_quotient_dimension_on_cover=2,
                E6_cover_critical_scheme_smooth=True,E6_deck_action_checked=True,
                E6_free_coordinates=['x2','x4'],
                G2_relative_curve_smooth=True,G2_generic_fiber_nonempty_finite_smooth=True)


def e6_root_trace_audit():
    """Original adjoint torus eigenvalues; does not evaluate the character polynomial."""
    C=CyclotomicField(12); eta=C.gen(); ii=eta**3; omega=eta**4
    eigen=[C(2),2*omega,2*omega**2,C(1)/2,omega/2,omega**2/2]
    assert prod(eigen)==1 and sum(eigen)==0 and sum(1/e for e in eigen)==0
    R=RootSystem('E6').ambient_space(); coroots=R.simple_coroots()
    roots=[vector(ZZ,[alpha.scalar(coroots[j+1]) for j in range(6)]) for alpha in R.roots()]
    H=vector(ZZ,[1,2,2,3,2,1]); order=[0,2,3,4,5]
    q=[prod(eigen[:j+1]) for j in range(5)]
    parameters=[ii**H[j] for j in range(6)]
    for j,value in zip(order,q): parameters[j]*=value
    values=[prod(parameters[j]**w[j] for j in range(6)) for w in roots]
    assert all(v!=1 for v in values)
    B=matrix(QQ,6,6,lambda i,j:sum(w[i]*w[j] for w in roots))
    gradient=vector(C,[sum(value*w[j] for value,w in zip(values,roots)) for j in range(6)])
    theta=B.change_ring(C).inverse()*gradient
    target=C(65)*ii/192
    assert theta==target*H and H*B*H==48
    cartan_hessian=matrix(C,6,6,lambda i,j:sum(value*w[i]*w[j] for value,w in zip(values,roots)))
    assert cartan_hessian.rank()==4
    outside=sum(w.dot_product(H)!=0 for w in roots)
    assert outside==42
    assert 6+sum(w.dot_product(H)==0 for w in roots)==36
    # At a regular torus point root directions are controlled by equivariance.
    root_derivatives=[w.dot_product(theta)/(value-1) for w,value in zip(roots,values)]
    assert sum(v!=0 for v in root_derivatives)==42
    assert cartan_hessian.rank()+42==46
    return dict(method='root-basis adjoint eigenvalues and equivariant root derivatives',
                root_eigenvalues=72,source_regular_semisimple=True,
                cyclotomic_field=12,SL6_eigenvalues=[str(e) for e in eigen],
                central_z=str(ii),target_scalar=str(target),trace_normalization=48,
                normalized_target_checked=True,target_centralizer_dimension=36,
                Cartan_differential_rank=4,root_direction_rank=42,
                full_differential_rank=46,fiber_tangent_dimension=32)


def main():
    grades=graded_characters();ideals=exact_ideals();point=e6_root_trace_audit()
    paths=['EXCEPTIONAL_ADJOINT_CONTACT_FIBERS.md','tools/verify_exceptional_adjoint_contact.py',
           'tools/verify_affine_branching_classification.py','GENERAL_CAYLEY_FIBER_CRITICAL.md']
    report=dict(result='passed',date='2026-10-23',sage_version=SAGE_VERSION,
                contact_gradings=grades,exact_ideals=ideals,E6_root_trace_audit=point,
                E6_nonzero_contact_fiber_regular_dimension=32,
                E6_actual_quotient_critical_scheme='A^2',
                all_exceptional_Levis_classified=False,
                classical_dimension_formula_extends_unchanged=False,
                public_API_changed=False,literature_priority_claim=False,
                source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in paths})
    (ROOT/'exceptional-adjoint-contact-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(result='passed',contact_gradings=len(grades),exact_ideals=ideals,
                         E6_root_trace_audit=point),indent=2))


if __name__=='__main__':main()
