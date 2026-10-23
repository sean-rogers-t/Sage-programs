"""Exact scheme, character and rank checks for the quadratic fiber theorem."""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json

from sage.all import (QQ, PolynomialRing, LaurentPolynomialRing, Partitions,
                      RootSystem, matrix, vector, prod)
from sage.env import SAGE_VERSION
from verify_affine_branching_classification import (character_weights, levi_data,
                                                    subsets, blocks_from_edges)

ROOT=Path(__file__).resolve().parents[1]


def scheme_controls():
    rows=[]
    def check(name,ell,s,phi,expected,dim,nonreduced=False):
        F=PolynomialRing(QQ,['c'+str(i) for i in range(s)]).fraction_field()
        c=F.gens()
        R=PolynomialRing(F,['x'+str(i) for i in range(ell)]+
                         ['z'+str(i) for i in range(s)]+['u'+str(i) for i in range(s)])
        x=R.gens()[:ell]; z=R.gens()[ell:ell+s]; u=R.gens()[ell+s:]
        chi=phi(x,z)
        generators=[chi.derivative(v) for v in x]
        generators += [z[j]*chi.derivative(z[j])-c[j] for j in range(s)]
        generators += [z[j]*u[j]-1 for j in range(s)]
        I=R.ideal(generators); J=R.ideal(expected(x,z,u,c))
        assert I==J,(name,I.groebner_basis(),J.groebner_basis())
        assert I.dimension()==dim,(name,I.dimension(),dim)
        if nonreduced: assert x[0] not in I and x[0]**2 in I
        rows.append(dict(name=name,exact_ideal_equality=True,dimension=int(dim),
                         nonreduced_control=nonreduced,
                         reduced_groebner_basis=[str(g) for g in I.groebner_basis()]))
    check('rank_drop_supplies_generic_fiber',1,1,
          lambda x,z:(z[0]-1)*x[0]**2/2,
          lambda x,z,u,c:[z[0]-1,x[0]**2-2*c[0],u[0]-1],0)
    check('nonreduced_incidence_excluded_generically',1,1,
          lambda x,z:(z[0]-1)**2*x[0]+z[0],
          lambda x,z,u,c:[1],-1)
    check('quadratic_double_rank_zero_excluded',1,1,
          lambda x,z:(z[0]-1)**2*x[0]**2/2+z[0],
          lambda x,z,u,c:[x[0],z[0]-c[0],u[0]-1/c[0]],0)
    check('unused_variable_positive_excess',2,1,
          lambda x,z:x[0]**2/2+z[0],
          lambda x,z,u,c:[x[0],z[0]-c[0],u[0]-1/c[0]],1)
    check('variable_kernel_square_completion',2,1,
          lambda x,z:(x[1]+(z[0]-1)*x[0])**2/2+z[0],
          lambda x,z,u,c:[x[1]+(c[0]-1)*x[0],z[0]-c[0],u[0]-1/c[0]],1)
    check('quadratic_and_affine_compatibility',2,1,
          lambda x,z:x[0]**2/2+(z[0]-1)*x[1],
          lambda x,z,u,c:[x[0],z[0]-1,x[1]-c[0],u[0]-1],0)
    check('two_parameter_rank_drop',2,2,
          lambda x,z:sum((z[i]-1)*x[i]**2/2 for i in range(2)),
          lambda x,z,u,c:[z[0]-1,z[1]-1,x[0]**2-2*c[0],x[1]**2-2*c[1],u[0]-1,u[1]-1],0)
    check('cubic_persistent_double_structure',1,1,
          lambda x,z:x[0]**3+z[0],
          lambda x,z,u,c:[x[0]**2,z[0]-c[0],u[0]-1/c[0]],0,True)
    return rows


def type_a_partitions():
    rows=[]
    for N in range(2,9):
        for partition in Partitions(N):
            b=list(map(int,partition))
            if len(b)==1: continue
            ell=sum(v-1 for v in b); count=len(b)
            names=['z'+str(i) for i in range(count)]+['x'+str(i) for i in range(ell)]
            R=LaurentPolynomialRing(QQ,names); z=R.gens()[:count]; x=R.gens()[count:]
            a=R(0); bb=R(0); u=[R(0)]*ell; v=[R(0)]*ell; offset=0
            images={}; eigenvalues=[]
            for i,size in enumerate(b):
                if size==1:
                    a+=z[i]; bb+=1/z[i]; eigenvalues.append(z[i]); continue
                coords=x[offset:offset+size-1]
                local=list(coords)+[1/prod(coords)]
                for j in range(1,size):
                    images[x[offset+j-1]]=sum(prod(c) for c in combinations(local,j))
                eigenvalues += [z[i]*w for w in local]
                u[offset]=z[i]; v[offset+size-2]=1/z[i]; offset+=size-1
            A=a+sum(u[i]*x[i] for i in range(ell))
            B=bb+sum(v[i]*x[i] for i in range(ell)); chi=A*B-1
            Q=matrix(R,ell,ell,lambda i,j:u[i]*v[j]+v[i]*u[j])
            assert matrix(R,ell,ell,lambda i,j:chi.derivative(x[i]).derivative(x[j]))==Q
            assert all(chi.derivative(x[i])==u[i]*B+v[i]*A for i in range(ell))
            # Independent root character: N-1 zero weights and every e_i-e_j.
            root_character=N-1+sum(eigenvalues[i]/eigenvalues[j] for i in range(N) for j in range(N) if i!=j)
            assert chi.subs(images)==root_character
            t=sum(size==2 for size in b); p=b.count(1)
            expected_rank=2 if any(size>=3 for size in b) or t>=2 else t
            actual_rank=Q.change_ring(R.fraction_field()).rank()
            assert actual_rank==expected_rank
            if all(size<=2 for size in b) and t:
                ones={zi:R(1) for zi in z}
                assert matrix(QQ,ell,ell,lambda i,j:Q[i,j].subs(ones)).rank()==1
                if p:
                    double=[i for i,size in enumerate(b) if size==2]
                    exponent_rows=[[2*int(j==i)-2*int(j==double[0]) for j in range(count)] for i in double[1:]]
                    exponent_rows.append(b)
                    assert matrix(QQ,exponent_rows).rank()==t
            rows.append(dict(N=N,blocks=b,semisimple_rank=ell,generic_quadratic_rank=int(actual_rank),
                             generic_regular_fiber_nonempty=all(size<=2 for size in b) and p>0,
                             dimension_if_nonempty=2*t if all(size<=2 for size in b) and p>0 else None,
                             independent_root_character_identity=True))
    return rows


def dominance_checks():
    rows=[]
    for t in range(1,5):
        for p in range(1,7):
            N=2*t+p; R=LaurentPolynomialRing(QQ,['v'+str(i) for i in range(p)])
            v=R.gens(); H=sum(1/y-y for y in v)/(4*N)
            nu=[H*(y+1/y) for y in v]; point={y:R(2) for y in v}
            J=matrix(QQ,p,p,lambda i,j:(v[j]*nu[i].derivative(v[j])).subs(point))
            expected=matrix(QQ,p,p,lambda i,j:-QQ(9*p*int(i==j)+25)/(16*N))
            assert J==expected and J.det()!=0
            assert J.det()==(-QQ(9*p)/(16*N))**(p-1)*(-QQ(34*p)/(16*N))
            rows.append(dict(double_blocks=t,singletons=p,determinant=str(J.det())))
    return rows


def convolve(A,B):
    out=Counter()
    for a,m in A.items():
        for b,n in B.items(): out[tuple(x+y for x,y in zip(a,b))]+=m*n
    return dict(out)


@lru_cache(None)
def fundamental_local(kind,I,pos):
    components,_=levi_data(kind,I)
    factor,order=next((f,o) for f,o in components if pos in o)
    hw=tuple(int(j==pos) for j in order); out={}
    for w,m in character_weights(factor,hw).items():
        full=[0]*len(I)
        for j,c in zip(order,w): full[j]=c
        out[tuple(full)]=m
    return out


def root_restriction_polynomial(kind,I):
    R=RootSystem(kind).ambient_space(); n=R.dimension()
    coroots=R.simple_coroots(); rank=len(coroots); weights=Counter({(0,)*len(I):rank})
    for root in R.roots():
        weights[tuple(int(root.scalar(coroots[i+1])) for i in I)]+=1
    if not I: return {():sum(weights.values())}
    _,height=levi_data(kind,I); polynomial={}
    while weights:
        dominant=[w for w in weights if all(v>=0 for v in w)]
        assert dominant
        nu=max(dominant,key=lambda w:(sum(h*v for h,v in zip(height,w)),w))
        multiplicity=weights[nu]; polynomial[nu]=polynomial.get(nu,0)+multiplicity
        assert sum(nu)<=3,(kind,I,nu)
        product_weights={(0,)*len(I):1}
        for pos,exponent in enumerate(nu):
            for _ in range(exponent): product_weights=convolve(product_weights,fundamental_local(kind,I,pos))
        assert product_weights[nu]==1
        for w,m in product_weights.items():
            weights[w]-=multiplicity*m
            if not weights[w]: del weights[w]
    return {w:m for w,m in polynomial.items() if m}


def residual_blocks(family,n,I):
    nodes={i+1 for i in I}
    if family in ('B','C'):
        m=0
        if n in nodes:
            m=1
            while n-m in nodes: m+=1
        return m,blocks_from_edges(n-m,nodes)
    if n-1 in nodes and n in nodes:
        m=2
        while n-m in nodes: m+=1
        return m,blocks_from_edges(n-m,nodes)
    edges={j for j in nodes if j<n-1}
    if n-1 in nodes or n in nodes: edges.add(n-1)
    return 0,blocks_from_edges(n,edges)


def classical_transfer():
    rows=[]
    for family,start in [('B',2),('C',2),('D',4)]:
        for n in range(start,7):
            kind=family+str(n)
            for I in subsets(n):
                polynomial=root_restriction_polynomial(kind,I)
                degree=max(map(sum,polynomial))
                m,blocks=residual_blocks(family,n,I); t=sum(b>=2 for b in blocks)
                excluded=(family=='B' and m==1 and t>0) or (family=='D' and m==2 and t>0)
                assert (degree==3)==excluded,(kind,I,m,blocks,polynomial)
                rows.append(dict(type=kind,levi_nodes=[i+1 for i in I],residual_rank=m,
                                 blocks=blocks,character_degree=degree,quadratic=not excluded,
                                 polynomial=[dict(exponent=list(w),coefficient=c) for w,c in sorted(polynomial.items())]))
            print('Root-restriction transfer checked:',kind,flush=True)
    return rows


def main():
    controls=scheme_controls(); print('Scheme controls passed',flush=True)
    a=type_a_partitions(); dominance=dominance_checks(); print('Type A formulas passed',flush=True)
    classical=classical_transfer()
    report=dict(result='passed',date='2026-10-23',sage_version=SAGE_VERSION,
                scheme_controls=controls,type_A_block_partitions=a,dominance_checks=dominance,
                classical_root_restrictions=classical,
                quadratic_theorem_requires_constant_rank=False,
                generic_smoothness_includes_rank_drop_loci=True,
                type_A_nonemptiness_proof_is_all_rank=True,
                classical_transfer_nonemptiness_classified=False,
                cubic_control_is_actual_representation=False,
                literature_priority_claim=False,public_API_changed=False,
                source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for name in ['QUADRATIC_BRANCHING_FIBERS.md',
                                            'tools/verify_quadratic_branching_fibers.py',
                                            'tools/verify_affine_branching_classification.py',
                                            'GENERAL_CAYLEY_FIBER_CRITICAL.md']})
    (ROOT/'quadratic-branching-fibers-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(result='passed',scheme_controls=len(controls),type_A_partitions=len(a),
                         dominance_checks=len(dominance),classical_Levis=len(classical),
                         cubic_exclusions=sum(not r['quadratic'] for r in classical)),indent=2))


if __name__=='__main__': main()
