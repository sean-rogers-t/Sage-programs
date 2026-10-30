"""Exact corroboration of the complete simple-group regularity table.

Run with Sage Python. Root-lattice/weight-lattice checks, abstract endpoint
patterns, independent character constructions, exact number-algebra witnesses,
and rank-one controls corroborate a written proof; they do not exhaust irreps.
"""
from collections import Counter, defaultdict, deque
from hashlib import sha256
from itertools import combinations, product
from math import comb, gcd
from pathlib import Path
import json

from sage.all import (QQ, ZZ, CartanMatrix, PolynomialRing, RootSystem,
                      WeylCharacterRing, matrix, prod, vector)
from sage.env import SAGE_VERSION
from verify_type_a_regularity import gt_weights

ROOT=Path(__file__).resolve().parents[1]


def root_data(kind):
    RS=RootSystem(kind)
    Q=RS.root_lattice()
    P=RS.weight_lattice()
    simples=[tuple(int(v) for v in P(a).to_vector()) for a in Q.simple_roots()]
    roots=[tuple(int(v) for v in P(a).to_vector()) for a in Q.roots()]
    positive=[tuple(int(v) for v in P(a).to_vector()) for a in Q.positive_roots()]
    root_coordinates={tuple(int(v) for v in P(a).to_vector()):tuple(int(v) for v in a.to_vector())
                      for a in Q.roots()}
    highest=max(Q.positive_roots(),key=lambda a:sum(a.to_vector()))
    return dict(kind=kind,rank=RS.cartan_type().rank(),simple=simples,roots=roots,
                positive=positive,root_coordinates=root_coordinates,
                highest=tuple(int(v) for v in P(highest).to_vector()))


def plus(*vectors):
    return tuple(sum(row) for row in zip(*vectors))


def lattice_checks():
    kinds=(['A'+str(n) for n in range(2,8)]+
           ['B'+str(n) for n in range(3,9)]+
           ['C'+str(n) for n in range(2,9)]+
           ['D'+str(n) for n in range(4,10)]+['G2','F4','E6','E7','E8','A1','B2'])
    rows=[]
    for kind in kinds:
        d=root_data(kind)
        nonprimitive=[]
        for a in d['roots']:
            g=0
            for v in a:
                g=gcd(g,abs(v))
            if g!=1:
                nonprimitive.append(a)
        exceptional=kind.startswith('C') or kind in ['A1','B2']
        assert bool(nonprimitive)==exceptional
        row=dict(type=kind,rank=d['rank'],roots=len(d['roots']),nonprimitive_roots=len(nonprimitive))
        if kind.startswith('C'):
            assert len(nonprimitive)==2*d['rank']
        if kind in ['A1','B2']:
            assert len(nonprimitive)==(2 if kind=='A1' else 4)
        if kind.startswith('B') and kind!='B2' or kind in ['G2','F4']:
            S=set(d['roots'])
            pair=next((i,j) for i in range(d['rank']) for j in range(d['rank'])
                      if i!=j and plus(d['simple'][i],d['simple'][j]) in S
                      and plus(d['simple'][i],d['simple'][i],d['simple'][j]) in S)
            row['root_string_simple_nodes']=[v+1 for v in pair]
        if kind.startswith('D') or kind in ['E6','E7','E8']:
            C=CartanMatrix(kind)
            center=next(i for i in range(d['rank']) if sum(C[i,j]!=0 for j in range(d['rank']) if j!=i)>=3)
            leaves=[j for j in range(d['rank']) if j!=center and C[center,j]!=0][:3]
            S=set(d['roots'])
            assert all(plus(d['simple'][center],d['simple'][j]) in S for j in leaves)
            assert all(plus(d['simple'][center],d['simple'][i],d['simple'][j]) in S
                       for i,j in combinations(leaves,2))
            row['fork_nodes']=[center+1]+[j+1 for j in leaves]
        rows.append(row)
    return rows


def binomial(*terms):
    result=Counter()
    for coefficient,monomial in terms:
        result[monomial]+=coefficient
    return {v:c for v,c in result.items() if c}


def linear_combination(*terms):
    return binomial(*[(a*c,v) for a,p in terms for v,c in p.items()])


def is_unit_binomial(p):
    return len(p)==2 and sorted(p.values())==[-1,1]


def endpoint_patterns():
    zero=(0,0)
    a=(1,0)
    b=(0,1)
    f=binomial((1,a),(-1,zero))
    branches=[binomial((1,zero),(-1,(0,-1))),binomial((1,plus(a,b)),(-1,a))]
    for g in branches:
        assert is_unit_binomial(linear_combination((1,f),(1,g)))
        assert not is_unit_binomial(linear_combination((2,f),(1,g)))
    zero=(0,0,0,0)
    center=(1,0,0,0)
    f0=binomial((1,center),(-1,zero))
    passing=0
    for orientations in product([0,1],repeat=3):
        leaves=[]
        for j,orientation in enumerate(orientations,1):
            b=tuple(int(i==j) for i in range(4))
            leaves.append(binomial((1,zero),(-1,tuple(-v for v in b))) if orientation==0
                          else binomial((1,plus(center,b)),(-1,center)))
        assert all(is_unit_binomial(linear_combination((1,f0),(1,g))) for g in leaves)
        passing+=all(is_unit_binomial(linear_combination((1,f0),(1,leaves[i]),(1,leaves[j])))
                     for i,j in combinations(range(3),2))
    assert passing==0
    return dict(root_string_paths=2,fork_orientations=8,passing_fork_orientations=0)


def character_weights(kind,labels,rings):
    if kind not in rings:
        rings[kind]=WeylCharacterRing(kind,style='coroots')
    c=rings[kind](tuple(labels))
    coroots=RootSystem(kind).ambient_space().simple_coroots()
    n=len(labels)
    w={tuple(int(mu.inner_product(coroots[i+1])) for i in range(n)):int(m)
       for mu,m in c.weight_multiplicities().items()}
    assert sum(w.values())==c.degree()
    return w


def add_weights(*terms):
    w=defaultdict(int)
    for a,D in terms:
        for mu,m in D.items():
            w[mu]+=a*m
    return {v:m for v,m in w.items() if m}


def orbit(highest,d):
    seen={tuple(highest)}
    queue=deque(seen)
    while queue:
        w=queue.popleft()
        for i,a in enumerate(d['simple']):
            v=tuple(x-w[i]*y for x,y in zip(w,a))
            if v not in seen:
                seen.add(v)
                queue.append(v)
    return seen


def gram(w):
    n=len(next(iter(w)))
    return matrix(QQ,n,n,lambda i,j:sum(m*mu[i]*mu[j] for mu,m in w.items()))


def positivity_walls(w,d):
    a=d['simple'][0]
    G=gram(w)
    row=vector(QQ,a)*G.inverse()
    assert row[0]>0 and all(row[j]==0 for j in range(1,d['rank']))
    checks=0
    for b in d['positive']:
        pairs=[(i,j) for i,j in combinations(range(d['rank']),2) if a[i]*b[j]!=a[j]*b[i]]
        if not pairs:
            continue
        i,j=pairs[0]
        h=[0]*d['rank']
        h[i],h[j]=b[j],-b[i]
        assert sum(x*y for x,y in zip(b,h))==0
        target=sum(x*y for x,y in zip(a,h))
        assert target!=0
        for sign in [-1,1]:
            value=sum(QQ(m*mu[0])*QQ(2)**(sign*sum(x*y for x,y in zip(mu,h))) for mu,m in w.items())
            assert value*target*sign>0
            checks+=1
    return checks


def root_witness(w,d):
    """An exact number-algebra point with every root !=1 and one output root 0."""
    n=d['rank']
    primes=[2,3,5,7,11,13,17,19,23,29]
    U=PolynomialRing(QQ,'X')
    X=U.gen()
    trials=0
    for moment in range(n):
        nonzero={mu:m*mu[moment] for mu,m in w.items() if mu[moment]}
        assert nonzero
        for variable in range(n):
            exponents=[mu[variable] for mu in nonzero]
            shift=min(exponents)
            if max(exponents)==shift:
                continue
            trials+=1
            fixed=[QQ(p) for p in primes[:n]]
            p=sum(c*prod(fixed[j]**mu[j] for j in range(n) if j!=variable)*X**(mu[variable]-shift)
                  for mu,c in nonzero.items())
            if p==0:
                continue
            q=U(p)
            while q.degree()>0 and q(0)==0:
                q=q//X
            for beta in d['positive']:
                c=prod(fixed[j]**beta[j] for j in range(n) if j!=variable)
                e=beta[variable]
                boundary=c*X**e-1 if e>=0 else c-X**(-e)
                assert boundary!=0
                while q.degree()>0:
                    common=q.gcd(boundary)
                    if common.degree()==0:
                        break
                    q=q//common
            if q.degree()<1:
                continue
            q=U(q.factor()[0][0]).monic()
            assert q(0)!=0
            A=U.quotient(q,'xx')
            t=[A(v) for v in fixed]
            t[variable]=A.gen()
            for beta in d['positive']:
                assert prod(v**e for v,e in zip(t,beta))!=1
            values={mu:prod(v**e for v,e in zip(t,mu)) for mu in w}
            rhs=vector(A,[sum(m*mu[i]*values[mu] for mu,m in w.items()) for i in range(n)])
            assert rhs[moment]==0
            image=gram(w).inverse().change_ring(A)*rhs
            root=d['simple'][moment]
            assert sum(v*c for v,c in zip(root,image))==0
            return dict(vanishing_simple_root=moment+1,variable_fundamental_character=variable+1,
                        fixed_fundamental_characters=[None if j==variable else str(v) for j,v in enumerate(fixed)],
                        root_polynomial=str(q),number_algebra_degree=int(q.degree()),
                        specializations_tried=trials,regular_input_roots_checked=len(d['positive']),
                        full_trace_Gram_solve_checked=True)
    raise AssertionError('No witness on finite prime-coordinate slices: '+d['kind'])


def A_checks(rings):
    rows=[]
    controls=[]
    mixed=[]
    for N in range(3,7):
        kind='A'+str(N-1)
        d=root_data(kind)
        for a in product(range(3),repeat=N-1):
            if not 0<sum(a)<=2:
                continue
            p=tuple(sum(a[i:]) for i in range(N-1))+(0,)
            gt={tuple(mu[i]-mu[i+1] for i in range(N-1)):int(m) for mu,m in gt_weights(p).items()}
            w=character_weights(kind,a,rings)
            assert w==gt
            support=[mu for mu in w if sum(mu)!=0]
            good=(len(support)==2 and tuple(x-y for x,y in zip(*support)) in
                  [d['highest'],tuple(-v for v in d['highest'])])
            endpoint=(sum(a)==1 and (a[0]==1 or a[-1]==1))
            assert good==endpoint
            rows.append(dict(N=N,labels=list(a),dimension=sum(w.values()),
                             highest_root_nonorthogonal_weight_support=len(support),endpoint=endpoint))
        defining=character_weights(kind,(1,)+(0,)*(N-2),rings)
        dual=character_weights(kind,(0,)*(N-2)+(1,),rings)
        w=add_weights((1,defining),(1,dual))
        mixed.append(dict(N=N,**root_witness(w,d)))
        z=[QQ(j+2) for j in range(N-1)]
        z.append(1/prod(z))
        t=[prod(z[:j+1]) for j in range(N-1)]
        for orientation,ww in [('defining',defining),('dual',dual)]:
            ww=add_weights((3,ww),(2,{(0,)*(N-1):1}))
            values={mu:prod(v**e for v,e in zip(t,mu)) for mu in ww}
            rhs=vector(QQ,[sum(m*mu[i]*values[mu] for mu,m in ww.items()) for i in range(N-1)])
            c=gram(ww).inverse()*rhs
            image=[c[0]]+[c[j]-c[j-1] for j in range(1,N-1)]+[-c[-1]]
            expected=[v-sum(z)/N for v in z] if orientation=='defining' else [-1/v+sum(1/u for u in z)/N for v in z]
            assert image==expected
            controls.append(dict(N=N,orientation=orientation,copy_multiplicity=3,trivial_multiplicity=2))
    return rows,controls,mixed


def remaining_type_checks(rings):
    cases=[]
    for n in [4,5,8]:
        kind='D'+str(n)
        d=root_data(kind)
        choices=[('vector',(1,)+(0,)*(n-1)),
                 ('half_spin_minus',(0,)*(n-2)+(1,0)),
                 ('half_spin_plus',(0,)*(n-1)+(1,)),('adjoint',d['highest'])]
        for name,a in choices:
            cases.append((kind,name,[(1,a)]))
        cases.append((kind,'vector_plus_half_spin',[(1,choices[0][1]),(1,choices[2][1])]))
    for kind,minimal in [('G2',(1,0)),('F4',(0,0,0,1)),('E6',(1,0,0,0,0,0)),
                         ('E7',(0,0,0,0,0,0,1))]:
        d=root_data(kind)
        cases.extend([(kind,'minimal',[(1,minimal)]),(kind,'adjoint',[(1,d['highest'])]),
                      (kind,'minimal_plus_adjoint',[(1,minimal),(1,d['highest'])])])
        if kind=='E6':
            cases.append((kind,'dual_minimal',[(1,(0,0,0,0,0,1))]))
    d=root_data('E8')
    cases.extend([('E8','adjoint',[(1,d['highest'])]),
                  ('E8','two_adjoint_plus_trivial',[(2,d['highest']),(1,(0,)*8)])])
    rows=[]
    for kind,name,terms in cases:
        d=root_data(kind)
        summands=[]
        for m,a in terms:
            w=character_weights(kind,a,rings)
            if a==d['highest']:
                assert w=={**{beta:1 for beta in d['roots']},(0,)*d['rank']:d['rank']}
            elif any(a):
                independent={mu:1 for mu in orbit(a,d)}
                if kind=='G2':
                    independent[(0,)*d['rank']]=1
                if kind=='F4':
                    independent[(0,)*d['rank']]=2
                assert w==independent
            summands.append((m,w))
        w=add_weights(*summands)
        walls=positivity_walls(w,d)
        witness=root_witness(w,d)
        rows.append(dict(type=kind,representation=name,dimension=sum(w.values()),
                         positive_wall_trace_checks=walls,**witness))
        print(kind,name,'exact regular-to-nonregular witness passed',flush=True)
    return rows


def rank_one_checks():
    R=PolynomialRing(QQ,'x')
    x=R.gen()
    characters=[R(1),x]
    for j in range(2,9):
        characters.append(x*characters[-1]-characters[-2])
    rows=[]
    for m in range(1,9):
        w={k:comb(2*m,m+k) for k in range(-m,m+1)}
        multiplicities={j:w.get(j,0)-w.get(j+2,0) for j in range(m+1)}
        assert all(a>=0 for a in multiplicities.values())
        assert sum(a*characters[j] for j,a in multiplicities.items())==(x+2)**m
        S=sum(QQ(a)*k*k for k,a in w.items())
        assert S==2*m*4**(m-1)
        for z in [QQ(2),QQ(3),QQ(-2),QQ(1)/3]:
            y=sum(QQ(a)*k*z**k for k,a in w.items())/S
            expected=QQ(2)/4**m*(z+2+1/z)**m*(z-1)/(z+1)
            assert y==expected and y!=0
        g=matrix(QQ,[[-1,1],[0,-1]])
        I=matrix.identity(QQ,2)
        image=QQ(2)/4**m*(g+I).det()**(m-1)*(g-I)*(g+I).adjugate()
        assert (image.rank()==1)==(m==1)
        assert (image==0)==(m>=2)
        rows.append(dict(tensor_power=m,character=str((x+2)**m),
                         trivial_multiplicity=multiplicities[0],boundary_image_rank=image.rank()))
    return rows


def main():
    roots=lattice_checks()
    patterns=endpoint_patterns()
    rings={}
    a,positive,mixed=A_checks(rings)
    other=remaining_type_checks(rings)
    rank_one=rank_one_checks()
    paths=['REGULARITY_CLASSIFICATION.md','tools/verify_regularity_classification.py',
           'tools/verify_type_a_regularity.py','MIXED_REGULARITY.md','mixed-regularity-verification.json']
    report=dict(result='passed',date='2026-10-22',sage_version=SAGE_VERSION,
        arbitrary_rank_proof='REGULARITY_CLASSIFICATION.md, Sections 2-7',
        scope='All connected complex simple groups and nontrivial actual finite-dimensional rational representations; finite checks corroborate a written proof, not an exhaustive representation search.',
        all_simply_connected_rows_classified=True,central_quotients_corollary=True,
        root_system_cases=len(roots),root_lattice_weight_lattice_checks=sum(r['roots'] for r in roots),
        abstract_endpoint_patterns=patterns,A_independent_GT_Freudenthal_cases=len(a),
        A_positive_reducible_trace_controls=len(positive),A_mixed_exact_witnesses=len(mixed),
        remaining_classical_exceptional_exact_witnesses=len(other),
        positive_root_wall_trace_checks=sum(r['positive_wall_trace_checks'] for r in other),
        rank_one_determinant_character_controls=len(rank_one),
        root_systems=roots,A_character_cases=a,A_positive_controls=positive,A_mixed_witnesses=mixed,
        other_exact_witnesses=other,rank_one_controls=rank_one,
        novelty_claim=False,externally_reviewed=False,runtime_API_changed=False,
        source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'regularity-classification-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in
        ('source_sha256','root_systems','A_character_cases','A_positive_controls','A_mixed_witnesses','other_exact_witnesses','rank_one_controls')},indent=2))


if __name__=='__main__':
    main()
