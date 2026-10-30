"""Exact restriction tests for joint affine Levi branching (research tool)."""
from collections import Counter, deque
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json
import argparse

from sage.all import (QQ, CartanType, CartanMatrix, RootSystem, WeylCharacterRing,
                      matrix, vector)
from sage.env import SAGE_VERSION

ROOT=Path(__file__).resolve().parents[1]


@lru_cache(None)
def character_weights(kind, highest):
    W=WeylCharacterRing(kind,style='coroots'); R=RootSystem(kind).ambient_space()
    coroots=R.simple_coroots(); indices=list(R.index_set())
    chi=W(tuple(highest))
    return {tuple(int(w.scalar(coroots[j])) for j in indices):int(mult)
            for w,mult in chi.weight_multiplicities().items()}


def subsets(r):
    for size in range(r):
        yield from combinations(range(r),size)


@lru_cache(None)
def levi_data(kind,I):
    A=CartanType(kind).cartan_matrix(); M=A.matrix_from_rows_and_columns(I,I)
    remaining=set(range(len(I))); components=[]
    while remaining:
        seed=min(remaining); comp={seed}; queue=[seed]; remaining.remove(seed)
        while queue:
            i=queue.pop()
            new=[j for j in remaining if M[i,j] or M[j,i]]
            remaining.difference_update(new); comp.update(new); queue.extend(new)
        comp=tuple(sorted(comp)); C=CartanMatrix(M.matrix_from_rows_and_columns(comp,comp))
        ct=C.cartan_type()
        # Sage's inferred Cartan type carries a verified standard-to-local map.
        base=getattr(ct,'_type',ct)
        relabel=getattr(ct,'_relabelling',dict(zip(base.index_set(),C.index_set())))
        local_positions={j:i for i,j in enumerate(C.index_set())}
        order=tuple(comp[local_positions[relabel[j]]] for j in base.index_set())
        standard=CartanType(base).cartan_matrix()
        assert all(M[order[i],order[j]]==standard[i,j] for i in range(len(comp)) for j in range(len(comp)))
        name=str(base.type())+str(base.rank())
        components.append((name,order))
    height=tuple(sum(M.inverse()[i,j] for i in range(len(I))) for j in range(len(I))) if I else ()
    return tuple(components),height


def restriction_test(kind,weights,I):
    I=tuple(I)
    if not I:
        return dict(affine=True,levi=[],constituents=[dict(highest=[],multiplicity=sum(weights.values()))])
    components,height=levi_data(kind,I)
    residual=Counter()
    for w,mult in weights.items(): residual[tuple(w[i] for i in I)]+=mult
    constituents=[]
    while residual:
        dominant=[w for w in residual if all(x>=0 for x in w)]
        assert dominant,'a positive character must have a dominant highest weight'
        nu=max(dominant,key=lambda w:(sum(h*x for h,x in zip(height,w)),w))
        mult=residual[nu]
        if sum(nu)>1:
            return dict(affine=False,levi=[c[0] for c in components],
                        forbidden_highest=list(nu),multiplicity=mult,
                        previously_removed=constituents)
        constituents.append(dict(highest=list(nu),multiplicity=mult))
        if not any(nu): small={(0,)*len(I):1}
        else:
            pos=nu.index(1)
            factor,order=next((f,o) for f,o in components if pos in o)
            hw=tuple(nu[j] for j in order)
            small={}
            for w,m in character_weights(factor,hw).items():
                full=[0]*len(I)
                for j,val in zip(order,w): full[j]=val
                small[tuple(full)]=m
        for w,m in small.items():
            residual[w]-=mult*m
            assert residual[w]>=0,(kind,I,nu,w,residual[w])
            if not residual[w]: del residual[w]
    return dict(affine=True,levi=[c[0] for c in components],constituents=constituents)


def blocks_from_edges(length,edges):
    blocks=[]; size=1
    for j in range(1,length):
        if j in edges: size+=1
        else: blocks.append(size); size=1
    if length: blocks.append(size)
    return blocks


def classical_prediction(family,n,k,I):
    nodes={j+1 for j in I}
    if family=='A':
        blocks=blocks_from_edges(n+1,nodes)
        return k in (1,n) or sum(m>=2 for m in blocks)<=1
    if family in ('B','C'):
        m=0
        if n in nodes:
            m=1
            while n-m in nodes: m+=1
        blocks=blocks_from_edges(n-m,nodes); t=sum(b>=2 for b in blocks)
        if family=='C': return k==1 or t==0
        if k==n: return t+int(m>0)<=1
        return (k==1 or t==0) and (m==0 or m>k)
    assert family=='D'
    if n-1 in nodes and n in nodes:
        m=2
        while n-m in nodes: m+=1
        blocks=blocks_from_edges(n-m,nodes); flip=0
    else:
        m=0; flip=int(n in nodes)
        edges={j for j in nodes if j<n-1}
        if n-1 in nodes or n in nodes: edges.add(n-1)
        blocks=blocks_from_edges(n,edges)
    t=sum(b>=2 for b in blocks)
    if k<=n-2: return (k==1 or t==0) and (m==0 or m>=k+2)
    if m>=2: return t==0
    if t<=1: return True
    parity=(n+int(k==n-1)+flip)%2
    # Only GL2 x GL2 has two active factors but no odd-parity allocation.
    return blocks==[2,2] and parity==1


def orbit_test(weights,I):
    # Exact for minuscule/quasi-minuscule representations only.
    witnesses=[list(w) for w in weights if any(w) and all(w[i]>=0 for i in I)
               and sum(w[i] for i in I)>1]
    return not witnesses, witnesses[:1]


def extremal_obstructions(kind,k):
    """Witness nonaffinity using Weyl words; no full character is needed."""
    A=CartanType(kind).cartan_matrix(); n=A.nrows()
    columns=[tuple(int(A[j,i]) for j in range(n)) for i in range(n)]
    start=tuple(int(j==k-1) for j in range(n)); parents={start:None}; queue=deque([start])
    todo={sum(1<<j for j in I):I for I in subsets(n) if I}; bad={}
    def word(w):
        out=[]
        while parents[w] is not None:
            prev,i=parents[w]; out.append(i+1); w=prev
        return list(reversed(out))
    processed=0
    while queue and todo:
        w=queue.popleft(); processed+=1
        negative=sum(1<<i for i,v in enumerate(w) if v<0)
        positive=sum(1<<i for i,v in enumerate(w) if v>0)
        large=sum(1<<i for i,v in enumerate(w) if v>=2)
        for mask,I in list(todo.items()):
            if not mask&negative and (mask&large or (mask&positive).bit_count()>=2):
                bad[I]=dict(ambient_extremal_weight=list(w),weyl_word=word(w),
                            forbidden_highest=[w[i] for i in I])
                del todo[mask]
        for i,col in enumerate(columns):
            if not w[i]: continue
            v=tuple(w[j]-w[i]*col[j] for j in range(n))
            if v not in parents: parents[v]=(w,i); queue.append(v)
    return bad,[()]+list(todo.values()),dict(processed_weights=processed,
                                           orbit_completed=not bool(queue))


def probe_exceptional_extremals():
    for kind in ['G2','F4','E6','E7','E8']:
        for k in range(1,CartanType(kind).rank()+1):
            bad,remaining,stats=extremal_obstructions(kind,k)
            print(kind,k,'negative',len(bad),'remaining',[[i+1 for i in I] for I in remaining],stats,flush=True)


def classical_checks(max_rank):
    rows=[]
    for family,start in [('A',2),('B',2),('C',2),('D',4)]:
        for n in range(start,max_rank+1):
            kind=family+str(n)
            for k in range(1,n+1):
                hw=tuple(int(j==k-1) for j in range(n)); weights=character_weights(kind,hw)
                for I in subsets(n):
                    actual=restriction_test(kind,weights,I)
                    expected=classical_prediction(family,n,k,I)
                    assert actual['affine']==expected,(kind,k,I,expected,actual)
                    rows.append(dict(type=kind,fundamental_node=k,levi_nodes=[i+1 for i in I],
                                     affine=expected,derived_factors=actual['levi'],
                                     forbidden_highest=actual.get('forbidden_highest')))
            print('Classical classification passed:',kind,flush=True)
    return rows


def exceptional_prediction(kind,k,I):
    if not I: return True
    factors=[f for f,_ in levi_data(kind,tuple(I))[0]]
    if kind=='G2': return k==1 and tuple(I)==(1,)
    if kind=='F4':
        if k==1: return tuple(I)==(0,1,2)
        if k==4: return tuple(I) in [(0,),(1,),(0,1),(1,2),(0,1,2),(1,2,3)]
        return False
    if kind=='E6':
        if k in (1,6): return len(factors)==1
        return k==2 and len(factors)==1 and factors[0] in ('D4','D5')
    if kind=='E7':
        if k==7: return len(factors)==1
        if k==2: return tuple(I)==tuple(range(6))
        return k==1 and len(factors)==1 and factors[0] in ('D4','D5','D6','E6')
    if kind=='E8':
        if k==1: return tuple(I)==tuple(range(7))
        return k==8 and len(factors)==1 and factors[0] in ('D4','D5','D6','D7','E6','E7')
    raise ValueError(kind)


def weyl_dimension(kind,k):
    R=RootSystem(kind).ambient_space(); fw=R.fundamental_weights(); rho=sum(fw.values())
    out=QQ(1)
    for root in R.positive_roots():
        coroot=root.associated_coroot(); out*=(rho+fw[k]).scalar(coroot)/rho.scalar(coroot)
    return int(out)


def verify_extremal_word(kind,k,cert):
    R=RootSystem(kind).ambient_space(); start=R.fundamental_weights()[k]; w=start
    roots=R.simple_roots(); coroots=R.simple_coroots()
    for j in cert['weyl_word']: w=w-w.scalar(coroots[j])*roots[j]
    labels=[int(w.scalar(coroots[j])) for j in R.index_set()]
    assert labels==cert['ambient_extremal_weight']
    assert w.scalar(w)==start.scalar(start)
    return labels


def norm_gap_certificate(kind,k,I):
    R=RootSystem(kind).ambient_space(); roots=R.simple_roots(); coroots=R.simple_coroots()
    n=CartanType(kind).rank(); A=matrix(QQ,[[roots[j+1].scalar(coroots[i+1]) for j in range(n)] for i in range(n)])
    ai=A.inverse(); lam=sum(ai[j,k-1]*roots[j+1] for j in range(n)); bound=lam.scalar(lam)
    M=A.matrix_from_rows_and_columns(I,I); mi=M.inverse()
    fundamental=[sum(mi[j,i]*roots[I[j]+1] for j in range(len(I))) for i in range(len(I))]
    assert all(fundamental[i].scalar(coroots[I[j]+1])==int(i==j) for i in range(len(I)) for j in range(len(I)))
    Gram=matrix(QQ,[[u.scalar(v) for v in fundamental] for u in fundamental])
    assert all(a>=0 for a in Gram.list())
    minimum=min(Gram[i,i]+2*Gram[i,j]+Gram[j,j] for i in range(len(I)) for j in range(i,len(I)))
    assert bound<minimum,(kind,k,I,bound,minimum)
    return dict(ambient_highest_norm_squared=str(bound),minimum_nonfundamental_norm_squared=str(minimum),
                levi_fundamental_Gram=[[str(x) for x in row] for row in Gram.rows()])


def complete_exceptional_checks():
    records=[]; summaries=[]; negative_total=0; positive_total=0
    small_orbits={('G2',1),('F4',4),('E6',1),('E6',2),('E6',6),('E7',1),('E7',7),('E8',8)}
    for kind in ['G2','F4','E6','E7','E8']:
        n=CartanType(kind).rank()
        for k in range(1,n+1):
            bad,remaining,stats=extremal_obstructions(kind,k)
            assert set(remaining)=={I for I in subsets(n) if exceptional_prediction(kind,k,I)}
            grouped={}
            for I,cert in bad.items():
                key=tuple(cert['ambient_extremal_weight'])
                if key not in grouped:
                    labels=verify_extremal_word(kind,k,cert)
                    grouped[key]=dict(ambient_extremal_weight=labels,weyl_word=cert['weyl_word'],
                                      covered_levi_node_sets=[])
                labels=grouped[key]['ambient_extremal_weight']
                assert all(labels[i]>=0 for i in I) and sum(labels[i] for i in I)>1
                grouped[key]['covered_levi_node_sets'].append([i+1 for i in I])
            dimension=weyl_dimension(kind,k); good=[]
            for I in remaining:
                if not I:
                    good.append(dict(levi_nodes=[],affine=True,levi=[],
                                     constituents=[dict(highest=[],multiplicity=dimension)]))
                    continue
                weights=character_weights(kind,tuple(int(j==k-1) for j in range(n)))
                assert sum(weights.values())==dimension
                actual=restriction_test(kind,weights,I); assert actual['affine']
                row=dict(levi_nodes=[i+1 for i in I],norm_gap=norm_gap_certificate(kind,k,I),**actual)
                if (kind,k) in small_orbits:
                    assert orbit_test(weights,I)[0]
                    row['independent_minuscule_or_quasiminuscule_test']=True
                good.append(row)
            negative_total+=len(bad); positive_total+=len(good)
            records.append(dict(type=kind,node=k,dimension=dimension,
                                negative_certificates=list(grouped.values()),
                                positive_restrictions=good,orbit_search=stats))
            summaries.append(dict(type=kind,node=k,dimension=dimension,affine_count=len(good),
                                  proper_standard_Levis=2**n-1,affine_node_sets=[x['levi_nodes'] for x in good]))
            print('Complete exceptional classification:',kind,k,'affine',len(good),'of',2**n-1,flush=True)
    return records,summaries,negative_total,positive_total


def independent_full_character_checks():
    count=0
    choices={'G2':[1,2],'F4':[1,2,3,4],'E6':list(range(1,7)),
             'E7':[1,2,6,7],'E8':[1,8]}
    for kind,ks in choices.items():
        n=CartanType(kind).rank()
        for k in ks:
            weights=character_weights(kind,tuple(int(j==k-1) for j in range(n)))
            for I in subsets(n):
                actual=restriction_test(kind,weights,I)
                assert actual['affine']==exceptional_prediction(kind,k,I),(kind,k,I,actual)
                count+=1
    return count


def built_in_branching_checks():
    rows=[]
    for source,k,target,terms,trivial in [
            ('F4',1,'B3',[(1,2),(2,1),(3,2)],1),
            ('E7',2,'E6',[(1,1),(2,2),(3,1),(5,1),(6,1)],0),
            ('E8',1,'E7',[(1,3),(2,2),(6,1),(7,2)],1)]:
        G=WeylCharacterRing(source,style='coroots'); S=WeylCharacterRing(target,style='coroots')
        n=CartanType(source).rank(); m=CartanType(target).rank()
        actual=G(tuple(int(i==k-1) for i in range(n))).branch(S,rule='levi')
        expected=S(trivial)
        for node,mult in terms: expected+=mult*S(tuple(int(i==node-1) for i in range(m)))
        assert actual==expected,(source,k,target,actual,expected)
        rows.append(dict(source=source,node=k,target=target,decomposition=str(actual)))
    return rows


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--max-rank',type=int,default=6)
    parser.add_argument('--explore',action='store_true'); args=parser.parse_args()
    classic=classical_checks(args.max_rank)
    exceptional,summary,negative,positive=complete_exceptional_checks()
    full_count=independent_full_character_checks(); native=built_in_branching_checks()
    report=dict(result='passed',date='2026-10-23',sage_version=SAGE_VERSION,
                classical_comparisons=classic,exceptional_certificates=exceptional,
                exceptional_summary=summary,all_rank_classical_proofs_in_note=True,
                exceptional_negative_pairs=negative,exceptional_affine_pairs=positive,
                independent_full_character_comparisons=full_count,built_in_branching_checks=native,
                E7_all_fundamentals_classified=True,E8_all_fundamentals_classified=True,
                nonaffinity_is_not_a_fiber_singularity_proof=True)
    if not args.explore:
        sources=['AFFINE_BRANCHING_CLASSIFICATION.md','tools/verify_affine_branching_classification.py']
        report['source_sha256']={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in sources}
    destination=ROOT/('affine-branching-classification-exploration.json' if args.explore else 'affine-branching-classification-verification.json')
    destination.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(result='passed',classical_pairs=len(classic),exceptional_pairs=negative+positive,
                         affine_pairs=positive,negative_pairs=negative,full_character_checks=full_count,
                         built_in_checks=len(native)),indent=2))


if __name__=='__main__': main()
