"""Fixed binary-arrangement tensor certificates by Sage polynomial substitution.

Internal proof tooling. Dimensions depend on the Frobenius offset, not on
representation rank. A finite certificate is not a general extension claim.
"""
from itertools import permutations, combinations, combinations_with_replacement
from sage.all import ZZ, PolynomialRing, matrix


def cubic_offset_obstruction():
    """A primitive compatible family that does not extend, k=r=3."""
    R = PolynomialRing(ZZ, names=['u0','u1','u2','u3','v0','v1','v2','v3'])
    u, v = R.gens()[:4], R.gens()[4:]
    subsets = [B for n in range(4) for B in combinations(range(1,4), n)]
    family = {B: R(0) for B in subsets}
    for B in subsets:
        if len(B) == 2:
            j, l = B
            i = next(i for i in range(1,4) if i not in B)
            family[B] = (u[j]*v[l]-u[l]*v[j])*(v[j]+v[l]-v[i])*(v[j]-v[l])
        elif len(B) == 3:
            family[B] = (sum(v[1:])*sum(u[i]*v[i]**2 for i in range(1,4))
                         -sum(u[1:])*sum(v[i]**3 for i in range(1,4)))
    for B, C in combinations(subsets, 2):
        coeff = {i: int(i in B)-int(i in C) for i in range(1,4)}
        pivot = next(i for i in coeff if coeff[i])
        sub = {u[pivot]: -sum(coeff[i]*u[i] for i in coeff if i != pivot)/coeff[pivot],
               v[pivot]: -sum(coeff[i]*v[i] for i in coeff if i != pivot)/coeff[pivot]}
        assert R((family[B]-family[C]).subs(sub)) == 0
    mons = [(i,jj) for i in range(4) for jj in combinations_with_replacement(range(4),3)]
    seen, ambient = set(), []
    for mon in mons:
        if mon in seen:
            continue
        orbit = set()
        for pm in permutations(range(1,4)):
            f = lambda x: 0 if x == 0 else pm[x-1]
            orbit.add((f(mon[0]),tuple(sorted(f(t) for t in mon[1]))))
        seen.update(orbit)
        ambient.append(sum(u[i]*R.prod(v[j] for j in jj) for i,jj in orbit))
    local = [u[i]*R.prod(v[j] for j in jj) for i in range(1,4)
             for jj in combinations_with_replacement(range(1,4),3)]
    records = []
    for proper in [False, True]:
        chosen = [B for B in subsets if not proper or len(B) < 3]
        A = matrix(ZZ, [[R(f.subs({u[0]:sum(u[i] for i in B),
                                  v[0]:sum(v[i] for i in B)})).monomial_coefficient(mon)
                         for f in ambient] for B in chosen for mon in local])
        W = matrix(ZZ, [[family[B].monomial_coefficient(mon)] for B in chosen for mon in local])
        diagonal = [abs(int(x)) for x in A.augment(W).elementary_divisors() if x]
        assert A.rank() == 20 and diagonal == [1]*21
        records.append(dict(proper_subsets_only=proper, ambient_rank=20,
                            augmented_rank=21, augmented_nonunits=[]))
    return dict(integer_pairwise_compatibility=True, records=records,
                obstruction_restricts_to_primitive_proper_obstruction=True)


def binary_tensor_certificate(k, r, proper=False):
    R=PolynomialRing(ZZ,names=['u'+str(i) for i in range(k+1)]+['v'+str(i) for i in range(k+1)])
    u=R.gens()[:k+1];v=R.gens()[k+1:]
    perms=list(permutations(range(1,k+1)))
    subsets=[tuple(B) for n in range(k if proper else k+1) for B in combinations(range(1,k+1),n)]
    mons=[(i,jj) for i in range(1,k+1) for jj in combinations_with_replacement(range(1,k+1),r)]
    def poly(mon):
     i,jj=mon
     return u[i]*R.prod(v[j] for j in jj)
    def permmon(mon,pm):return (pm[mon[0]-1],tuple(sorted(pm[j-1] for j in mon[1])))
    def orbits(items,ps):
     seen=set();out=[]
     for mon in items:
      if mon in seen:continue
      orb=set(permmon(mon,pm) for pm in ps);seen.update(orb);out.append(orb)
     return out
    families=[]
    for size in range(k if proper else k+1):
     B0=tuple(range(1,size+1))
     stab=[pm for pm in perms if set(pm[:size])==set(B0)]
     for orb in orbits(mons,stab):
      family={B:R(0) for B in subsets}
      for B in subsets:
       if len(B)!=size:continue
       pm=next(pm for pm in perms if set(pm[:size])==set(B))
       family[B]=sum(poly(permmon(mon,pm)) for mon in orb)
      families.append(family)
    constraints={}
    for B,C in combinations(subsets,2):
     coeff={i:int(i in B)-int(i in C) for i in range(1,k+1)}
     pivot=next(i for i in coeff if coeff[i])
     sub={u[pivot]:-sum(coeff[i]*u[i] for i in coeff if i!=pivot)/coeff[pivot],v[pivot]:-sum(coeff[i]*v[i] for i in coeff if i!=pivot)/coeff[pivot]}
     for col,family in enumerate(families):
      dif=R((family[B]-family[C]).subs(sub))
      for exp,c in dif.dict().items():constraints.setdefault((B,C,tuple(exp)),[0]*len(families))[col]+=c
    C=matrix(ZZ,list(constraints.values()))
    # Ambient invariant orbit sums, including fixed index 0.
    ambientmons=[(i,jj) for i in range(k+1) for jj in combinations_with_replacement(range(k+1),r)]
    def apm(mon,pm):
     f=lambda x:0 if x==0 else pm[x-1]
     return(f(mon[0]),tuple(sorted(f(j) for j in mon[1])))
    seen=set();ambient=[]
    for mon in ambientmons:
     if mon in seen:continue
     orb=set(apm(mon,pm) for pm in perms);seen.update(orb);ambient.append(sum(poly(t) for t in orb))
    # Expand all family vectors in the common direct sum of local monomials.
    rowkeys=[(B,mon) for B in subsets for mon in mons]
    expof={mon:next(iter(poly(mon).dict())) for mon in mons}
    F=matrix(ZZ,[[fam[B].monomial_coefficient(poly(mon)) for fam in families] for B,mon in rowkeys])
    Araw=matrix(ZZ,[[R(a.subs({u[0]:sum(u[i] for i in B),v[0]:sum(v[i] for i in B)})).monomial_coefficient(poly(mon)) for a in ambient] for B,mon in rowkeys])
    A=matrix(ZZ,F.solve_right(Araw));assert C*A==0
    K=C.right_kernel_matrix();coords=K.transpose().solve_right(A)


    ideal_rows={}
    for B in subsets:
     if len(B)==k:continue
     for j,a in enumerate(ambient):
      f=R(a.subs({u[0]:sum(u[i] for i in B),v[0]:sum(v[i] for i in B)}))
      for mon in mons:
       if mon[0] in B and all(i in B for i in mon[1]):continue
       c=f.monomial_coefficient(poly(mon))
       if c:ideal_rows.setdefault((B,mon),[0]*len(ambient))[j]=c
    I=matrix(ZZ,list(ideal_rows.values()))
    cd = [abs(int(x)) for x in C.elementary_divisors() if x]
    ad = [abs(int(x)) for x in A.elementary_divisors() if x]
    idg = [abs(int(x)) for x in I.elementary_divisors() if x]
    return dict(k=k, offset=r, proper_subsets_only=proper, family_dimension=len(families),
                ambient_dimension=len(ambient), compatibility_shape=list(map(int,C.dimensions())),
                compatibility_rank=len(cd), compatibility_nonunits=[x for x in cd if x>1],
                restriction_rank=len(ad), restriction_nonunits=[x for x in ad if x>1],
                gluing_obstruction_dimension=len(families)-len(cd)-len(ad),
                coefficient_rank=len(idg), coefficient_nonunits=[x for x in idg if x>1],
                coefficient_kernel_dimension=len(ambient)-len(idg))
