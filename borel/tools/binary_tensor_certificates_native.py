"""Independent fixed tensor certificates via coefficient dictionaries.

No Sage or polynomial-substitution certificate is imported.
"""
from itertools import permutations, combinations, combinations_with_replacement, product
from math import prod
import sympy as s
from sympy.matrices.normalforms import smith_normal_form


def cubic_offset_obstruction_native():
    """Independent dictionary expansion of the explicit three-coordinate family."""
    subsets = [B for n in range(4) for B in combinations(range(1,4), n)]

    def term(result, i, indices, value):
        key = i, tuple(sorted(indices))
        result[key] = result.get(key, 0)+value

    def transform(poly, substitutions):
        out = {}
        for (i,jj), value in poly.items():
            for chosen in product(*(substitutions.get(t,[(t,1)]) for t in (i,)+jj)):
                term(out, chosen[0][0], [t[0] for t in chosen[1:]],
                     value*prod(t[1] for t in chosen))
        return {key:value for key,value in out.items() if value}

    family = {B:{} for B in subsets}
    for B in subsets:
        if len(B) == 2:
            j,l = B
            i = next(i for i in range(1,4) if i not in B)
            for u,v,sign in [(j,l,1),(l,j,-1)]:
                for t,c in [(j,1),(l,1),(i,-1)]:
                    for h,e in [(j,1),(l,-1)]:
                        term(family[B],u,[v,t,h],sign*c*e)
        elif len(B) == 3:
            for i in range(1,4):
                for j in range(1,4):
                    term(family[B],i,[i,i,j],1)
                    term(family[B],i,[j,j,j],-1)
    for B,C in combinations(subsets,2):
        coeff = {i:int(i in B)-int(i in C) for i in range(1,4)}
        pivot = next(i for i in coeff if coeff[i])
        substitution = {pivot:[(i,-coeff[i]//coeff[pivot]) for i in coeff if i != pivot and coeff[i]]}
        assert transform(family[B],substitution) == transform(family[C],substitution)
    mons = [(i,jj) for i in range(4) for jj in combinations_with_replacement(range(4),3)]
    ambient, seen = [],set()
    for mon in mons:
        if mon in seen:
            continue
        orbit = set()
        for pm in permutations(range(1,4)):
            f = lambda x: 0 if x == 0 else pm[x-1]
            orbit.add((f(mon[0]),tuple(sorted(f(t) for t in mon[1]))))
        seen.update(orbit)
        ambient.append({key:1 for key in orbit})
    local = [(i,jj) for i in range(1,4) for jj in combinations_with_replacement(range(1,4),3)]
    records = []
    for proper in [False,True]:
        chosen = [B for B in subsets if not proper or len(B) < 3]
        restricted = {(B,i):transform(poly,{0:[(j,1) for j in B]})
                      for B in chosen for i,poly in enumerate(ambient)}
        A = s.Matrix([[restricted[B,i].get(mon,0) for i in range(20)]
                      for B in chosen for mon in local])
        W = s.Matrix([family[B].get(mon,0) for B in chosen for mon in local])
        smith = smith_normal_form(A.row_join(W),domain=s.ZZ)
        diagonal = [abs(int(smith[i,i])) for i in range(min(smith.shape)) if smith[i,i]]
        assert A.rank() == 20 and diagonal == [1]*21
        records.append(dict(proper_subsets_only=proper,ambient_rank=20,
                            augmented_rank=21,augmented_nonunits=[]))
    return dict(integer_pairwise_compatibility=True,records=records,
                obstruction_restricts_to_primitive_proper_obstruction=True)


def binary_tensor_certificate_native(k, r, proper=False):
    """Build tensor substitutions as integer coefficient dictionaries.

    No polynomial-substitution builder or Sage certificate is imported.
    Indices encode one u coordinate and an ordered multiset of r v coordinates.
    """
    perms = list(permutations(range(1,k+1)))
    subsets = [tuple(B) for n in range(k if proper else k+1) for B in combinations(range(1,k+1),n)]
    local = [(i,jj) for i in range(1,k+1) for jj in combinations_with_replacement(range(1,k+1),r)]
    ambient = [(i,jj) for i in range(k+1) for jj in combinations_with_replacement(range(k+1),r)]
    def permute(mon, pm):
        f = lambda x: 0 if x==0 else pm[x-1]
        return f(mon[0]),tuple(sorted(f(j) for j in mon[1]))
    def orbits(items, group):
        seen, result = set(), []
        for mon in items:
            if mon in seen:
                continue
            orbit = {permute(mon,pm) for pm in group}
            seen.update(orbit)
            result.append(orbit)
        return result
    def transform(polynomial, substitutions):
        result = {}
        for (i,jj), coefficient in polynomial.items():
            choices = [substitutions.get(t,[(t,1)]) for t in (i,)+jj]
            for selected in product(*choices):
                mon = selected[0][0],tuple(sorted(t[0] for t in selected[1:]))
                result[mon] = result.get(mon,0)+coefficient*prod(t[1] for t in selected)
        return {mon:c for mon,c in result.items() if c}
    families, representatives = [], []
    for size in range(k if proper else k+1):
        B0 = tuple(range(1,size+1))
        group = [pm for pm in perms if set(pm[:size])==set(B0)]
        for orbit in orbits(local,group):
            family = {B:{} for B in subsets}
            for B in subsets:
                if len(B)==size:
                    pm = next(pm for pm in perms if set(pm[:size])==set(B))
                    family[B] = {permute(mon,pm):1 for mon in orbit}
            families.append(family)
            representatives.append((B0,next(mon for mon in local if mon in orbit)))
    constraints = {}
    for B,C in combinations(subsets,2):
        coefficients = {i:int(i in B)-int(i in C) for i in range(1,k+1)}
        pivot = next(i for i,c in coefficients.items() if c)
        substitution = {pivot:[(i,-c//coefficients[pivot]) for i,c in coefficients.items()
                               if i!=pivot and c]}
        for column,family in enumerate(families):
            difference = family[B].copy()
            for mon,c in family[C].items():
                difference[mon] = difference.get(mon,0)-c
            for mon,c in transform(difference,substitution).items():
                constraints.setdefault((B,C,mon),[0]*len(families))[column] += c
    global_polynomials = [{mon:1 for mon in orbit} for orbit in orbits(ambient,perms)]
    restricted = {(B,i):transform(poly,{0:[(j,1) for j in B]})
                  for B in subsets for i,poly in enumerate(global_polynomials)}
    restriction = s.Matrix([[restricted[B,i].get(mon,0) for i in range(len(global_polynomials))]
                            for B,mon in representatives])
    compatibility = s.Matrix(list(constraints.values()))
    ideal_constraints = {}
    for B in subsets:
        if len(B)==k:
            continue
        for i in range(len(global_polynomials)):
            for mon,c in restricted[B,i].items():
                if mon[0] not in B or any(j not in B for j in mon[1]):
                    ideal_constraints.setdefault((B,mon),[0]*len(global_polynomials))[i] = c
    ideal = s.Matrix(list(ideal_constraints.values()))
    assert compatibility*restriction == s.zeros(compatibility.rows,restriction.cols)
    def nonzero_smith(mat):
        diagonal = smith_normal_form(mat,domain=s.ZZ)
        return [abs(int(diagonal[i,i])) for i in range(min(diagonal.shape)) if diagonal[i,i]]
    cd,ad,idg = map(nonzero_smith,(compatibility,restriction,ideal))
    return dict(k=k, offset=r, proper_subsets_only=proper, family_dimension=len(families),
                ambient_dimension=len(global_polynomials), compatibility_shape=list(compatibility.shape),
                compatibility_rank=len(cd), compatibility_nonunits=[x for x in cd if x>1],
                restriction_rank=len(ad), restriction_nonunits=[x for x in ad if x>1],
                gluing_obstruction_dimension=len(families)-len(cd)-len(ad),
                coefficient_rank=len(idg), coefficient_nonunits=[x for x in idg if x>1],
                coefficient_kernel_dimension=len(global_polynomials)-len(idg))
