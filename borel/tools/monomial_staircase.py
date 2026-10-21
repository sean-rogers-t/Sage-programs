"""Count standard monomials using only integer tuples and recursion.

At a fixed exponent of one variable, divisibility by a monomial generator
reduces to divisibility in the remaining variables. The active generators
change only at their exponent levels. Summing these constant slices gives
the colength, including multiplicities, without a computer algebra system.
"""
from functools import lru_cache
from itertools import groupby


def staircase_length(generators):
    """Return the finite colength of the given monomial ideal.

    Each generator is an exponent tuple of the same length. A pure power
    of every variable must belong to the ideal; otherwise raise ValueError.
    """
    generators=tuple(tuple(e) for e in generators)
    if not generators or not generators[0]:
        raise ValueError('provide monomial generators in at least one variable')
    n=len(generators[0])
    if any(len(e)!=n or any(a<0 or int(a)!=a for a in e) for e in generators):
        raise ValueError('use nonnegative integer exponent tuples of equal length')
    generators=tuple(tuple(int(a) for a in e) for e in generators)

    @lru_cache(None)
    def minimal(gens):
        answer=[]
        ordered=sorted(set(gens),key=lambda e:(sum(e),e))
        for _,group in groupby(ordered,key=sum):
            # Distinct monomials of equal total degree cannot divide one
            # another. Compare a whole degree only with lower degrees.
            survivors=[e for e in group
                       if not any(all(a<=b for a,b in zip(f,e)) for f in answer)]
            answer.extend(survivors)
        return tuple(sorted(answer))

    @lru_cache(None)
    def count(gens):
        if any(not any(e) for e in gens):
            return 0
        n=len(gens[0])
        bounds=[min((e[j] for e in gens if e[j] and sum(e)==e[j]),default=None)
                for j in range(n)]
        if any(b is None for b in bounds):
            raise ValueError('the monomial ideal does not have finite colength')
        if n==1:
            return bounds[0]
        choices=[]
        for j,bound in enumerate(bounds):
            levels=sorted({0,bound}|{e[j] for e in gens if e[j]<bound})
            choices.append((len(levels),j,levels))
        _,j,levels=min(choices)
        answer=0
        for a,b in zip(levels,levels[1:]):
            sliced=minimal(tuple(e[:j]+e[j+1:] for e in gens if e[j]<=a))
            answer+=(b-a)*count(sliced)
        return answer

    return count(minimal(generators))
