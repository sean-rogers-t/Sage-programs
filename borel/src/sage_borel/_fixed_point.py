"""Exact fixed-point interpolation for a small maximal-parabolic quotient.

Billey's reduced-subword restriction formula is evaluated at a regular
rational point. Homogeneity extracts ordinary Schubert coefficients in the
matching degree. No full Weyl-group enumeration or polynomial products.
Private ordinary-image backend for the certified E6/P1 and E7/P7 factories.
Promoted from the independently checked exceptional localization experiment.
This does not expose an equivariant API.
"""
from collections import deque
from functools import lru_cache
from sage.all import QQ, WeylGroup, matrix, vector


class FixedPointBorel:
    def __init__(self, kind, node, point_weights=None):
        self.W = WeylGroup(kind, prefix='loc_s')
        self.rank = int(self.W.cartan_type().rank())
        self.node = node
        self.s = self.W.simple_reflections()
        domain = self.W.domain()
        self.fundamental = domain.fundamental_weights()
        roots = domain.simple_roots()
        weight = self.fundamental[node]
        representatives = {weight: self.W.one()}
        queue = deque([weight])
        while queue:
            current = queue.popleft()
            w = representatives[current]
            for i in range(1, self.rank+1):
                target = self.s[i].action(current)
                if target not in representatives:
                    representatives[target] = self.s[i]*w
                    queue.append(target)
                    if len(representatives) > 10000:
                        raise ValueError('fixed-point experiment limited to 10000 cosets')
        self.basis = tuple(sorted(representatives.values(),
                                 key=lambda w: (w.length(), tuple(w.reduced_word()))))
        self.lengths = tuple(int(w.length()) for w in self.basis)
        self.top = max(self.lengths)
        coefficients = point_weights or tuple(range(1, self.rank+1))
        point = sum(QQ(c)*self.fundamental[i]
                    for i,c in enumerate(coefficients,1))
        self.values = tuple(tuple(-QQ(w.action(self.fundamental[i]).dot_product(point))
                                  for i in range(1,self.rank+1)) for w in self.basis)
        state_ids, state_lengths, transitions = {}, [], []
        def state(v):
            if v in state_ids:
                return state_ids[v]
            index = len(state_ids)
            state_ids[v] = index
            state_lengths.append(int(v.length()))
            transitions.append({})
            for i in range(1,self.rank+1):
                if v.has_descent(i, side='right'):
                    transitions[index][i] = state(v*self.s[i])
            return index
        targets = [state(v) for v in self.basis]
        self.weak_states = len(state_ids)
        rows = []
        for w in self.basis:
            word = tuple(w.reduced_word())
            prefix = self.W.one()
            root_values = []
            for i in word:
                root_values.append(QQ(prefix.action(roots[i]).dot_product(point)))
                prefix *= self.s[i]
            assert all(v > 0 for v in root_values)
            @lru_cache(None)
            def restriction(index, length):
                if state_lengths[index] == 0:
                    return QQ.one()
                if length < state_lengths[index]:
                    return QQ.zero()
                answer = restriction(index, length-1)
                i = word[length-1]
                if i in transitions[index]:
                    answer += root_values[length-1]*restriction(transitions[index][i], length-1)
                return answer
            rows.append([restriction(index,len(word)) for index in targets])
            restriction.cache_clear()
        self.restrictions = matrix(QQ,rows)
        assert all(not self.restrictions[i,j] for i in range(len(rows))
                   for j in range(i+1,len(rows)))
        assert all(self.restrictions[i,i] for i in range(len(self.basis)))
        self.inverse = self.restrictions.inverse()

    def polynomial_values(self, polynomial):
        return vector(QQ,[polynomial(*point) for point in self.values])

    def coefficients(self, values, degree):
        equivariant = self.inverse*vector(QQ,values)
        assert all(not equivariant[i] for i,d in enumerate(self.lengths) if d > degree)
        return {w: equivariant[i] for i,w in enumerate(self.basis)
                if self.lengths[i] == degree and equivariant[i]}

    def expand(self, polynomial):
        return self.coefficients(self.polynomial_values(polynomial), int(polynomial.total_degree()))
