"""Experiment-only sparse BGG frontier in orthonormal type-D coordinates.

Adjacent roots act by coordinate swaps. The fork root x_(n-1)+x_n
sends the final two coordinates to their negatives in reverse order.
Its divided difference is a finite alternating geometric sum.
No complete Weyl group or ambient Schubert basis is enumerated.
"""

from sage.all import ZZ

class CartesianDFrontier:
    def __init__(self, borel):
        if borel.cartan_type.type() != "D" or tuple(borel.nodes) != tuple(range(1, borel.rank + 1)):
            raise ValueError("standard-numbered type D is required")
        self.borel = borel
        self.rank = borel.rank
        self._ascents = {}
        self._minimal = {}
        self.calls = 0
        self.divided_differences = 0

    def _difference(self, coefficients, node):
        self.divided_differences += 1
        result = {}
        i = node - 1
        for exponent, coefficient in coefficients.items():
            if node == self.rank:
                i = self.rank-2
                first, second = exponent[i], exponent[i+1]
                if first > second:
                    for offset in range(first-second):
                        output = list(exponent)
                        output[i], output[i+1] = first-1-offset, second+offset
                        output = tuple(output)
                        result[output] = result.get(output, 0) + coefficient*((-1)**offset)
                elif first < second:
                    for offset in range(second-first):
                        output = list(exponent)
                        output[i], output[i+1] = first+offset, second-1-offset
                        output = tuple(output)
                        result[output] = result.get(output, 0) + coefficient*((-1)**offset)
                continue
            first, second = exponent[i], exponent[i + 1]
            if first > second:
                for offset in range(first - second):
                    output = list(exponent)
                    output[i], output[i + 1] = first - 1 - offset, second + offset
                    output = tuple(output)
                    result[output] = result.get(output, 0) + coefficient
            elif first < second:
                for offset in range(second - first):
                    output = list(exponent)
                    output[i], output[i + 1] = first + offset, second - 1 - offset
                    output = tuple(output)
                    result[output] = result.get(output, 0) - coefficient
        return {exponent: coefficient for exponent, coefficient in result.items() if coefficient}

    def expand(self, polynomial):
        self.calls += 1
        coefficients = {tuple(map(int, exponent)): int(ZZ(coefficient))
                        for exponent, coefficient in polynomial.dict().items() if coefficient}
        if not coefficients:
            return self.borel._expansion()
        degrees = {sum(exponent) for exponent in coefficients}
        if len(degrees) != 1:
            raise ValueError("this experiment accepts homogeneous inputs only")
        degree, = degrees
        for node in self.borel.retained_nodes:
            reflected = {}
            for exponent, coefficient in coefficients.items():
                image = list(exponent)
                if node == self.rank:
                    image[-2], image[-1] = image[-1], image[-2]
                    coefficient *= (-1)**(exponent[-2]+exponent[-1])
                else:
                    image[node - 1], image[node] = image[node], image[node - 1]
                reflected[tuple(image)] = coefficient
            if reflected != coefficients:
                raise ValueError("Cartesian input is not invariant under the retained Levi nodes")
        frontier = {self.borel.W.one(): coefficients}
        for depth in range(degree):
            following, visited = {}, set()
            for word, value in frontier.items():
                if word not in self._ascents:
                    self._ascents[word] = tuple(node for node in self.borel.nodes
                                                if not word.has_descent(node, side="left"))
                for node in self._ascents[word]:
                    target = self.borel.simple_reflections[node] * word
                    if target in visited:
                        continue
                    visited.add(target)
                    if target not in self._minimal:
                        self._minimal[target] = self.borel._is_minimal(target)
                    if not self._minimal[target]:
                        continue
                    difference = self._difference(value, node)
                    if difference:
                        following[target] = difference
            frontier = following
            if not frontier:
                return self.borel._expansion()
        zero = (0,) * self.rank
        if any(set(value) != {zero} for value in frontier.values()):
            raise ArithmeticError("degree-matching Cartesian frontier is not constant")
        return self.borel._expansion({word: value[zero] for word, value in frontier.items()})
