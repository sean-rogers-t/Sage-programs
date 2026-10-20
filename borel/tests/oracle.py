"""Independent BGG oracle used by the regression suite.

This code uses explicit orthonormal coordinate reflections, rather than the
library's fundamental-weight reflection formulas or sparse search.  Expected
Schubert representatives are obtained from the normalized root discriminant.
Only Sage's finite Weyl group supplies word multiplication and length.
"""

from sage.all import PolynomialRing, QQ, WeylGroup, prod


class OrthonormalBGG:
    def __init__(self, kind):
        self.kind = kind
        self.family, self.rank = kind[0], int(kind[1:])
        self.W = WeylGroup(kind, prefix="oracle_s")
        self.w0 = self.W.long_element()
        size = self.rank + 1 if self.family == "A" else self.rank
        self.ring = PolynomialRing(QQ, ["x%d" % i for i in range(1, size + 1)])
        self.x = self.ring.gens()
        roots = [self.x[i] - self.x[j] for i in range(size)
                 for j in range(i + 1, size)]
        if self.family != "A":
            roots += [self.x[i] + self.x[j] for i in range(size)
                      for j in range(i + 1, size)]
        if self.family == "B":
            roots += list(self.x)
        if self.family == "C":
            roots += [2 * x for x in self.x]
        self.top = self.ring(prod(roots) / self.W.cardinality())
        self._schubert = {}
        self._levels = None

    def reflect(self, f, node):
        images = list(self.x)
        if node < self.rank or self.family == "A":
            images[node - 1], images[node] = images[node], images[node - 1]
        elif self.family == "D":
            images[-2], images[-1] = -self.x[-1], -self.x[-2]
        else:
            images[-1] = -self.x[-1]
        return self.ring(f)(*images)

    def divided(self, f, node):
        if node < self.rank or self.family == "A":
            alpha = self.x[node - 1] - self.x[node]
        elif self.family == "D":
            alpha = self.x[-2] + self.x[-1]
        else:
            alpha = (2 if self.family == "C" else 1) * self.x[-1]
        quotient, remainder = (self.ring(f) - self.reflect(f, node)).quo_rem(alpha)
        if remainder:
            raise ArithmeticError("Oracle encountered a nonexact divided difference")
        return quotient

    def operate(self, f, element):
        for node in reversed(element.reduced_word()):
            f = self.divided(f, node)
        return f

    def omega(self, node):
        if self.family == "A":
            return sum(self.x[:node]) - QQ(node) / len(self.x) * sum(self.x)
        if self.family == "B" and node == self.rank:
            return sum(self.x) / 2
        if self.family == "D" and node >= self.rank - 1:
            sign = 1 if node == self.rank else -1
            return (sum(self.x[:-1]) + sign * self.x[-1]) / 2
        return sum(self.x[:node])

    def from_fundamental(self, f):
        return self.ring(f(*[self.omega(i) for i in range(1, self.rank + 1)]))

    def levels(self):
        if self._levels is None:
            self._levels = {}
            for w in self.W:
                self._levels.setdefault(w.length(), []).append(w)
        return self._levels

    def schubert(self, w):
        if w not in self._schubert:
            self._schubert[w] = self.operate(self.top, w.inverse() * self.w0)
        return self._schubert[w]

    def expansion(self, f, degree=None):
        f = self.ring(f)
        degrees = [degree] if degree is not None else range(max(0, f.degree()) + 1)
        answer = {}
        zeros = [0] * len(self.x)
        for d in degrees:
            for w in self.levels().get(d, []):
                coefficient = QQ(self.operate(f, w)(*zeros))
                if coefficient:
                    answer[w] = coefficient
        return answer


def elementary(variables, degree):
    """Coefficient recursion independent of the library's implementation."""
    if degree < 0:
        return 0
    coefficients = [1] + [0] * degree
    for x in variables:
        for i in range(degree, 0, -1):
            coefficients[i] += x * coefficients[i - 1]
    return coefficients[degree]
