"""Complete dominance-bounded fundamental-character monomial sets for E6."""
from sage.all import CartanType, QQ, ZZ, matrix, vector, floor


def candidates(highest, degree):
    C = matrix(QQ, CartanType('E6').cartan_matrix())
    inv = C.inverse()
    target = inv*vector(highest)*degree
    columns = inv.columns()
    answer = []
    def walk(i, remaining, prefix):
        if i == 6:
            if all(a.denominator() == 1 for a in remaining):
                answer.append(tuple(prefix))
            return
        bound = min(ZZ(floor(a/b)) for a,b in zip(remaining, columns[i]))
        for k in range(int(bound)+1):
            walk(i+1, remaining-k*columns[i], prefix+[k])
    walk(0, target, [])
    return answer
