"""Run with: sage -python reproduce_springer.py (tested with SageMath 9.5).

Bourbaki numbering; H_i = alpha_i^vee and t_i = omega_i(t).
The Springer coordinates satisfy theta(t) = sum_i Theta_i(t) H_i.
No published formula or previous computation is imported.
"""
from sage.all import (QQ, RootSystem, WeylCharacterRing, matrix, vector,
                      LaurentPolynomialRing, identity_matrix, prod)

cases = [('F4', (0, 0, 0, 1)),
         ('E6', (1, 0, 0, 0, 0, 0)),
         ('E7', (0, 0, 0, 0, 0, 0, 1))]

for kind, highest in cases:
    r = len(highest)
    ambient = RootSystem(kind).ambient_space()
    roots = ambient.simple_roots()
    coroots = ambient.simple_coroots()
    alpha = [tuple(roots[i+1].inner_product(coroots[j+1])
                   for j in range(r)) for i in range(r)]

    # Reflect a weight mu by s_i(mu) = mu - mu(H_i) alpha_i.
    orbit, pending = {highest}, [highest]
    while pending:
        mu = pending.pop()
        for i in range(r):
            nu = tuple(mu[j] - mu[i]*alpha[i][j] for j in range(r))
            if nu not in orbit:
                orbit.add(nu)
                pending.append(nu)
    multiplicities = {mu: 1 for mu in orbit}
    if kind == 'F4':
        multiplicities[(0,)*r] = 2

    # Cross-check the orbit construction against the irreducible character.
    character = WeylCharacterRing(kind, style='coroots')(*highest)
    check = {tuple(w.inner_product(coroots[i+1]) for i in range(r)): int(m)
             for w, m in character.weight_multiplicities().items()}
    assert check == multiplicities

    M = matrix(QQ, r, r, lambda i, j: sum(
        m*mu[i]*mu[j] for mu, m in multiplicities.items()))
    L = LaurentPolynomialRing(QQ, names=['t%d' % (i+1) for i in range(r)])
    t = L.gens()
    monomial = {mu: prod(t[i]**int(mu[i]) for i in range(r))
                for mu in multiplicities}
    b = vector(L, [sum(m*mu[j]*monomial[mu]
                      for mu, m in multiplicities.items()) for j in range(r)])
    Theta = M.inverse()*b
    D = matrix(QQ, r, r, lambda i, j:
               Theta[i].derivative(t[j])(*([1]*r)))
    assert D == identity_matrix(QQ, r)
    assert all(f(*([1]*r)) == 0 for f in Theta)

    print('\n%s: dimension %s, nonzero orbit size %s' %
          (kind, sum(multiplicities.values()), len(orbit)))
    print('M =\n%s\nM^-1 =\n%s' % (M, M.inverse()))
    for i, f in enumerate(Theta, 1):
        print('Theta_%d = %s' % (i, f))
    print('d(theta)_1 = I: verified exactly')
