"""Independent projective, divisor, BGG and original-matrix checks."""
import unittest
from itertools import combinations
from math import prod as integer_product
from sage.all import ZZ, QQ, PolynomialRing, matrix, prod, gcd
from sage_borel import (BorelMap, type_a_equivariant, symplectic_projective_k_image)
from sage_borel.equivariant import _word


def original_projective_matrix(values):
    """Original 2n inclusion, with no even/odd reduction."""
    ring = PolynomialRing(ZZ,"L")
    L = ring.gen()
    P = prod(L*L-s*L+1 for s in values)
    n = len(values)
    inverse = -sum(P[k]*L**(k-1) for k in range(1,2*n+1))
    b = (L-inverse)%P
    columns = []
    power = ring.one()
    for j in range(2*n):
        columns.append([power[i] for i in range(2*n)])
        power = power*b%P
    return matrix(ZZ,columns).transpose()


class TypeAChecks(unittest.TestCase):
    def test_projective_restrictions(self):
        for m in range(2,7):
            for theory in ("H","K"):
                a = type_a_equivariant(m,1,theory=theory)
                M = a.restriction_matrix()
                for row,w in enumerate(a.basis):
                    p = a.diagonal_parameters
                    value = p[w[0]-1]
                    for j in range(m):
                        expected = prod(value-p[k] for k in range(j)) if theory=="H" else prod(1-p[k]/value for k in range(j))
                        self.assertEqual(M[row,j],expected)

    def test_full_flag_K_divisors(self):
        for m in (3,4):
            a = type_a_equivariant(m,theory="K")
            M = a.restriction_matrix()
            p = a.diagonal_parameters
            for node in range(1,m):
                v = list(range(1,m+1))
                v[node-1],v[node] = v[node],v[node-1]
                column = a.basis.index(tuple(v))
                for row,w in enumerate(a.basis):
                    self.assertEqual(M[row,column],1-prod(p[:node])/prod(p[w[i]-1] for i in range(node)))

    def test_H_specialization_matches_BGG(self):
        for m,omitted in ((3,None),(4,None),(4,(2,)),(5,(1,3))):
            a = type_a_equivariant(m,omitted,theory="H")
            old = BorelMap("A"+str(m-1),omitted_nodes=omitted)
            omega = old.weights
            eps = (omega[0],)+tuple(omega[i]-omega[i-1] for i in range(1,m-1))+(-omega[-1],)
            f = m*sum(a.weights[:a.block_sizes[0]])
            for polynomial in (f,f**2,f**3):
                actual = a.expand(polynomial).specialize()
                expected = old.expand(polynomial(*eps[:-1]))
                for w in a.basis:
                    element = old.W.from_reduced_word(_word(w))
                    self.assertEqual(actual.get(w,0),expected.coefficient(element))

    def test_membership_and_operations(self):
        a = type_a_equivariant(4,1,theory="K")
        x = a.weights[0]
        with self.assertRaises(ValueError):
            a.expand(x)
        generators = a.source_generators()
        self.assertTrue(all(a.is_integral_source(f) for f in generators))
        f,g = generators[:2]
        self.assertEqual(a.expand(f*g),a.expand(f)*a.expand(g))
        self.assertEqual(a.expand(f**3),a.expand(f)**3)
        for f in generators:
            self.assertEqual(a.expand(f).restrictions(),a.restrictions(f))
        rational = type_a_equivariant(4,1,theory="K",coefficients=QQ)
        rational.expand(rational.weights[0])
        partial = type_a_equivariant(5,(1,3))
        with self.assertRaises(NotImplementedError):
            partial.source_generators()
        with self.assertRaises(ValueError):
            partial.expand(partial.weights[1])
        self.assertEqual(type_a_equivariant(10,1).schubert_rank,10)
        with self.assertRaises(ValueError):
            type_a_equivariant(10)
        with self.assertRaises(ValueError):
            type_a_equivariant(3,True)

    def test_integral_H_projective_example(self):
        a = type_a_equivariant(4,1,theory="H")
        h = a.weights[0]
        u = a.diagonal_parameters
        e2 = sum(u[i]*u[j] for i in range(4) for j in range(i+1,4))
        E1,E2,E3 = a.source_generators()
        for w,value in a.restrictions(E2).items():
            self.assertEqual(value,6*u[w[0]-1]**2+e2)
        self.assertEqual(a.expand(E1).specialize(),{(2,1,3,4):4})

    def test_K_nonreduced_subword_sum(self):
        a = type_a_equivariant(4,theory="K")
        w = (4,3,2,1)
        word = _word(w)
        identity = (1,2,3,4)
        prefix = identity
        factors = []
        for node in word:
            p = a.diagonal_parameters
            factors.append(1-p[prefix[node-1]-1]/p[prefix[node]-1])
            v = list(prefix)
            v[node-1],v[node] = v[node],v[node-1]
            prefix = tuple(v)
        expected = {v:a.coefficient_ring.zero() for v in a.basis}
        for mask in range(1<<len(word)):
            v = list(identity)
            weight = a.coefficient_ring.one()
            taken = length = 0
            for j,node in enumerate(word):
                if mask & (1<<j):
                    taken += 1
                    weight *= factors[j]
                    if v[node-1]<v[node]:
                        v[node-1],v[node] = v[node],v[node-1]
                        length += 1
            expected[tuple(v)] += (-1)**(taken-length)*weight
        row = a.basis.index(w)
        self.assertEqual(tuple(a.restriction_matrix()[row]),tuple(expected[v] for v in a.basis))
        bad = {v:0 for v in a.basis}
        bad[w] = 1
        with self.assertRaises(ValueError):
            a.from_restrictions(bad)


class SymplecticChecks(unittest.TestCase):
    def test_original_matrix_smith(self):
        for n in range(1,7):
            c = symplectic_projective_k_image(n,spectral=True)
            for values in ([2]*n,list(range(1,n+1)),[(-1)**i*(i+1) for i in range(n)]):
                reduced = matrix(ZZ,n,lambda i,j:c.presentation_matrix()[i,j](*values))
                original = original_projective_matrix(values)
                expected = (1,)*n+tuple(abs(int(s)) for s in reduced.elementary_divisors())
                self.assertEqual(tuple(abs(int(s)) for s in original.elementary_divisors()),expected)
            self.assertEqual(c.ordinary_smith_factors(),tuple(abs(int(s)) for s in original_projective_matrix([2]*n).elementary_divisors()))

    def test_symbolic_determinants_and_fitting(self):
        for n in range(1,4):
            c = symplectic_projective_k_image(n,spectral=True)
            self.assertEqual(c.presentation_matrix().det(),c.determinant_generator())
            self.assertEqual(c.fitting_generators(n),(c.coefficient_ring.one(),))
            factors = c.cyclic_cokernel_factors()
            self.assertEqual(prod(factors),c.determinant_generator())
        c = symplectic_projective_k_image(3,spectral=True)
        self.assertEqual(c.fitting_generators(2),(2,))
        self.assertTrue(all(cyclic in c.fitting_ideal(1) for cyclic in c.fitting_generators(1)))
        c = symplectic_projective_k_image(4,spectral=True)
        s = c.s
        self.assertEqual(c.fitting_generators(3),(2*sum(s),2*sum(prod(x) for x in combinations(s,3))))
        with self.assertRaises(NotImplementedError):
            c.cyclic_cokernel_factors()
        with self.assertRaises(ValueError):
            symplectic_projective_k_image(8,spectral=True).fitting_generators(2,max_minors=1)

    def test_doubled_blocks_against_original_inclusion(self):
        for n in range(1,9):
            c = symplectic_projective_k_image(n,spectral=True)
            self.assertEqual(c.half_presentation_matrix().nrows(),n//2)
            self.assertIsNone(c._matrix)
            block = c.block_presentation_matrix()
            self.assertEqual(block.nrows(),n)
            self.assertTrue(block.is_immutable())
            # Includes singular fibers where the resultant vanishes.
            for values in ([2]*n,list(range(1,n+1)),[(-1)**i*(i+1) for i in range(n)],
                           [0]*n,[1,-1]+[2]*(n-2) if n>=2 else [0]):
                reduced = matrix(ZZ,n,lambda i,j:block[i,j](*values))
                actual = tuple(abs(int(s)) for s in original_projective_matrix(values).elementary_divisors())
                expected = (1,)*n+tuple(abs(int(s)) for s in reduced.elementary_divisors())
                self.assertEqual(actual,expected)

    def test_rank_four_closed_fitting(self):
        c = symplectic_projective_k_image(4,spectral=True)
        a = sum(c.s)
        d = sum(prod(x) for x in combinations(c.s,3))
        e2 = sum(prod(x) for x in combinations(c.s,2))
        e4 = prod(c.s)
        F = prod(c.s[i]+c.s[j] for i in range(4) for j in range(i+1,4))
        self.assertEqual(c.half_presentation_matrix(),
                         matrix(c.coefficient_ring,[[-2*d,2*a*e4],[-2*a,2*(a*e2-d)]]))
        self.assertEqual(d*d-a*e2*d+a*a*e4,-F)
        self.assertEqual(c.fitting_generators(1,max_minors=1),(8*F*a,8*F*d))
        self.assertEqual(c.fitting_generators(2,max_minors=1),(4*a*a,4*a*d,4*d*d))
        self.assertIsNone(c._matrix)
        self.assertEqual(c.ordinary_smith_factors(),(1,1,1,1,16,16,1024,1024))

    def test_half_minor_convolution(self):
        for n in (5,6):
            c = symplectic_projective_k_image(n,spectral=True)
            ideals = {i:c.fitting_generators(i) for i in range(1,n)}
            points = ([2]*n,list(range(1,n+1)),[1,-1]+[2]*(n-2)) if n==5 else ([2]*n,)
            for values in points:
                factors = tuple(abs(int(s)) for s in original_projective_matrix(values).elementary_divisors())[n:]
                for i in range(1,n):
                    generators = ideals[i]
                    # Sparse integer evaluation avoids Sage's slow generic
                    # callable substitution on expanded six-variable products.
                    evaluated = [sum(int(coefficient)*integer_product(v**power for v,power in zip(values,exponent))
                                     for exponent,coefficient in f.dict().items()) for f in generators]
                    self.assertEqual(abs(gcd(evaluated)),prod(factors[:n-i]))
            self.assertIsNone(c._matrix)

    def test_initial_column_half_minors(self):
        c = symplectic_projective_k_image(6,spectral=True)
        for j,number in ((0,1),(1,3),(2,3),(3,1)):
            # The full minor counts are 1,9,9,1. The compressed
            # routine must succeed with the strictly smaller guard.
            generators = c.half_minor_generators(j,max_minors=number)
            self.assertLessEqual(len(generators),number)
            for values in ([2]*6,[1,-1]+[2]*4):
                H = c.half_presentation_matrix()
                numeric = matrix(ZZ,3,lambda a,b:H[a,b](*values))
                expected = prod(numeric.elementary_divisors()[:j])
                evaluated = [sum(int(coefficient)*integer_product(v**power for v,power in zip(values,exponent))
                                 for exponent,coefficient in f.dict().items()) for f in generators]
                self.assertEqual(abs(gcd(evaluated)),abs(expected))
        self.assertEqual(c.half_minor_generators(4),(0,))
        with self.assertRaises(ValueError):
            c.half_minor_generators(2,max_minors=2)
        self.assertIsNone(c._matrix)
        # Exercise the general even/odd convolution beyond the closed ranks.
        for n,i,guard in ((8,6,32),(9,7,8)):
            c = symplectic_projective_k_image(n,spectral=True)
            generators = c.fitting_generators(i,max_minors=guard)
            evaluated = [sum(int(coefficient)*2**sum(exponent) for exponent,coefficient in f.dict().items())
                         for f in generators]
            factors = original_projective_matrix([2]*n).elementary_divisors()[n:]
            self.assertEqual(abs(gcd(evaluated)),abs(prod(factors[:n-i])))
            self.assertIsNone(c._matrix)

    def test_rank_five_and_cubic_closed_formulas(self):
        c = symplectic_projective_k_image(5,spectral=True)
        e = [c.coefficient_ring.one()]+[sum(prod(x) for x in combinations(c.s,j)) for j in range(1,6)]
        A,B = e[1]*e[4]-e[5],e[1]*e[2]-e[3]
        self.assertEqual(c.fitting_generators(2,max_minors=1),(8*A*A,8*A*B,8*B*B))
        self.assertEqual(c.fitting_generators(3,max_minors=1),(4*A,4*B))
        self.assertEqual(c.ordinary_smith_factors(),(1,1,1,1,1,2,128,128,32768,32768))
        # Rank six is independently tested in every degree above. Check
        # both cubic remainder-coefficient and square-ideal cases in rank seven.
        c = symplectic_projective_k_image(7,spectral=True)
        factors = original_projective_matrix([2]*7).elementary_divisors()[7:]
        for i in (4,5):
            generators = c.fitting_generators(i,max_minors=1)
            evaluated = [sum(int(coefficient)*2**sum(exponent) for exponent,coefficient in f.dict().items())
                         for f in generators]
            self.assertEqual(abs(gcd(evaluated)),abs(prod(factors[:7-i])))
        self.assertIsNone(c._matrix)

    def test_torus_base_change(self):
        for n in (1,2,3):
            spectral = symplectic_projective_k_image(n,spectral=True)
            torus = symplectic_projective_k_image(n)
            M = spectral.presentation_matrix()
            converted = matrix(torus.coefficient_ring,n,lambda i,j:M[i,j](*torus.s))
            self.assertEqual(converted,torus.presentation_matrix())
            r = n//2
            H = spectral.half_presentation_matrix()
            self.assertEqual(matrix(torus.coefficient_ring,r,r,lambda i,j:H[i,j](*torus.s)),
                             torus.half_presentation_matrix())
            self.assertEqual(spectral.ordinary_smith_factors(),torus.ordinary_smith_factors())

    def test_rank_nine_augmentation_without_symbolic_matrix(self):
        c = symplectic_projective_k_image(9,spectral=True)
        self.assertIsNone(c._matrix)
        factors = c.ordinary_smith_factors()
        self.assertIsNone(c._matrix)
        self.assertEqual(factors,tuple(abs(int(s)) for s in original_projective_matrix([2]*9).elementary_divisors()))


if __name__=="__main__":
    unittest.main()
