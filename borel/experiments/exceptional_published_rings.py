"""Independent ordinary-cohomology engine for E6/P1 and E7/P7.

Uses Duan--Zhao, math/0511332, Theorems 4 and 5, and positive-root
Chevalley multiplication. No fixed-point interpolation is used here.
The E6/P6 theorem is transported by the diagram involution to P1.
"""
from collections import defaultdict
from pathlib import Path
import sys
from sage.all import PolynomialRing, QQ, matrix, vector
from sage_borel import exceptional_borel

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from tests.test_exceptional import RootChevalleyOracle


class PublishedMinusculeRing:
    def __init__(self, kind):
        self.kind = kind
        self.node = 1 if kind == 'E6' else 7
        self.top = 16 if kind == 'E6' else 27
        self.borel = b = exceptional_borel(kind,self.node)
        self.oracle = oracle = RootChevalleyOracle(kind)
        self.R = PolynomialRing(QQ, ['h','s'] if kind == 'E6' else ['h','s','t'])
        if kind == 'E6':
            h,s = self.R.gens()
            relations = [2*h**9+3*h*s**2-6*h**5*s, s**3-6*h**4*s**2+h**12]
            word = (2,4,3,1)
            # A degree-four polynomial lift of the specified Schubert class.
            candidates = [b.weights[0]**4,b.generator_polynomial('e2')]
            labels = b.basis(4)
            expansions = [b.expand(p) for p in candidates]
            M = matrix(QQ,[[e.get(w,0) for e in expansions] for w in labels])
            coefficients = M.solve_right(vector(QQ,[int(w==b.W.from_reduced_word(word)) for w in labels]))
            lift = sum(c*p for c,p in zip(coefficients,candidates))
            squared = b.expand(lift**2)
            self.e6_square_anchor = squared.to_terms()
            anchors = [(0,self.R.one(),{oracle.W.one():QQ.one()}),
                       (4,s,{oracle.W.from_reduced_word(word):QQ.one()}),
                       (8,s*s,{oracle.W.from_reduced_word(w.reduced_word()):c for w,c in squared.items()})]
        elif kind == 'E7':
            h,s,t = self.R.gens()
            relations = [s*s-2*h*t,
                         2*s*t-9*h**4*s*s+6*h**9*s-h**14,
                         t*t+10*h**3*s**3-9*h**8*s*s+2*h**13*s]
            anchors = [(0,self.R.one(),{oracle.W.one():QQ.one()}),
                       (5,s,{oracle.W.from_reduced_word((2,4,5,6,7)):QQ.one()}),
                       (9,t,{oracle.W.from_reduced_word((1,5,4,2,3,4,5,6,7)):QQ.one()})]
        else:
            raise ValueError('published engine only covers E6/P1 and E7/P7')
        self.relations = tuple(relations)
        self.groebner = self.R.ideal(relations).groebner_basis()
        self.labels = {d: tuple(oracle.W.from_reduced_word(w.reduced_word()) for w in b.basis(d))
                       for d in range(self.top+1)}
        self.anchor_polynomials = {d: [] for d in self.labels}
        anchor_vectors = {d: [] for d in self.labels}
        for birth,p,expansion in anchors:
            for d in range(birth,self.top+1):
                assert set(expansion) <= set(self.labels[d])
                self.anchor_polynomials[d].append(h**(d-birth)*p)
                anchor_vectors[d].append(vector(QQ,[expansion.get(w,0) for w in self.labels[d]]))
                next_expansion = defaultdict(QQ)
                for w,c in expansion.items():
                    for target,weights in oracle.covers(w):
                        if weights[self.node-1]:
                            next_expansion[target] += c*weights[self.node-1]
                expansion = {w:c for w,c in next_expansion.items() if c}
        self.schubert_matrices = {}
        self.normal_matrices = {}
        self.normal_exponents = {}
        for d in self.labels:
            V = matrix(QQ,anchor_vectors[d]).transpose()
            normal = [p.reduce(self.groebner) for p in self.anchor_polynomials[d]]
            exponents = sorted({tuple(e) for p in normal for e in p.dict()})
            dictionaries = [{tuple(e):c for e,c in p.dict().items()} for p in normal]
            N = matrix(QQ,[[p.get(e,0) for p in dictionaries] for e in exponents])
            assert V.rank() == N.rank() == len(self.labels[d]), (kind,d)
            assert all(not V*v for v in N.right_kernel().basis())
            self.schubert_matrices[d] = V
            self.normal_matrices[d] = N
            self.normal_exponents[d] = exponents
        # Determine each input's expression using BGG generator images and
        # independently propagated ordinary Chevalley columns only.
        self.generator_expressions = {}
        for g in b.generators():
            d = g.polynomial_degree
            image = {oracle.W.from_reduced_word(w.reduced_word()):c
                     for w,c in b.generator_image(g.name).items()}
            v = vector(QQ,[image.get(w,0) for w in self.labels[d]])
            coefficients = self.schubert_matrices[d].solve_right(v)
            expression = self.R(sum(c*p for c,p in zip(coefficients,self.anchor_polynomials[d])))
            assert self.expand(expression,d) == image
            self.generator_expressions[g.name] = expression

    def expand(self, polynomial, degree):
        if degree > self.top:
            assert not self.R(polynomial).reduce(self.groebner)
            return {}
        normal = self.R(polynomial).reduce(self.groebner)
        coefficients = {tuple(e):c for e,c in normal.dict().items()}
        exponents = self.normal_exponents[degree]
        assert set(coefficients) <= set(exponents)
        coordinates = self.normal_matrices[degree].solve_right(vector(QQ,[coefficients.get(e,0) for e in exponents]))
        answer = self.schubert_matrices[degree]*coordinates
        return {w:c for w,c in zip(self.labels[degree],answer) if c}
