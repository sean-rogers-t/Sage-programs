"""Complete E7 quotient maps, full root products, and independent fixtures."""
import unittest
from sage.all import QQ, GF, matrix, vector, prod, CartanType, RootSystem
from sage_borel import (e7_adjoint_quotient, e7_cayley_determinant,
                       e7_fundamental_characters, SpringerGeometry,
                       springer_character_map)


class E7GeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases={}
        for name,weight in (('minimal',(0,0,0,0,0,0,1)),('adjoint',(1,0,0,0,0,0,0))):
            cls.cases[name]=(e7_adjoint_quotient(representation=name),
                             e7_cayley_determinant(representation=name),
                             springer_character_map('E7',weight))

    def test_all_degrees_and_quadratic_normalizations(self):
        expected={'minimal':[4,41,105,234,501,972,3101], 'adjoint':[5,78,223,535,1231,2515,8697]}
        for name,(q,det,s) in self.cases.items():
            self.assertEqual(q.degrees,(2,6,8,10,12,14,18))
            counts=[len(f.dict()) for f in q.character_polynomials]
            self.assertEqual(counts,expected[name])
            w=vector(s.ring.gens())
            self.assertEqual(q.invariants[0],w*matrix(QQ,CartanType('E7').cartan_matrix())*w/2)
            self.assertEqual(q.character_polynomials[0],
                SpringerGeometry(s).character_coordinates(s.pullback(q.invariants[0])))
            for degree,f in zip(q.degrees,q.character_polynomials):
                exponent=tuple(degree if i==(6 if name=='minimal' else 0) else 0 for i in range(7))
                numerator=2*3**degree+54 if name=='minimal' else 24
                denominator=(24 if name=='minimal' else 36)**degree*(24 if degree==2 else 1)
                monomial=prod(x**a for x,a in zip(f.parent().gens(),exponent))
                self.assertEqual(f.monomial_coefficient(monomial),QQ(numerator)/denominator)

    def test_every_coordinate_at_fresh_exact_points(self):
        for field in (QQ,GF(65537)):
            point=tuple(map(field,[-2,3,5,-7,11,13,17]))
            chars=e7_fundamental_characters(point)
            for q,det,s in self.cases.values():
                for f,pullback in zip(q.character_polynomials,q.torus_pullbacks):
                    polynomial=f if field is QQ else f.change_ring(field)
                    self.assertEqual(polynomial(*chars),pullback.evaluate(point))

    def test_basic_invariant_jacobian_certificate(self):
        q=self.cases['minimal'][0]
        C=matrix(QQ,CartanType('E7').cartan_matrix())
        point=C.inverse()*vector(QQ,[2,3,5,7,11,13,17])
        J=matrix(QQ,[[f.derivative(w)(*point) for w in f.parent().gens()] for f in q.invariants])
        roots=[C*vector(a.to_vector()) for a in RootSystem('E7').root_space().positive_roots()]
        discriminant=prod(sum(a*b for a,b in zip(root,point)) for root in roots)
        self.assertEqual(J.det()/discriminant,2676292792892522496000)

    def test_all_root_factors_have_explicit_weight_string_forms(self):
        C=matrix(QQ,CartanType('E7').cartan_matrix())
        metric=C.inverse()
        root_weights=[tuple(C*vector(a.to_vector())) for a in RootSystem('E7').root_space().roots()]
        weights=self.cases['minimal'][2].weight_multiplicities
        for name,(q,det,s) in self.cases.items():
            self.assertEqual(len(det.root_factors),63)
            monomial=lambda w:prod(t**a for t,a in zip(s.torus_ring.gens(),w))
            for alpha,r in det.root_factors:
                inner=lambda w:vector(QQ,w)*metric*vector(QQ,alpha)
                if name=='minimal':
                    selected=[w for w,m in weights if inner(w)==-1]
                    self.assertEqual(len(selected),12)
                    expected=sum(monomial(w) for w in selected)/12
                    self.assertEqual(len(r.dict()),12)
                else:
                    selected=[w for w in root_weights if inner(w)==-1]
                    self.assertEqual(len(selected),32)
                    expected=(2+2*monomial(tuple(-a for a in alpha))+
                              sum(monomial(w) for w in selected))/36
                    self.assertEqual(len(r.dict()),34)
                self.assertEqual(r,expected)

    def test_full_determinants_and_cartan_character_polynomials(self):
        for name,(q,det,s) in self.cases.items():
            self.assertEqual(det.evaluate([1]*7),1)
            central=tuple(-1 if i in (1,4,6) else 1 for i in range(7))
            self.assertEqual(det.evaluate(central),-1 if name=='minimal' else 1)
            self.assertEqual(len(det.torus_character_polynomial.dict()),24 if name=='minimal' else 68)
            for field in (QQ,GF(65537)):
                point=tuple(map(field,[2,3,5,7,11,13,17]))
                chars=e7_fundamental_characters(point)
                polynomial=det.torus_character_polynomial
                if field is not QQ:
                    polynomial=polynomial.change_ring(field)
                evaluate=lambda f:sum(field(c)*prod(t**a for t,a in zip(point,e)) for e,c in f.dict().items())
                direct=matrix(field,[[evaluate(f) for f in row] for row in det.torus_jacobian_matrix]).det()
                self.assertEqual(polynomial(*chars),direct)
                if field is not QQ:
                    # Chain rule through T/W and t/W, using the Weyl
                    # denominator for the fundamental-character Jacobian.
                    DF=matrix(field,[[f.change_ring(field)(*chars) for f in row]
                                     for row in q.character_jacobian()])
                    c=field(2676292792892522496000)
                    self.assertEqual(det.evaluate(point)*direct*c*c,DF.det()**2)
                    fixture=10310 if name=='minimal' else 51977
                    self.assertEqual(det.evaluate(point),field(fixture))

    def test_lazy_expansion_and_central_parity(self):
        odd=(1,4,6)
        for name,(q,det,s) in self.cases.items():
            self.assertEqual(q.torus_pullbacks[0].expand(),s.pullback(q.invariants[0]))
            with self.assertRaises(ValueError):
                q.torus_pullbacks[-1].expand(max_terms=1)
            for f in q.character_polynomials:
                self.assertTrue(all(sum(e[i] for i in odd)%2==0 for e in f.dict()))
            parity=1 if name=='minimal' else 0
            self.assertTrue(all(sum(e[i] for i in odd)%2==parity
                                for e in det.torus_character_polynomial.dict()))

    def test_validation_guards_and_exact_integer_points(self):
        for helper in (e7_adjoint_quotient,e7_cayley_determinant):
            with self.assertRaises(ValueError):
                helper(representation='other')
        self.assertEqual(tuple(map(int,e7_fundamental_characters([1]*7))),
                         (133,912,8645,365750,27664,1539,56))
        dimensions=e7_fundamental_characters([1]*7)
        for q,det,s in self.cases.values():
            self.assertTrue(all(f(*dimensions)==0 for f in q.character_polynomials))
        for point in ([1]*6,[0]*7):
            with self.assertRaises(ValueError):
                e7_fundamental_characters(point)
            for q,det,s in self.cases.values():
                with self.assertRaises(ValueError):
                    det.evaluate(point)
                with self.assertRaises(ValueError):
                    q.torus_pullbacks[0].evaluate(point)
        with self.assertRaises(ValueError):
            self.cases['minimal'][0].generic_degree()


if __name__=='__main__':
    unittest.main()
