"""Higher E6 coordinates, exact product determinants, and independent fixtures."""
import unittest
from sage.all import QQ, GF, matrix, vector, prod, CartanType
from sage_borel import (e6_adjoint_quotient, e6_cayley_determinant,
                       e6_fundamental_characters, SpringerGeometry,
                       springer_character_map)


class E6GeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases={}
        for name,weight in (('minimal',(1,0,0,0,0,0)),('adjoint',(0,1,0,0,0,0))):
            cls.cases[name]=(e6_adjoint_quotient(representation=name),
                             e6_cayley_determinant(representation=name),
                             springer_character_map('E6',weight))

    def test_all_degrees_and_quadratic_normalizations(self):
        expected={'minimal':[3,16,34,82,106,370], 'adjoint':[5,20,93,270,272,1459]}
        for name,(q,det,s) in self.cases.items():
            self.assertEqual(q.degrees,(2,5,6,8,9,12))
            self.assertEqual([len(f.dict()) for f in q.character_polynomials],expected[name])
            w=vector(s.ring.gens())
            self.assertEqual(q.invariants[0],w*matrix(QQ,CartanType('E6').cartan_matrix())*w/2)
            self.assertEqual(q.character_polynomials[0],
                SpringerGeometry(s).character_coordinates(s.pullback(q.invariants[0])))

    def test_every_coordinate_at_fresh_exact_points(self):
        for field in (QQ,GF(65537)):
            point=tuple(map(field,[-2,3,5,-7,11,13]))
            chars=e6_fundamental_characters(point)
            for q,det,s in self.cases.values():
                for f,pullback in zip(q.character_polynomials,q.torus_pullbacks):
                    polynomial=f if field is QQ else f.change_ring(field)
                    self.assertEqual(polynomial(*chars),pullback.evaluate(point))

    def test_basic_invariant_jacobian_certificate(self):
        q=self.cases['minimal'][0]
        C=matrix(QQ,CartanType('E6').cartan_matrix())
        point=C.inverse()*vector(QQ,[2,3,5,7,11,13])
        J=matrix(QQ,[[f.derivative(w)(*point) for w in f.parent().gens()] for f in q.invariants])
        from sage.all import RootSystem
        roots=[C*vector(a.to_vector()) for a in RootSystem('E6').root_space().positive_roots()]
        discriminant=prod(sum(a*b for a,b in zip(root,point)) for root in roots)
        self.assertEqual(J.det()/discriminant,423263232000)

    def test_all_root_factors_have_explicit_weight_string_forms(self):
        from sage.all import RootSystem
        C=matrix(QQ,CartanType('E6').cartan_matrix())
        metric=C.inverse()
        root_weights=[tuple(C*vector(a.to_vector())) for a in RootSystem('E6').root_space().roots()]
        weights=self.cases['minimal'][2].weight_multiplicities
        for name,(q,det,s) in self.cases.items():
            self.assertEqual(len(det.root_factors),36)
            monomial=lambda w:prod(t**a for t,a in zip(s.torus_ring.gens(),w))
            for alpha,r in det.root_factors:
                inner=lambda w:vector(QQ,w)*metric*vector(QQ,alpha)
                if name=='minimal':
                    selected=[w for w,m in weights if inner(w)==-1]
                    self.assertEqual(len(selected),6)
                    expected=sum(monomial(w) for w in selected)/6
                    self.assertEqual(len(r.dict()),6)
                else:
                    selected=[w for w in root_weights if inner(w)==-1]
                    self.assertEqual(len(selected),20)
                    expected=(2+2*monomial(tuple(-a for a in alpha))+
                              sum(monomial(w) for w in selected))/24
                    self.assertEqual(len(r.dict()),22)
                self.assertEqual(r,expected)

    def test_full_determinants_against_independent_matrix_fixtures(self):
        finite=GF(65537)
        point=tuple(map(finite,[2,3,5,7,11,13]))
        fixtures={'minimal':51226,'adjoint':39220}
        for name,(q,det,s) in self.cases.items():
            self.assertEqual(det.evaluate(point),fixtures[name])
            self.assertEqual(det.evaluate([1]*6),1)
            self.assertEqual(len(det.torus_character_polynomial.dict()),12 if name=='minimal' else 49)
            rational=tuple(map(QQ,[2,3,5,7,11,13]))
            chars=e6_fundamental_characters(rational)
            direct=matrix(QQ,[[f(*rational) for f in row] for row in det.torus_jacobian_matrix]).det()
            self.assertEqual(det.torus_character_polynomial(*chars),direct)

    def test_lazy_expansion_and_generic_factored_interface(self):
        for q,det,s in self.cases.values():
            self.assertEqual(q.torus_pullbacks[0].expand(),s.pullback(q.invariants[0]))
            with self.assertRaises(ValueError):
                q.torus_pullbacks[-1].expand(max_terms=1)
        s=springer_character_map('G2',(1,0))
        g=SpringerGeometry(s)
        self.assertEqual(g.factored_cayley_determinant().expand(),g.cayley_determinant())

    def test_validation_guards_and_exact_integer_point_evaluation(self):
        for helper in (e6_adjoint_quotient,e6_cayley_determinant):
            with self.assertRaises(ValueError):
                helper(representation='other')
        self.assertEqual(tuple(map(int,e6_fundamental_characters([1]*6))),
                         (27,78,351,2925,351,27))
        for point in ([1]*5,[0]*6):
            with self.assertRaises(ValueError):
                e6_fundamental_characters(point)
            for q,det,s in self.cases.values():
                with self.assertRaises(ValueError):
                    det.evaluate(point)
                with self.assertRaises(ValueError):
                    q.torus_pullbacks[0].evaluate(point)
        with self.assertRaises(ValueError):
            self.cases['minimal'][0].generic_degree()


if __name__=='__main__':
    unittest.main()
