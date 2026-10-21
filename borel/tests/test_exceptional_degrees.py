"""Degree certificates, independent volumes and protection of modified maps."""
from dataclasses import replace
import unittest

from sage.all import Polyhedron, matrix, factorial
from sage_borel import (exceptional_springer_degree, e6_adjoint_quotient,
                        e7_adjoint_quotient, g2_adjoint_quotient,
                        springer_character_map)


class ExceptionalDegreeTests(unittest.TestCase):
    def test_lattice_volume_and_center(self):
        for kind, name, weight in (
            ('E6','minimal',(1,0,0,0,0,0)),
            ('E6','adjoint',(0,1,0,0,0,0)),
            ('E7','minimal',(0,0,0,0,0,0,1)),
            ('E7','adjoint',(1,0,0,0,0,0,0))):
            with self.subTest(kind=kind, representation=name):
                c=exceptional_springer_degree(kind,representation=name)
                s=springer_character_map(kind,weight)
                C=matrix(s.cartan_type.cartan_matrix())
                vertices=[mu for mu,m in s.weight_multiplicities if any(mu)]
                volume=factorial(C.nrows())*Polyhedron(vertices=vertices).volume()
                if name=='adjoint':
                    volume/=C.det()
                self.assertEqual(volume,c.normalized_volume)
                self.assertEqual(c.central_cover,C.det() if name=='adjoint' else 1)

    def test_packaged_and_modified_coordinates(self):
        for helper in (e6_adjoint_quotient,e7_adjoint_quotient):
            for name in ('minimal','adjoint'):
                q=helper(representation=name)
                self.assertEqual(q.generic_degree(),q.degree_certificate.degree)
                with self.assertRaises(ValueError):
                    q.generic_degree(algorithm='elimination')
                changed=replace(q,character_polynomials=(
                    q.character_polynomials[0]+1,)+q.character_polynomials[1:])
                with self.assertRaisesRegex(ValueError,'not bound'):
                    changed.generic_degree()
        for kwargs in ({'algorithm':'numerical'},{'max_rank':0}):
            with self.assertRaises(ValueError):
                q.generic_degree(**kwargs)

    def test_uncached_elimination_and_validation(self):
        q=g2_adjoint_quotient()
        self.assertEqual(q.generic_degree(algorithm='elimination'),6)
        for kind,name in (('E8','minimal'),('E6','other')):
            with self.assertRaises(ValueError):
                exceptional_springer_degree(kind,representation=name)


if __name__=='__main__':
    unittest.main()
