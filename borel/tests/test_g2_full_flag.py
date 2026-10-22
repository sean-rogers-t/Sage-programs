"""Full G2/B literal character image and its monic integral presentation."""
import unittest
from sage.all import ZZ, matrix
from sage_borel import g2_character_algebra


class G2FullFlagTests(unittest.TestCase):
    def test_complete_minimal_image_and_integer_kernel(self):
        algebra = g2_character_algebra()
        x,y = algebra.domain.gens()
        relations = (x*x-3*x*y+3*y*y,y**6)
        self.assertEqual(set(algebra.presentation().relations),set(relations))
        expected = ((1,),(2,6),(12,12),(24,144),(288,288),(576,1728),(3456,),())
        normal = [(0,j) for j in range(6)]+[(1,j) for j in range(6)]
        for d,factors in enumerate(expected):
            row = algebra.integral_degree(d,include_ambient=True)
            self.assertEqual(row.smith_factors,factors)
            self.assertTrue(row.full_ambient)
            selected = [p for p in normal if sum(p)==d]
            columns = [[algebra.evaluate(x**a*y**b).get(w,0) for a,b in selected]
                       for w in row.schubert_indices]
            self.assertEqual(matrix(ZZ,columns).rank(),len(selected))
        for f in relations:
            self.assertFalse(algebra.evaluate(f))

    def test_adjoint_scaling_and_original_maximal_parabolics(self):
        for node,top in ((None,6),(1,5),(2,5)):
            minimal = g2_character_algebra(node)
            adjoint = g2_character_algebra(node,representation='adjoint')
            for d in range(top+2):
                self.assertEqual(adjoint.integral_degree(d).smith_factors,
                    tuple(4**d*c for c in minimal.integral_degree(d).smith_factors))
        for bad in ('Minimal',7,None,True):
            with self.assertRaises(ValueError):
                g2_character_algebra(representation=bad)
        for bad in (0,3,True):
            with self.assertRaises(ValueError):
                g2_character_algebra(bad)


if __name__ == '__main__':
    unittest.main()
