"""Closed lattice formulas compared with the full Laurent construction."""
import unittest
from sage.all import CartanMatrix, CartanType, QQ, identity_matrix, matrix
from sage_borel import (IrreducibleCharacterLattice, SpringerCharacterMap,
                       irreducible_character_lattice, springer_character_map)


class CharacterLatticeTests(unittest.TestCase):
    def same_lattice(self, A, B):
        transition = A.change_ring(QQ).inverse()*B
        self.assertEqual(transition.denominator(), 1)
        self.assertEqual(abs(transition.det()), 1)

    def test_exceptional_minimal_and_adjoint_normalizations(self):
        cases = [("G2",(1,0),7,2,"Q_long"), ("G2",(0,1),14,8,"Q_long"),
                 ("F4",(0,0,0,1),26,6,"Q_long"), ("F4",(1,0,0,0),52,18,"Q_long"),
                 ("E6",(1,0,0,0,0,0),27,6,"Q"), ("E6",(0,1,0,0,0,0),78,24,"P"),
                 ("E7",(0,0,0,0,0,0,1),56,12,"Q"), ("E7",(1,0,0,0,0,0,0),133,36,"P"),
                 ("E8",(0,0,0,0,0,0,0,1),248,60,"P")]
        for kind,highest,dimension,index,name in cases:
            with self.subTest(kind=kind, highest=highest):
                fast = irreducible_character_lattice(kind,highest)
                full = springer_character_map(kind,highest)
                self.assertIsInstance(fast,IrreducibleCharacterLattice)
                self.assertEqual((fast.dimension,fast.trace_index,fast.lattice_type),
                                 (dimension,index,name))
                self.assertEqual(fast.trace_gram,full.trace_gram)
                self.same_lattice(fast.linear_basis_matrix,full.linear_basis_matrix)
                self.assertTrue(all(full.is_integral(f) for f in fast.linear_generators))

    def test_composite_center_and_nonsimply_laced_controls(self):
        for kind,highest in [("A3",(0,1,0)), ("B3",(0,0,2)),
                             ("C3",(0,0,2)), ("D4",(0,0,1,0)), ("A1",(4,))]:
            fast = irreducible_character_lattice(kind,highest)
            full = springer_character_map(kind,highest)
            self.assertEqual(fast.trace_gram,full.trace_gram)
            self.same_lattice(fast.linear_basis_matrix,full.linear_basis_matrix)
        a3 = irreducible_character_lattice("A3",(0,1,0))
        self.assertEqual(a3.central_character_order,2)
        self.assertEqual(abs(a3.moment_lattice_basis.det()),2)

    def test_highest_weight_only_calculation_for_large_modules(self):
        fast = irreducible_character_lattice("E8",(100,0,0,0,0,0,0,0))
        self.assertGreater(fast.dimension,10**30)
        self.assertEqual(fast.central_character_order,1)
        self.same_lattice(fast.linear_basis_matrix,
                          fast.trace_index*identity_matrix(QQ,8))

    def test_irreducibility_is_essential(self):
        triple = SpringerCharacterMap("F4",[(0,0,0,1)]*3)
        adjoint = springer_character_map("F4",(1,0,0,0))
        self.assertEqual(triple.trace_gram,adjoint.trace_gram)
        self.same_lattice(3*triple.linear_basis_matrix,adjoint.linear_basis_matrix)
        self.assertEqual(abs((triple.linear_basis_matrix.inverse()
                             *adjoint.linear_basis_matrix).det()),3**4)

    def test_reject_invalid_irreducible_inputs(self):
        for weight in [(0,0),(-1,0),(1,),(True,0),(0.5,0),[(1,0),(0,1)]]:
            with self.assertRaises(ValueError):
                irreducible_character_lattice("G2",weight)
        with self.assertRaises(ValueError):
            irreducible_character_lattice(["A",2,1],(1,0,0))
        with self.assertRaises(ValueError):
            irreducible_character_lattice(CartanType("G2").relabel({1:2,2:1}),(1,0))


if __name__ == "__main__":
    unittest.main()
