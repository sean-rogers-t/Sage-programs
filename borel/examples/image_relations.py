"""Run from borel/: PYTHONPATH=src sage -python examples/image_relations.py."""

from sage_borel import BorelMap, GeneratorAlgebra, two_step_symplectic


def main():
    projective = BorelMap("A2", omitted_nodes=1)
    A = GeneratorAlgebra(projective, {"h": projective.weights[0]})
    assert A.hilbert_function(4) == (1, 1, 1, 0, 0)
    assert A.first_relation_degree(4) == 3
    print("Projective-plane relation:", A.degree(3).relations)

    flags = two_step_symplectic(n=3, b=1, k=1)
    joint = flags.joint_algebra()
    c1, e1 = joint.domain.gens()
    report = joint.degree(4)
    assert (report.domain_dimension, report.image_dimension, report.kernel_dimension) == (3, 2, 1)
    assert not joint.evaluate(c1**4 + c1**2*e1 + e1**2)
    assert not joint.evaluate(e1**3)
    print("Two-step flag nodes:", flags.omitted_nodes)
    print("Degree-4 mixed relation:", report.relations)
    print("Complete small-case image Hilbert polynomial:", joint.hilbert_function(8))
    print("Presentation agrees:", joint.hilbert_function(8) == flags.predicted_hilbert_prefix(8))


if __name__ == "__main__":
    main()
