"""Complete image presentations, runnable with SageMath.

Run from borel/: PYTHONPATH=src sage -python examples/kernel_presentations.py
Add --output PATH to save exact source/basis/certificate data as JSON.
"""

import argparse
import json
from pathlib import Path
from time import perf_counter

from sage_borel import BorelMap, GeneratorAlgebra, two_step_symplectic


def examples():
    projective = BorelMap("A2", omitted_nodes=1)
    yield "projective plane", GeneratorAlgebra(projective, {"h": projective.weights[0]})

    grassmannian = BorelMap("A3", omitted_nodes=2)
    w1, w2, _ = grassmannian.weights
    yield "Gr(2,4)", GeneratorAlgebra(
        grassmannian, {"c1": w2, "c2": w1 * (w2 - w1)})

    yield "C3 two-step joint image", two_step_symplectic(3, 1, 1).joint_algebra()

    orthogonal = BorelMap("D4", omitted_nodes=2)
    _, w2, w3, w4 = orthogonal.weights
    x3, x4 = w3 + w4 - w2, w4 - w3
    yield "D4 tail including Euler", GeneratorAlgebra(
        orthogonal, {"q": x3**2 + x4**2, "euler": x3*x4})

    exceptional = BorelMap("G2")
    yield "G2 full flag", GeneratorAlgebra(
        exceptional, dict(zip(("u", "v"), exceptional.weights)))

    exceptional = BorelMap("F4", omitted_nodes=4)
    yield "F4/P4 divisor image", GeneratorAlgebra(
        exceptional, {"h": exceptional.weights[3]})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    reports = []
    for name, algebra in examples():
        started = perf_counter()
        presentation = algebra.presentation()
        assert presentation.verify().certified
        print("{}: dimension {}, kernel {}, Hilbert {}".format(
            name, presentation.dimension, presentation.relations,
            presentation.hilbert_series))
        reports.append({"case": name, "seconds": perf_counter() - started,
                        "presentation": presentation.as_dict()})
    if arguments.output:
        arguments.output.write_text(json.dumps(reports, indent=2) + "\n")


if __name__ == "__main__":
    main()
