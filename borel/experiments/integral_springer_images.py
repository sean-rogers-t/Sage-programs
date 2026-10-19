"""Bounded exact integral Springer-image calculations; run inside SageMath.

From borel/: PYTHONPATH=src sage -python experiments/integral_springer_images.py
The inputs are fixed normalized catalogue/classical generators. This computes
degreewise Z-spans, not the intersection with the literal integral
representation ring. Finite cases do not establish a uniform theorem.
"""

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from sage.all import lcm
from sage.env import SAGE_VERSION
from sage_borel import Expansion, GeneratorAlgebra, classical_borel, exceptional_borel


def classical_case(family, n, k, bound):
    borel = classical_borel(family, n, k)
    generators = {"e%d" % j: borel.tail_polynomial(k, j) for j in range(1, k + 1)}
    if family == "D":
        generators["u"] = borel.euler_polynomial(k)
    return "%s(n=%d,k=%d) tail%s" % (family, n, k, " plus Euler" if family == "D" else ""), borel, generators, bound


def cases():
    for args in (("A", 5, 2, 6), ("C", 3, 1, 6), ("C", 4, 2, 6),
                 ("C", 5, 2, 6), ("B", 3, 1, 6), ("B", 4, 2, 6),
                 ("D", 4, 2, 6), ("D", 4, 3, 6), ("D", 5, 2, 6)):
        yield classical_case(*args)
    for kind, node, bound in (("G2", 1, 6), ("G2", 2, 6), ("F4", 1, 6),
                              ("F4", 4, 6), ("E6", 1, 6), ("E7", 7, 4)):
        borel = exceptional_borel(kind, node)
        generators = {info.name: borel.generator_polynomial(info.name)
                      for info in borel.generators()}
        yield "%s/P%d catalogue" % (kind, node), borel, generators, bound


def column_class(borel, report, position, saturated):
    basis = report.saturation_lattice_basis if saturated else report.image_lattice_basis
    return Expansion(borel, {w: basis[row, position]
                             for row, w in enumerate(report.schubert_indices)
                             if basis[row, position]})


def run_case(label, borel, generators, bound):
    started = perf_counter()
    algebra = GeneratorAlgebra(borel, generators)
    rows = []
    witness = None
    preimages_checked = 0
    for degree in range(bound + 1):
        result = algebra.integral_degree(degree, include_ambient=True)
        for position in range(result.rank):
            image = column_class(borel, result, position, False)
            assert result.contains(image)
            assert algebra.evaluate(result.preimage(image)) == image
            preimages_checked += 1
        if witness is None and not result.is_saturated:
            position = next(j for j, factor in enumerate(result.smith_factors) if factor > 1)
            factor = result.smith_factors[position]
            primitive = column_class(borel, result, position, True)
            assert not result.contains(primitive)
            assert result.contains(factor * primitive)
            preimage = result.preimage(factor * primitive)
            assert algebra.evaluate(preimage) == factor * primitive
            witness = {"polynomial_degree": degree, "primitive_class": primitive.to_terms(),
                       "least_positive_multiple_in_image": int(factor),
                       "integer_preimage_of_multiple": str(preimage)}
        rows.append({"polynomial_degree": degree, "cohomological_degree": 2 * degree,
                     "rank": result.rank, "ambient_dimension": result.ambient_dimension,
                     "free_cokernel_rank": result.free_cokernel_rank,
                     "smith_factors": list(map(int, result.smith_factors)),
                     "torsion_factors": list(map(int, result.torsion_factors)),
                     "saturation_index": int(result.saturation_index),
                     "is_saturated": result.is_saturated})
        print("%s d=%d: rank=%d/%d, Smith=%s" %
              (label, degree, result.rank, result.ambient_dimension, result.smith_factors),
              flush=True)
    return {"case": label, "cartan_type": str(borel.cartan_type),
            "omitted_nodes": list(borel.omitted_nodes),
            "generator_names": list(algebra.generator_names),
            "generator_polynomial_degrees": list(algebra.generator_degrees),
            "generator_polynomial_sha256": [
                sha256(str(p).encode("utf-8")).hexdigest() for p in algebra.generator_polynomials],
            "maximum_polynomial_degree": bound, "degrees": rows,
            "first_recorded_nonsaturation_witness": witness,
            "integer_preimages_checked": preimages_checked,
            "seconds": round(perf_counter() - started, 6)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("integral-springer-results.json"))
    args = parser.parse_args()
    report = {"sage_version": SAGE_VERSION,
              "generated_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "Bounded degreewise images of explicitly normalized generator Z-algebras. No uniform theorem or novelty claim follows from this finite experiment.",
              "command": "PYTHONPATH=src sage -python experiments/integral_springer_images.py",
              "rejected_original_normalizations": [], "cases": []}
    for label, borel, generators, bound in cases():
        images = {name: borel.expand(polynomial) for name, polynomial in generators.items()}
        scales = {name: int(lcm([coefficient.denominator() for coefficient in image.values()]))
                  for name, image in images.items()}
        original_label = label
        if any(scale != 1 for scale in scales.values()):
            rejected = GeneratorAlgebra(borel, generators)
            try:
                rejected.integral_degree(0)
            except ValueError as error:
                rejection = str(error)
            else:
                raise AssertionError("a nonintegral original normalization was accepted")
            report["rejected_original_normalizations"].append({
                "case": original_label, "reason": rejection,
                "nonintegral_generator_images": {name: image.to_terms()
                                                for name, image in images.items()
                                                if scales[name] != 1}})
            # This is an explicitly DIFFERENT generator Z-form. It is not
            # silently substituted for the original form or identified with
            # its intersection with the integral representation ring.
            generators = {name: scales[name] * polynomial
                          for name, polynomial in generators.items()}
            label += " EXPLICITLY RESCALED IMAGE-INTEGRAL FORM"
        result = run_case(label, borel, generators, bound)
        result["original_case"] = original_label
        result["explicit_generator_multipliers"] = scales
        result["original_normalization_used"] = all(scale == 1 for scale in scales.values())
        report["cases"].append(result)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        if original_label.startswith(("B(", "C(", "D(")):
            # In the defining/vector Springer setting, the literal integral
            # polynomial-character form has these explicitly different scales.
            # See INTEGRAL_LITERATURE.md for the character-intersection proof.
            literal_scales = {name: 2 ** int(polynomial.total_degree())
                              for name, polynomial in generators.items()}
            literal = run_case(original_label + " LITERAL POLYNOMIAL-CHARACTER FORM",
                               borel, {name: literal_scales[name] * polynomial
                                       for name, polynomial in generators.items()}, bound)
            for normalized, scaled in zip(result["degrees"], literal["degrees"]):
                d = normalized["polynomial_degree"]
                assert scaled["smith_factors"] == [2 ** d * factor
                                                   for factor in normalized["smith_factors"]]
            literal["original_case"] = original_label
            literal["explicit_generator_multipliers"] = literal_scales
            literal["original_normalization_used"] = False
            literal["integral_form"] = "Literal polynomial-character intersection for the defining/vector Springer map of the Sp/SO tail factor; not an arbitrary representation or full Levi assertion."
            report["cases"].append(literal)
            args.output.write_text(json.dumps(report, indent=2) + "\n")
    print("All %d integral Springer-image cases passed." % len(report["cases"]), flush=True)


if __name__ == "__main__":
    main()
