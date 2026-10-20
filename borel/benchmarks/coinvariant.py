"""Reproduce the formerly slow rank-nine complete-kernel calculation.

Run from borel/ with a fresh Sage process:
PYTHONPATH=src sage -python benchmarks/coinvariant.py --output benchmarks/coinvariant-results.json
"""

import argparse
import json
from pathlib import Path
from time import perf_counter

SCRIPT_STARTED = perf_counter()

from sage.env import SAGE_VERSION
from sage_borel import BorelMap, GeneratorAlgebra, __version__

IMPORT_SECONDS = perf_counter() - SCRIPT_STARTED


def run(rank):
    started = perf_counter()
    borel = BorelMap("C{}".format(rank), omitted_nodes=rank - 1)
    algebra = GeneratorAlgebra(borel, {"q": borel.tail_polynomial(1, 1)})
    q, = algebra.domain.gens()
    setup_finished = perf_counter()
    presentation = algebra.presentation(backend="coinvariant")
    presentation_finished = perf_counter()
    rechecked = presentation.verify()
    verify_finished = perf_counter()

    assert presentation.relations == (q**rank,)
    assert presentation.standard_monomials == tuple(q**i for i in range(rank))
    assert presentation.dimension == rank
    assert presentation.hilbert_coefficients == tuple(
        int(degree % 2 == 0) for degree in range(2*rank - 1))
    assert presentation.certificate.certified and rechecked.certified
    assert presentation.backend == "coinvariant"
    assert not algebra._image_cache  # No Schubert image was requested.
    assert len(borel._basis_levels) == 1  # No ambient Schubert basis enumeration.

    return {
        "library_version": __version__,
        "sage_version": SAGE_VERSION,
        "case": "C{}/P{} single squared tail generator".format(rank, rank - 1),
        "input": "q=x{}^2".format(rank),
        "backend": presentation.backend,
        "result": "passed",
        "expected_and_observed_kernel": "(q^{})".format(rank),
        "dimension": presentation.dimension,
        "source_generators": 1,
        "schubert_images_computed": len(algebra._image_cache),
        "schubert_basis_levels_generated": len(borel._basis_levels) - 1,
        "timings_seconds": {
            "library_import": IMPORT_SECONDS,
            "map_and_generator_setup": setup_finished - started,
            "presentation_including_quotient_setup_and_certificate":
                presentation_finished - setup_finished,
            "fresh_quotient_certificate_recheck": verify_finished - presentation_finished,
            "total_calculation_including_setup_and_recheck": verify_finished - started,
            "script_including_imports_through_recheck": verify_finished - SCRIPT_STARTED,
        },
        "timing_scope": "Fresh process and fresh map/algebra; no precomputed image or presentation caches. Interpreter startup before this script is excluded. Hardware-dependent single observation.",
        "previous_observation": "The 0.4.0 BGG version exceeded a 45-second process limit during concurrent validation; that censored observation is not a precise speedup measurement.",
        "presentation": presentation.as_dict(),
        "rechecked_certificate": rechecked.as_dict(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rank", type=int, default=9)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    if arguments.rank < 2:
        parser.error("rank must be at least two")
    report = run(arguments.rank)
    if arguments.output:
        arguments.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ("presentation", "rechecked_certificate")}, indent=2))


if __name__ == "__main__":
    main()
