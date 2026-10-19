"""Independent formal checks of the relative D compensated-parity lemma.

Run with ordinary Python 3 from any directory. No Sage or Borel engine is
imported. The base Euler symbol Z is retained, rather than set to zero.
This checks the algebraic projection used in the incidence proof; it is
not a computational proof of equivariant Giambelli or the geometric steps.
"""

from datetime import datetime, timezone
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from time import perf_counter


def add_term(polynomial, exponent, coefficient):
    if coefficient:
        polynomial[exponent] = polynomial.get(exponent, 0) + coefficient
        if not polynomial[exponent]:
            del polynomial[exponent]


def multiply_monomial(polynomial, exponent, coefficient=1):
    return {tuple(a+b for a, b in zip(key, exponent)): value*coefficient
            for key, value in polynomial.items() if value*coefficient}


def average(polynomial, heads):
    """Average all actual double sign substitutions xi,p -> -xi,-p."""
    total = {}
    for signs in product((0, 1), repeat=heads):
        for exponent, coefficient in polynomial.items():
            parity = sum(signs[i]*(exponent[i]+exponent[heads])
                         for i in range(heads))
            add_term(total, exponent, coefficient*(-1)**parity)
    order = 2**heads
    assert all(value % order == 0 for value in total.values())
    return {key: value//order for key, value in total.items()}


def coefficient_projection(polynomial, heads):
    return {exponent: coefficient for exponent, coefficient in polynomial.items()
            if all((exponent[i]-exponent[heads]) % 2 == 0
                   for i in range(heads))}


def normalize_target(exponent, heads):
    """Use (u1...ut)q=Z^2 in Z[u1,...,ut,q,Z]."""
    exponent = list(exponent)
    count = min(exponent[:heads]+[exponent[heads]])
    for i in range(heads+1):
        exponent[i] -= count
    exponent[-1] += 2*count
    return tuple(exponent)


def invariant_to_target(polynomial, heads):
    """Use p^2=q and (x1...xt)p=Z, then the top squared relation."""
    result = {}
    for exponent, coefficient in polynomial.items():
        x, p, q, z = exponent[:heads], *exponent[heads:]
        assert all((degree-p) % 2 == 0 for degree in x)
        odd = p % 2
        target = tuple((degree-odd)//2 for degree in x)+(q+p//2, z+odd)
        add_term(result, normalize_target(target, heads), coefficient)
    return result


def run_parameter(heads):
    # Source coordinates are x1,...,xt,p,q,Z; target coordinates u1,...,ut,q,Z.
    length = heads+3
    zero = (0,)*length
    euler_product = (1,)*heads+(1, 0, 0)
    base_euler = (0,)*heads+(0, 0, 1)
    euler_relation = {euler_product: 1, base_euler: -1}
    p_square = (0,)*heads+(2, 0, 0)
    tail_square = (0,)*heads+(0, 1, 0)
    square_relation = {p_square: 1, tail_square: -1}
    count = 0
    for powers in product(range(4), repeat=heads+1):
        for q, z in product(range(2), repeat=2):
            exponent = powers+(q, z)
            polynomial = {exponent: 1}
            projected = coefficient_projection(polynomial, heads)
            assert average(polynomial, heads) == projected
            assert coefficient_projection(projected, heads) == projected
            for relation in (euler_relation, square_relation):
                multiple = multiply_monomial(relation, exponent)
                projected_multiple = average(multiple, heads)
                assert projected_multiple == coefficient_projection(multiple, heads)
                assert invariant_to_target(projected_multiple, heads) == {}
            count += 1

    # The Euler invariant stays nonzero: rho((product xi)p)=rho(Z)=Z.
    assert invariant_to_target(average({euler_product: 1}, heads), heads) == {
        (0,)*(heads+1)+(1,): 1}
    assert invariant_to_target(average({base_euler: 1}, heads), heads) == {
        (0,)*(heads+1)+(1,): 1}

    # A fixed head-square polynomial is preserved with no division by two.
    target_checks = 0
    for powers in product(range(3), repeat=heads):
        exponent = tuple(2*degree for degree in powers)+(0, 0, 0)
        assert average({exponent: 1}, heads) == {exponent: 1}
        assert invariant_to_target({exponent: 1}, heads) == {powers+(0, 0): 1}
        target_checks += 1

    # Uncompensated head averaging kills (product xi)p, but fixes Z.
    # Its image of the Euler relation would be -Z, which is nonzero.
    bad = {base_euler: -1}
    assert invariant_to_target(bad, heads) == {(0,)*(heads+1)+(1,): -1}
    assert average({zero: 1}, heads) == {zero: 1}
    return {"head_size": heads, "monomial_checks": count,
            "Euler_relation_multiples": count,
            "p_square_relation_multiples": count,
            "head_square_monomials_fixed": target_checks,
            "nonzero_base_Euler_preserved": True,
            "uncompensated_projection_counterexample_confirmed": True}


def main():
    started = perf_counter()
    here = Path(__file__).resolve()
    proof = here.parent.parent/"TYPE_D_SQUARE_THEOREM.md"
    sources = (here, proof)
    hashes = {path.relative_to(here.parent.parent).as_posix(): sha256(path.read_bytes()).hexdigest()
              for path in sources}
    cases = [run_parameter(heads) for heads in range(1, 5)]
    report = {"schema_version": 1, "status": "passed",
              "generated_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "Formal compensated-parity/Euler algebra only; no geometric or all-rank Smith theorem inferred from finite checks.",
              "independence": "Pure Python integer polynomial dictionaries; no Sage, Borel, or previous verifier imports.",
              "source_sha256": hashes, "cases": cases,
              "summary": {"head_sizes": len(cases),
                          "monomial_checks": sum(row["monomial_checks"] for row in cases),
                          "Euler_relation_multiples": sum(row["Euler_relation_multiples"] for row in cases),
                          "p_square_relation_multiples": sum(row["p_square_relation_multiples"] for row in cases),
                          "head_square_monomials_fixed": sum(row["head_square_monomials_fixed"] for row in cases),
                          "nonzero_base_Euler_checks": len(cases),
                          "false_uncompensated_projection_checks": len(cases)},
              "seconds": round(perf_counter()-started, 6)}
    assert all(sha256(path.read_bytes()).hexdigest() == hashes[path.relative_to(here.parent.parent).as_posix()]
               for path in sources)
    output = here.with_name("type_d_incidence_projection_verification.json")
    output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"status": "passed", **report["summary"], "seconds": report["seconds"]}))


if __name__ == "__main__":
    main()
