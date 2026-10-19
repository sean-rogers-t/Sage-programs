"""Check the proved full D two-tail partition rule and recorded exact matrices.

Run with ordinary Python: python borel/experiments/verify_type_d_two_tail_gluing.py
The partition enumeration is compared with the independent two-row B rule,
the unit correction, and previously recorded full integer Smith matrices.
No Sage, Borel engine, or other verifier is imported. Finite checks are
corroboration of TYPE_D_TWO_TAIL_GLUING.md, not its all-rank proof.
"""

from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter


def full_partition_counts(a, degree):
    if degree < 0 or degree % 2:
        return Counter()
    answer = Counter()
    for second in range(degree // 2 + 1):
        first = degree - second
        if first > 2*a+1 or first == second > a:
            continue
        answer[int(first > a)+int(second > a)] += 1
    return answer


def two_tail_b_counts(a, degree):
    """The separate B quadric/two-row partition rules, including point case."""
    if degree < 0 or degree % 2:
        return Counter()
    m = degree // 2
    if a == 0:
        return Counter({0: 1}) if m == 0 else Counter()
    if a == 1:
        return Counter({int(m > 1): 1}) if m <= 2 else Counter()
    upper, lower = (a+1)//2, a//2
    answer = Counter()
    for second in range(m//2+1):
        first = m-second
        if first > a:
            continue
        exponent = (0 if first <= upper and second <= lower
                    else 2 if second > upper else 1)
        answer[exponent] += 1
    return answer


def expected_counts(a, degree):
    answer = two_tail_b_counts(a, degree)+two_tail_b_counts(a-1, degree-2)
    overlap = int(degree % 2 == 0 and a//2+1 <= degree//2 <= a)
    answer[0] -= overlap
    answer[1] += overlap
    assert all(value >= 0 for value in answer.values())
    return +answer, overlap


def exponents(counts):
    return sorted(exponent for exponent, count in counts.items()
                  for _ in range(count))


def main():
    start = perf_counter()
    here = Path(__file__).resolve()
    borel = here.parent.parent
    inputs = [here, borel/'TYPE_D_TWO_TAIL_GLUING.md',
              here.with_name('orthogonal_sector_search_results.json'),
              here.with_name('type_d_frobenius_search_results.json'),
              here.with_name('type_d_even_gluing_verification.json')]
    hashes = {p.relative_to(borel).as_posix(): sha256(p.read_bytes()).hexdigest()
              for p in inputs}
    rows = []
    for a in range(1, 41):
        total_rank = 0
        for degree in range(0, 4*a+1, 2):
            actual = full_partition_counts(a, degree)
            expected, overlap = expected_counts(a, degree)
            assert actual == expected, (a, degree, actual, expected)
            sector = two_tail_b_counts(a, degree)+two_tail_b_counts(a-1, degree-2)
            assert actual[2] == sector[2]
            index_difference = (sum(e*c for e, c in actual.items())
                                -sum(e*c for e, c in sector.items()))
            assert index_difference == overlap
            total_rank += sum(actual.values())
            rows.append({'a': a, 'k': 2, 'polynomial_degree': degree,
                         'partition_exponents': exponents(actual),
                         'sector_formula_matches': True,
                         'saturated_gluing_dimension': overlap})
        assert total_rank == (a+1)**2, (a, total_rank)

    recorded = []

    def check(a, degree, actual, source):
        predicted = exponents(full_partition_counts(a, degree))
        assert predicted == actual, (a, degree, predicted, actual, source)
        recorded.append({'a': a, 'k': 2, 'polynomial_degree': degree,
                         'source': source, 'actual_exponents': actual,
                         'partition_rule_matches': True})

    historical = json.loads(inputs[2].read_text())
    for case in historical['cases']:
        if case['k'] != 2:
            continue
        for row in case['degrees']:
            if row['status'] == 'computed':
                check(case['a'], row['polynomial_degree'],
                      row['combined']['power_of_two_exponents'], inputs[2].name)
    fresh = json.loads(inputs[3].read_text())
    for row in fresh['fresh_D_checks']:
        if row['k'] == 2 and row['status'] == 'computed':
            check(row['a'], row['polynomial_degree'],
                  row['actual_combined_exponents'], inputs[3].name)
    latest = json.loads(inputs[4].read_text())
    for row in latest['three_head_components']+latest['higher_Cartesian_components']:
        if row['k'] == 2:
            check(row['a'], row['polynomial_degree'], row['full_exponents'], inputs[4].name)
    assert all(sha256(p.read_bytes()).hexdigest() == hashes[p.relative_to(borel).as_posix()]
               for p in inputs)
    report = {'schema_version': 1, 'status': 'passed',
              'generated_utc': datetime.now(timezone.utc).isoformat(),
              'scope': 'Finite checks of the all-rank two-tail proof; no theorem inferred from finite data.',
              'source_sha256': hashes, 'partition_checks': rows,
              'recorded_exact_matrix_checks': recorded,
              'summary': {'head_sizes_checked': 40,
                          'degree_components_checked': len(rows),
                          'recorded_exact_full_matrix_comparisons': len(recorded),
                          'counterexamples': 0},
              'seconds': round(perf_counter()-start, 6)}
    here.with_name('type_d_two_tail_gluing_verification.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode('utf-8'))
    print(json.dumps({'status': 'passed', **report['summary'], 'seconds': report['seconds']}))


if __name__ == '__main__':
    main()
