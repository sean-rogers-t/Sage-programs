"""Run from borel/: PYTHONPATH=src sage -python examples/basic.py."""

from sage_borel import BorelMap, classical_borel, to_type_a


def main():
    c = classical_borel("C", 12, 10)
    answer = c.tail_elementary(10, 5)
    assert len(answer) == 6
    print("C12, e_5 of the ten squared tail coordinates:")
    print(answer.to_terms())
    print("Frontier sizes:", c.last_stats["frontier_sizes"])

    d = classical_borel("D", 4, 3)
    assert d.tail_elementary(3, 1).to_terms() == [
        {"word": [2, 1], "coefficient": "-1"}]
    print("Corrected D4 elementary image:", d.tail_elementary(3, 1))
    print("D4 Euler image:", d.euler(3))

    a = BorelMap("A3", omitted_nodes=(1, 3))
    w1, w2, w3 = a.weights
    print("A3 with two omitted nodes:", a.expand((w1 + 2*w3)**3 + 7))
    print("A3 top Schubert index:", a.longest_representative().reduced_word())

    e = BorelMap("E7", omitted_nodes=7)
    assert dict(e.divisor_product([7, 7, 7])) == dict(e.expand(e.weights[6]**3))
    print("E7 divisor cube:", e.divisor_product([7, 7, 7]))
    print("D4 fork embedded in type A:",
          to_type_a(d.W.simple_reflections()[4]).reduced_word())


if __name__ == "__main__":
    main()
