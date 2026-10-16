"""Run from borel/: PYTHONPATH=src sage -python examples/exceptional.py."""

from sage_borel import (exceptional_borel, exceptional_generators,
                        exceptional_image)


def main():
    print("F4/P1 generators:",
          [(g.name, g.polynomial_degree) for g in exceptional_generators("F4", 1)])
    print("F4/P1 e2:", exceptional_image("F4", 1, "e2").to_terms())

    e6 = exceptional_borel("E6", 6)
    print("E6/P6 supplementary Euler input:", e6.generator_polynomial("euler"))
    print("E6/P6 supplementary Euler image:", e6.generator_image("euler").to_terms())

    e7 = exceptional_borel("E7", 7)
    print("E7/P7 generators:", [(g.name, g.polynomial_degree) for g in e7.generators()])
    print("E7/P7 psi5:", e7.generator_image("psi5").to_terms())

    e8 = exceptional_borel("E8", 2)
    print("E8 catalogue coverage:", dict(e8.catalogue_coverage))
    print("E8/P2 e4 (corrected dissertation label):", e8.generator_image("e4").to_terms())


if __name__ == "__main__":
    main()
