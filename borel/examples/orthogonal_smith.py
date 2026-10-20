"""Fast formula reports; run inside Sage with PYTHONPATH=src."""

from sage_borel import orthogonal_smith


def main():
    result = orthogonal_smith("D", 4, 2, 4)
    assert result.multiplicities == {1: 1, 2: 2}
    assert result.unit_gluing_dimension == 1
    print("D4, k=2, H^8 factors:", result.smith_factors())
    print("Index in selected saturation:", result.saturation_index())

    character = orthogonal_smith("D", 4, 2, 4, normalization="character")
    assert character.smith_factors() == (16, 32, 32)
    print("Explicit character normalization:", character.smith_factors())

    for degree in (14, 16):
        rank_nine = orthogonal_smith("D", 9, 4, degree)
        print("D9, k=4, polynomial degree", degree,
              "multiplicities:", rank_nine.multiplicities)

    large = orthogonal_smith("D", 200, 80, 1000)
    print("D200 selected image rank:", large.rank)
    print("Its saturation index is 2^E, with E =", large.saturation_exponent)
    print("Only", len(large.exponent_multiplicities), "nonzero multiplicities stored.")


if __name__ == "__main__":
    main()
