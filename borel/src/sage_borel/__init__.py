"""Exact, sparse Borel characteristic maps for SageMath.

Run inside SageMath (``sage -python``), rather than ordinary Python.  The
historical scripts elsewhere in this repository are not imported or changed.
"""

from .core import BorelMap
from .expansion import Expansion
from .classical import classical_borel, to_type_a
from .algebra import DegreeAnalysis, GeneratorAlgebra
from .integral import IntegralDegreeAnalysis
from .lattices import (IntegralImageLattice, IntegralImageComparison,
                       IntegralLatticeQuotient, compare_integral_images)
from .orthogonal import OrthogonalSmithSummary, orthogonal_smith
from .presentation import KernelCertificate, KernelPresentation
from .families import TwoStepSymplectic, two_step_symplectic
from .characters import (SpringerCharacterMap, LeadingTermCertificate, TraceSourceCertificate, CharacterObstruction,
                         springer_character_map, polynomial_source_certificate,
                         modular_character_relation, g2_character_algebra)
from .character_invariants import (CharacterInvariantDegree, IntegralCharacterSource,
                                    f4_character_source)
from .minuscule_characters import e6_character_source, e7_character_source
from .character_lattices import IrreducibleCharacterLattice, irreducible_character_lattice
from .geometry import (AdjointQuotient, SpringerGeometry, g2_adjoint_quotient,
                       PowerSumPullback, FactoredCayleyDeterminant,
                       quadratic_character_pullback)
from .e6_geometry import (e6_adjoint_quotient, e6_cayley_determinant,
                          e6_fundamental_characters)
from .e7_geometry import (e7_adjoint_quotient, e7_cayley_determinant,
                          e7_fundamental_characters)
from .degrees import SpringerDegreeCertificate, exceptional_springer_degree
from .relative_images import (RelativeSmithSummary, ExceptionalSourceComparison,
                              ExceptionalRelativeDegree, ExceptionalRootWeightComparison,
                              exceptional_root_weight_comparison)
from .exceptional import (ExceptionalBorel, GeneratorInfo, exceptional_borel, exceptional_algebra,
                          exceptional_generators, exceptional_image,
                          exceptional_polynomial)

__version__ = "0.22.0"
__all__ = ["BorelMap", "Expansion", "classical_borel", "to_type_a",
           "ExceptionalBorel", "GeneratorInfo", "exceptional_borel",
           "exceptional_generators", "exceptional_polynomial", "exceptional_image",
           "exceptional_algebra",
           "DegreeAnalysis", "IntegralDegreeAnalysis", "OrthogonalSmithSummary", "orthogonal_smith",
           "IntegralImageLattice", "IntegralImageComparison", "IntegralLatticeQuotient",
           "compare_integral_images",
           "GeneratorAlgebra", "KernelCertificate", "KernelPresentation",
           "TwoStepSymplectic", "two_step_symplectic",
           "SpringerCharacterMap", "LeadingTermCertificate", "TraceSourceCertificate", "CharacterObstruction",
           "springer_character_map", "polynomial_source_certificate",
           "modular_character_relation", "g2_character_algebra",
           "CharacterInvariantDegree", "IntegralCharacterSource", "f4_character_source",
           "e6_character_source", "e7_character_source",
           "IrreducibleCharacterLattice", "irreducible_character_lattice",
           "AdjointQuotient", "SpringerGeometry", "g2_adjoint_quotient",
           "quadratic_character_pullback", "PowerSumPullback", "FactoredCayleyDeterminant",
           "e6_adjoint_quotient", "e6_cayley_determinant", "e6_fundamental_characters",
           "e7_adjoint_quotient", "e7_cayley_determinant", "e7_fundamental_characters",
           "SpringerDegreeCertificate", "exceptional_springer_degree",
           "RelativeSmithSummary", "ExceptionalSourceComparison",
           "ExceptionalRelativeDegree", "ExceptionalRootWeightComparison",
           "exceptional_root_weight_comparison"]
