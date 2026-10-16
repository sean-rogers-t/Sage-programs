"""Exact, sparse Borel characteristic maps for SageMath.

Run inside SageMath (``sage -python``), rather than ordinary Python.  The
historical scripts elsewhere in this repository are not imported or changed.
"""

from .core import BorelMap
from .expansion import Expansion
from .classical import classical_borel, to_type_a
from .algebra import DegreeAnalysis, GeneratorAlgebra
from .families import TwoStepSymplectic, two_step_symplectic
from .exceptional import (ExceptionalBorel, GeneratorInfo, exceptional_borel,
                          exceptional_generators, exceptional_image,
                          exceptional_polynomial)

__version__ = "0.3.0"
__all__ = ["BorelMap", "Expansion", "classical_borel", "to_type_a",
           "ExceptionalBorel", "GeneratorInfo", "exceptional_borel",
           "exceptional_generators", "exceptional_polynomial", "exceptional_image",
           "DegreeAnalysis", "GeneratorAlgebra", "TwoStepSymplectic", "two_step_symplectic"]
