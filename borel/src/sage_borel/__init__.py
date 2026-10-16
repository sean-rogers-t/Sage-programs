"""Exact, sparse Borel characteristic maps for SageMath.

Run inside SageMath (``sage -python``), rather than ordinary Python.  The
historical scripts elsewhere in this repository are not imported or changed.
"""

from .core import BorelMap
from .expansion import Expansion
from .classical import classical_borel, to_type_a

__version__ = "0.1.0"
__all__ = ["BorelMap", "Expansion", "classical_borel", "to_type_a"]
