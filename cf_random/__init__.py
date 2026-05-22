"""CF-Random: Protein fold-switching and alternative conformation prediction.

A comprehensive package for identifying and analyzing protein fold-switching
events and alternative conformations using AlphaFold predictions and
structural analysis tools.

Key Features:
    - Predict fold-switching regions in proteins
    - Analyze alternative conformations
    - Assess prediction quality via TM-scores
    - Visualize structural variations
    - Blind structural screening without reference structures

Main entry point:
    Use cf_random.main() to execute the prediction workflow.
"""

__version__ = "0.2.0"
__author__ = "Myeongsang (Samuel) Lee, Pramesh Sharma"
__all__ = ["main"]

import logging

# Configure package-level logging
logging.getLogger(__name__).addHandler(logging.NullHandler())

# Import main modules for easier access
try:
    from .core.main import (
        main,
    )
except ImportError as e:
    import warnings

    warnings.warn(f"Failed to import main module: {e}", ImportWarning)
