"""
CF-random: Predicting alternative conformations and fold-switching proteins

A package for identifying and analyzing protein fold-switching and alternative conformations
using AlphaFold predictions and structural analysis tools.
"""

__version__ = "0.1.0"
__author__ = "Myeongsang (Samuel) Lee"
__all__ = [
    "main",
]

# Import main modules for easier access
try:
    from .core import main
except ImportError:
    pass
