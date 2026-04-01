"""ChronoDecon - Modified Cross-Correlation Algorithm for MS Single-Component Reconstruction.

Based on Anal. Chem. 2021 modification of homologue factor algorithm.
"""

__version__ = "0.1.0"
__author__ = "ChronoDecon Contributors"

from .decon import homologue_factor, deconvolute_mzml

__all__ = ["homologue_factor", "deconvolute_mzml"]
