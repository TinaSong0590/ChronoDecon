"""
ChronoDecon - Chronological Deconvolution for Mass Spectrometry

A Python package implementing the modified correlation algorithm from the 2021 Analytical Chemistry paper
for single component reconstruction in mass spectrometry data.
"""

from .decon import homologue_factor, deconvolute_mzml

__version__ = "0.1.0"
__all__ = ["homologue_factor", "deconvolute_mzml"]
