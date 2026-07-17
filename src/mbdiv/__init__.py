"""
mbdiv — Microbiome Diversity Analysis Pipeline

A single-command pipeline for species-level metagenomic diversity analysis.

    raw_data  ->  merge species  ->  normalize  ->  alpha diversity
                                          |
              beta diversity  <-  relative abundance
                                          |
              filter viruses/fungi  ->  Top-N stacked bar

Quick start:
    mbdiv species.xlsx metadata.xlsx
    mbdiv species.xlsx metadata.xlsx -o results --top 10
"""

__version__ = "1.0.0"
__author__ = "WorkBuddy"

from .config import PipelineConfig

__all__ = ["PipelineConfig", "__version__"]
