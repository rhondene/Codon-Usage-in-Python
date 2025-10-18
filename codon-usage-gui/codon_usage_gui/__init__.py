"""
Codon Usage GUI - A Streamlit application for codon usage bias analysis

This package provides tools for analyzing codon usage bias in coding sequences.
It includes both a web-based GUI (using Streamlit) and command-line interface.

Main components:
- core.analysis: Core analysis functions for codon usage calculations
- app: Streamlit web interface
- cli: Command-line interface for launching the GUI
"""

__version__ = "0.1.0"
__author__ = "Rhondene Wint"

# Make core functions available at package level for easier importing
try:
    from .core.analysis import (
        parse_fasta_from_text,
        parse_fasta_file,
        compute_codon_frequencies,
        compute_rscu_weights,
        compute_amino_acid_usage,
        compute_per_gene_rscu,
        compute_codon_usage_per_1000,
        compute_relative_codon_frequencies
    )
except ImportError:
    # Handle import errors gracefully during installation
    pass

# Package metadata
__all__ = [
    'parse_fasta_from_text',
    'parse_fasta_file', 
    'compute_codon_frequencies',
    'compute_rscu_weights',
    'compute_amino_acid_usage',
    'compute_per_gene_rscu',
    'compute_codon_usage_per_1000',
    'compute_relative_codon_frequencies'
]