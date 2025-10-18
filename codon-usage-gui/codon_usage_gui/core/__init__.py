"""
Core analysis modules for codon usage analysis

This module contains all the core computational functions for:
- FASTA file parsing
- Codon frequency calculations
- RSCU (Relative Synonymous Codon Usage) analysis
- Amino acid usage analysis
- Various codon usage metrics
"""

from .analysis import (
    parse_fasta_from_text,
    parse_fasta_file,
    compute_codon_frequencies,
    compute_rscu_weights,
    compute_amino_acid_usage,
    compute_per_gene_rscu,
    compute_codon_usage_per_1000,
    compute_relative_codon_frequencies,
    codon_to_aa
)

__all__ = [
    'parse_fasta_from_text',
    'parse_fasta_file',
    'compute_codon_frequencies',
    'compute_rscu_weights',
    'compute_amino_acid_usage',
    'compute_per_gene_rscu',
    'compute_codon_usage_per_1000',
    'compute_relative_codon_frequencies',
    'codon_to_aa'
]