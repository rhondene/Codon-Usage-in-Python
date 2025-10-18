"""
Unit tests for core codon usage analysis functions

These tests verify the correctness of RSCU calculations, 
codon frequency analysis, and amino acid usage computations.
"""

import unittest
import pandas as pd
import numpy as np
from codon_usage_gui.core.analysis import (
    compute_codon_frequencies,
    compute_rscu_weights,
    compute_amino_acid_usage,
    compute_per_gene_rscu,
    compute_codon_usage_per_1000,
    compute_relative_codon_frequencies,
    codon_to_aa,
    EmptyDataError,
    InvalidSequenceError
)


class TestCodonAnalysis(unittest.TestCase):
    """Test core codon usage analysis functions."""

    def setUp(self):
        """Set up test fixtures."""
        # Simple test sequences (all valid CDS)
        self.headers = [">gene1", ">gene2", ">gene3"]
        self.sequences = [
            "AUGGCUAAGUAG",  # ATG-GCT-AAG-TAG (Met-Ala-Lys-STOP)
            "AUGUUUGCCUAG",  # ATG-UUU-GCC-TAG (Met-Phe-Ala-STOP)
            "AUGGCCAAAUAG"   # ATG-GCC-AAA-TAG (Met-Ala-Lys-STOP)
        ]
        
        # Test sequence with length not divisible by 3
        self.invalid_length_headers = [">gene1", ">gene2"]
        self.invalid_length_sequences = [
            "AUGGCUAAGUAG",  # Valid (12 bp)
            "AUGGCUAAGUA"    # Invalid (11 bp)
        ]
        
        # Empty sequences
        self.empty_headers = []
        self.empty_sequences = []

    def test_compute_codon_frequencies_valid(self):
        """Test codon frequency computation with valid sequences."""
        df_codcount, skipped = compute_codon_frequencies(self.headers, self.sequences)
        
        # Check return types
        self.assertIsInstance(df_codcount, pd.DataFrame)
        self.assertIsInstance(skipped, list)
        
        # Check DataFrame structure
        expected_columns = ['Codon', 'Obs_Freq', 'Amino_Acid']
        self.assertListEqual(list(df_codcount.columns), expected_columns)
        
        # Check that we have all codons
        self.assertEqual(len(df_codcount), len(codon_to_aa))
        
        # Check that no sequences were skipped
        self.assertEqual(len(skipped), 0)
        
        # Check specific codon counts
        # ATG should appear 3 times (once per sequence)
        atg_count = df_codcount[df_codcount['Codon'] == 'AUG']['Obs_Freq'].iloc[0]
        self.assertEqual(atg_count, 3)
        
        # TAG should appear 3 times (stop codon in all sequences)
        tag_count = df_codcount[df_codcount['Codon'] == 'UAG']['Obs_Freq'].iloc[0]
        self.assertEqual(tag_count, 3)

    def test_compute_codon_frequencies_invalid_length(self):
        """Test codon frequency computation with invalid sequence lengths."""
        df_codcount, skipped = compute_codon_frequencies(
            self.invalid_length_headers, 
            self.invalid_length_sequences
        )
        
        # One sequence should be skipped
        self.assertEqual(len(skipped), 1)
        self.assertIn(">gene2", skipped)
        
        # Should still process the valid sequence
        atg_count = df_codcount[df_codcount['Codon'] == 'AUG']['Obs_Freq'].iloc[0]
        self.assertEqual(atg_count, 1)

    def test_compute_codon_frequencies_empty_input(self):
        """Test codon frequency computation with empty input."""
        with self.assertRaises(EmptyDataError):
            compute_codon_frequencies(self.empty_headers, self.empty_sequences)

    def test_compute_codon_frequencies_mismatched_lengths(self):
        """Test codon frequency computation with mismatched header/sequence lengths."""
        with self.assertRaises(InvalidSequenceError):
            compute_codon_frequencies(self.headers, self.sequences[:2])

    def test_compute_rscu_weights_valid(self):
        """Test RSCU computation with valid codon frequency data."""
        df_codcount, _ = compute_codon_frequencies(self.headers, self.sequences)
        df_rscu = compute_rscu_weights(df_codcount)
        
        # Check DataFrame structure
        expected_columns = ['Codon', 'Obs_Freq', 'Amino_Acid', 'RSCU', 
                           'Relative_Adaptive_Weights', 'optimal']
        self.assertListEqual(list(df_rscu.columns), expected_columns)
        
        # Check that RSCU values are reasonable
        # For codons that appear, RSCU should be > 0
        rscu_values = df_rscu[df_rscu['Obs_Freq'] > 0]['RSCU']
        self.assertTrue(all(rscu_values > 0))
        
        # Check that optimal codons are marked correctly
        # Each amino acid should have at least one optimal codon
        aa_groups = df_rscu.groupby('Amino_Acid')
        for aa, group in aa_groups:
            if group['Obs_Freq'].sum() > 0:  # Only check amino acids with observations
                optimal_count = group['optimal'].sum()
                self.assertGreaterEqual(optimal_count, 1)

    def test_compute_rscu_weights_empty_input(self):
        """Test RSCU computation with empty DataFrame."""
        empty_df = pd.DataFrame()
        with self.assertRaises(EmptyDataError):
            compute_rscu_weights(empty_df)

    def test_compute_rscu_weights_invalid_columns(self):
        """Test RSCU computation with missing required columns."""
        invalid_df = pd.DataFrame({'Wrong_Column': [1, 2, 3]})
        with self.assertRaises(InvalidSequenceError):
            compute_rscu_weights(invalid_df)

    def test_compute_amino_acid_usage_valid(self):
        """Test amino acid usage computation."""
        df_codcount, _ = compute_codon_frequencies(self.headers, self.sequences)
        df_rscu = compute_rscu_weights(df_codcount)
        aa_df = compute_amino_acid_usage(df_rscu)
        
        # Check DataFrame structure
        expected_columns = ['Amino_acid', 'Expected_Freq(%)', 'Obs_Freq(%)', 
                           'Abs_Freq', 'Num_Codons']
        self.assertListEqual(list(aa_df.columns), expected_columns)
        
        # Check that frequencies are reasonable (percentages should be 0-100)
        self.assertTrue(all(aa_df['Expected_Freq(%)'] >= 0))
        self.assertTrue(all(aa_df['Expected_Freq(%)'] <= 100))
        self.assertTrue(all(aa_df['Obs_Freq(%)'] >= 0))
        self.assertTrue(all(aa_df['Obs_Freq(%)'] <= 100))
        
        # Check that observed frequencies sum to 100% (approximately)
        total_obs_freq = aa_df['Obs_Freq(%)'].sum()
        self.assertAlmostEqual(total_obs_freq, 100.0, places=1)

    def test_compute_per_gene_rscu_valid(self):
        """Test per-gene RSCU computation."""
        df_per_gene, skipped = compute_per_gene_rscu(self.headers, self.sequences)
        
        # Check return types
        self.assertIsInstance(df_per_gene, pd.DataFrame)
        self.assertIsInstance(skipped, list)
        
        # Check that we have results for all valid sequences
        unique_genes = df_per_gene['Gene_ID'].unique()
        self.assertEqual(len(unique_genes), 3)
        
        # Check DataFrame structure
        expected_columns = ['Codon', 'Obs_Freq', 'Amino_Acid', 'RSCU', 
                           'Relative_Adaptive_Weights', 'optimal', 'Gene_ID']
        self.assertListEqual(list(df_per_gene.columns), expected_columns)

    def test_compute_codon_usage_per_1000_valid(self):
        """Test codon usage per 1000 computation."""
        df_cu1000, skipped = compute_codon_usage_per_1000(self.headers, self.sequences)
        
        # Check DataFrame structure
        expected_columns = ['Codon', 'Obs_Freq', 'Amino_Acid', 'Usage_per_1000']
        self.assertListEqual(list(df_cu1000.columns), expected_columns)
        
        # Check that usage per 1000 sums to 1000 (approximately)
        total_usage = df_cu1000['Usage_per_1000'].sum()
        self.assertAlmostEqual(total_usage, 1000.0, places=1)
        
        # Check that values are reasonable (0-1000)
        self.assertTrue(all(df_cu1000['Usage_per_1000'] >= 0))
        self.assertTrue(all(df_cu1000['Usage_per_1000'] <= 1000))

    def test_compute_relative_codon_frequencies_valid(self):
        """Test relative codon frequency computation."""
        df_rel_freq, skipped = compute_relative_codon_frequencies(self.headers, self.sequences)
        
        # Check return types
        self.assertIsInstance(df_rel_freq, pd.DataFrame)
        self.assertIsInstance(skipped, list)
        
        # Check that we have one row per gene
        self.assertEqual(len(df_rel_freq), 3)
        
        # Check that Gene_ID column exists
        self.assertIn('Gene_ID', df_rel_freq.columns)
        
        # Check that all codon columns exist
        for codon in codon_to_aa.keys():
            self.assertIn(codon, df_rel_freq.columns)
        
        # Check that frequencies for each gene sum to 1.0 (approximately)
        for _, row in df_rel_freq.iterrows():
            codon_freqs = [row[codon] for codon in codon_to_aa.keys()]
            total_freq = sum(codon_freqs)
            self.assertAlmostEqual(total_freq, 1.0, places=5)

    def test_compute_relative_codon_frequencies_empty_input(self):
        """Test relative codon frequency computation with empty input."""
        with self.assertRaises(EmptyDataError):
            compute_relative_codon_frequencies(self.empty_headers, self.empty_sequences)

    def test_edge_case_single_sequence(self):
        """Test analysis with a single sequence."""
        single_header = [">single_gene"]
        single_sequence = ["AUGGCUAAGUAG"]
        
        df_codcount, skipped = compute_codon_frequencies(single_header, single_sequence)
        self.assertEqual(len(skipped), 0)
        
        df_rscu = compute_rscu_weights(df_codcount)
        self.assertFalse(df_rscu.empty)
        
        aa_df = compute_amino_acid_usage(df_rscu)
        self.assertFalse(aa_df.empty)

    def test_all_stop_codons_present(self):
        """Test that all stop codons are correctly identified."""
        # Sequence with all three stop codons
        headers = [">stop1", ">stop2", ">stop3"]
        sequences = [
            "AUGUAA",  # ATG-UAA
            "AUGUAG",  # ATG-UAG  
            "AUGUGA"   # ATG-UGA
        ]
        
        df_codcount, skipped = compute_codon_frequencies(headers, sequences)
        
        # Check that all stop codons are counted
        stop_codons = ['UAA', 'UAG', 'UGA']
        for stop_codon in stop_codons:
            count = df_codcount[df_codcount['Codon'] == stop_codon]['Obs_Freq'].iloc[0]
            self.assertEqual(count, 1)


if __name__ == '__main__':
    unittest.main()