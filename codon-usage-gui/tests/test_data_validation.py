"""
Unit tests for data validation and error handling

These tests ensure that the package handles edge cases gracefully
and provides appropriate error messages for biologists.
"""

import unittest
import tempfile
import os
import pandas as pd
from codon_usage_gui.core.analysis import (
    parse_fasta_from_text,
    compute_codon_frequencies,
    compute_rscu_weights,
    validate_dna_sequence,
    InvalidSequenceError,
    EmptyDataError,
    CodonUsageError
)


class TestDataValidation(unittest.TestCase):
    """Test data validation and error handling."""

    def test_sequence_validation_comprehensive(self):
        """Test comprehensive sequence validation scenarios."""
        
        # Test valid sequences
        valid_cases = [
            ("ATCG", "AUCG"),
            ("atcg", "AUCG"),  # lowercase
            ("AUCG", "AUCG"),  # already RNA
            ("ATCGATCG", "AUCGAUCG"),
        ]
        
        for input_seq, expected in valid_cases:
            with self.subTest(input_seq=input_seq):
                result = validate_dna_sequence(input_seq, "test")
                self.assertEqual(result, expected)

    def test_sequence_validation_errors(self):
        """Test sequence validation error cases."""
        
        # Test invalid characters
        invalid_cases = [
            "ATCGX",     # Invalid character X
            "ATCG123",   # Numbers
            "ATCG!@#",   # Special characters
            "ATCGRYSWKM", # IUPAC ambiguity codes (not N)
        ]
        
        for invalid_seq in invalid_cases:
            with self.subTest(invalid_seq=invalid_seq):
                with self.assertRaises(InvalidSequenceError):
                    validate_dna_sequence(invalid_seq, "test")

    def test_sequence_validation_edge_cases(self):
        """Test sequence validation edge cases."""
        
        # Empty sequence
        with self.assertRaises(InvalidSequenceError):
            validate_dna_sequence("", "test")
        
        # None input
        with self.assertRaises(InvalidSequenceError):
            validate_dna_sequence(None, "test")
        
        # Sequence with N characters (should be removed)
        result = validate_dna_sequence("ATCGNATCGN", "test")
        self.assertEqual(result, "AUCGAUCG")
        
        # Only N characters (should result in empty sequence)
        result = validate_dna_sequence("NNNN", "test")
        self.assertEqual(result, "")

    def test_fasta_parsing_malformed_input(self):
        """Test FASTA parsing with malformed input."""
        
        # No header before sequence
        malformed_fasta = """ATCGATCG
>gene1
ATCGATCG"""
        with self.assertRaises(InvalidSequenceError):
            parse_fasta_from_text(malformed_fasta)
        
        # Header without sequence followed by another header
        empty_seq_fasta = """>gene1
>gene2
ATCGATCG"""
        headers, seqs = parse_fasta_from_text(empty_seq_fasta)
        # Should only get gene2 since gene1 has no sequence
        self.assertEqual(len(headers), 1)
        self.assertEqual(headers[0], ">gene2")

    def test_fasta_parsing_whitespace_handling(self):
        """Test FASTA parsing with various whitespace scenarios."""
        
        # Leading/trailing whitespace
        whitespace_fasta = """  >gene1  
  ATCGATCG  
  >gene2  
  ATCGATCG  """
        headers, seqs = parse_fasta_from_text(whitespace_fasta)
        self.assertEqual(len(headers), 2)
        self.assertEqual(headers[0], ">gene1")
        self.assertEqual(seqs[0], "AUCGAUCG")
        
        # Sequences split across multiple lines
        multiline_fasta = """>gene1
ATCG
ATCG
ATCG
>gene2
ATCGATCGATCG"""
        headers, seqs = parse_fasta_from_text(multiline_fasta)
        self.assertEqual(len(seqs), 2)
        self.assertEqual(seqs[0], "AUCGAUCGAUCG")
        self.assertEqual(seqs[1], "AUCGAUCGAUCG")

    def test_codon_analysis_error_propagation(self):
        """Test that errors propagate correctly through analysis pipeline."""
        
        # Empty input should raise EmptyDataError
        with self.assertRaises(EmptyDataError):
            compute_codon_frequencies([], [])
        
        # Mismatched headers and sequences
        with self.assertRaises(InvalidSequenceError):
            compute_codon_frequencies([">gene1"], ["ATCG", "ATCG"])

    def test_rscu_computation_edge_cases(self):
        """Test RSCU computation with edge cases."""
        
        # DataFrame with all zero frequencies
        zero_freq_df = pd.DataFrame({
            'Codon': ['AUG', 'UUU'],
            'Obs_Freq': [0, 0],
            'Amino_Acid': ['Met', 'Phe']
        })
        
        # Should not raise an error, but RSCU should be 0
        result = compute_rscu_weights(zero_freq_df)
        self.assertTrue(all(result['RSCU'] == 0))
        self.assertTrue(all(result['Relative_Adaptive_Weights'] == 0))
        self.assertTrue(all(result['optimal'] == False))

    def test_file_operations_error_handling(self):
        """Test file operation error handling."""
        
        # Test with non-existent file
        with self.assertRaises(FileNotFoundError):
            from codon_usage_gui.core.analysis import parse_fasta_file
            parse_fasta_file("/this/path/does/not/exist.fasta")
        
        # Test with file that has permission issues (simulated)
        # Create a temporary file and then try to read it after deletion
        with tempfile.NamedTemporaryFile(mode='w', delete=True) as f:
            temp_path = f.name
        
        # File is now deleted, should raise FileNotFoundError
        with self.assertRaises(FileNotFoundError):
            from codon_usage_gui.core.analysis import parse_fasta_file
            parse_fasta_file(temp_path)

    def test_robust_header_parsing(self):
        """Test robust header parsing for gene ID extraction."""
        
        # Headers with various formats
        complex_headers = [
            ">gene1 description text",
            ">gene2|annotation|more_info",
            ">gene3",
            ">complex_gene_name_123 [organism=E.coli] length=300"
        ]
        sequences = ["ATGATGATG"] * 4  # Valid 9bp sequences
        
        df_codcount, skipped = compute_codon_frequencies(complex_headers, sequences)
        self.assertEqual(len(skipped), 0)  # All should be processed

    def test_mixed_dna_rna_input(self):
        """Test handling of mixed DNA/RNA input."""
        
        mixed_fasta = """>dna_seq
ATCGATCG
>rna_seq
AUCGAUCG
>mixed_seq
ATCGAUCG"""
        
        headers, seqs = parse_fasta_from_text(mixed_fasta)
        
        # All should be converted to RNA format
        for seq in seqs:
            self.assertNotIn('T', seq)
            self.assertTrue(all(base in 'AUCG' for base in seq))

    def test_very_short_sequences(self):
        """Test handling of very short sequences."""
        
        short_sequences = [
            ("ATG", True),    # 3bp - valid
            ("AT", False),    # 2bp - invalid
            ("ATGC", False),  # 4bp - invalid
            ("", False),      # empty - invalid
        ]
        
        for seq, should_be_valid in short_sequences:
            headers = [f">seq_{len(seq)}bp"]
            sequences = [seq] if seq else []
            
            if should_be_valid and seq:
                try:
                    df_codcount, skipped = compute_codon_frequencies(headers, [seq])
                    # Should process successfully
                    self.assertEqual(len(skipped), 0)
                except EmptyDataError:
                    # Empty sequences will raise this error
                    pass
            else:
                if seq:  # Only test non-empty sequences for skipping
                    df_codcount, skipped = compute_codon_frequencies(headers, [seq])
                    self.assertEqual(len(skipped), 1)

    def test_large_dataset_simulation(self):
        """Test with a larger simulated dataset."""
        
        # Create 100 sequences of varying lengths (all multiples of 3)
        headers = [f">gene_{i}" for i in range(100)]
        sequences = []
        
        for i in range(100):
            # Create sequences of length 12, 15, or 18
            length = 12 + (i % 3) * 3
            seq = "ATG" * (length // 3)
            sequences.append(seq)
        
        # Should process all sequences without error
        df_codcount, skipped = compute_codon_frequencies(headers, sequences)
        self.assertEqual(len(skipped), 0)
        
        # ATG should appear length/3 times per sequence * 100 sequences
        expected_atg_count = sum(len(seq) // 3 for seq in sequences)
        actual_atg_count = df_codcount[df_codcount['Codon'] == 'AUG']['Obs_Freq'].iloc[0]
        self.assertEqual(actual_atg_count, expected_atg_count)

    def test_error_message_quality(self):
        """Test that error messages are informative for users."""
        
        try:
            validate_dna_sequence("ATCGXYZ", "test_gene")
        except InvalidSequenceError as e:
            error_msg = str(e)
            # Check that error message contains useful information
            self.assertIn("test_gene", error_msg)
            self.assertIn("invalid characters", error_msg.lower())
            self.assertIn("X", error_msg)
        
        try:
            compute_codon_frequencies([], [])
        except EmptyDataError as e:
            error_msg = str(e)
            self.assertIn("No sequences", error_msg)


if __name__ == '__main__':
    unittest.main()