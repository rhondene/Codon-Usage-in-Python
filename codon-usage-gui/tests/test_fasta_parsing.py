"""
Unit tests for FASTA parsing functions

These tests ensure that FASTA files and text are correctly parsed
and that appropriate errors are raised for invalid inputs.
"""

import unittest
import tempfile
import os
from codon_usage_gui.core.analysis import (
    parse_fasta_from_text,
    parse_fasta_file,
    validate_dna_sequence,
    InvalidSequenceError,
    EmptyDataError
)


class TestFastaParsing(unittest.TestCase):
    """Test FASTA parsing functionality."""

    def setUp(self):
        """Set up test fixtures."""
        self.valid_fasta_text = """>gene1
ATGGCTAAGTAG
>gene2
ATGTTTGCCTAG
>gene3
ATGGCCAAATAG"""

        self.valid_fasta_with_whitespace = """>gene1
ATG GCT AAG TAG
>gene2 description here
ATG TTT 
GCC TAG
>gene3
ATGGCCAAATAG"""

        self.invalid_fasta_text = """>gene1
ATGGCTXXXTAG
>gene2
ATGTTTGCCTAG"""

        self.mixed_valid_invalid = """>gene1
ATGGCTAAGTAG
>gene2_invalid
ATGGCTAAGTA
>gene3
ATGTTTGCCTAG"""

    def test_parse_valid_fasta_text(self):
        """Test parsing of valid FASTA text."""
        headers, seqs = parse_fasta_from_text(self.valid_fasta_text)
        
        self.assertEqual(len(headers), 3)
        self.assertEqual(len(seqs), 3)
        self.assertEqual(headers[0], ">gene1")
        self.assertEqual(seqs[0], "AUGGCUAAGUAG")  # Should be converted to RNA
        self.assertEqual(seqs[1], "AUGUUUGCCUAG")
        self.assertEqual(seqs[2], "AUGGCCAAAUAG")

    def test_parse_fasta_with_whitespace(self):
        """Test parsing FASTA with whitespace in sequences."""
        headers, seqs = parse_fasta_from_text(self.valid_fasta_with_whitespace)
        
        self.assertEqual(len(headers), 3)
        self.assertEqual(len(seqs), 3)
        self.assertEqual(seqs[0], "AUGGCUAAGUAG")  # Whitespace should be removed
        self.assertEqual(seqs[1], "AUGUUUGCCUAG")  # Newlines should be handled

    def test_parse_empty_fasta_text(self):
        """Test parsing of empty FASTA text."""
        with self.assertRaises(EmptyDataError):
            parse_fasta_from_text("")
        
        with self.assertRaises(EmptyDataError):
            parse_fasta_from_text(None)

    def test_parse_fasta_no_sequences(self):
        """Test FASTA with headers but no sequences."""
        fasta_text = """>gene1
>gene2"""
        with self.assertRaises(EmptyDataError):
            parse_fasta_from_text(fasta_text)

    def test_parse_fasta_sequence_before_header(self):
        """Test FASTA with sequence before any header."""
        fasta_text = """ATGGCTAAGTAG
>gene1
ATGTTTGCCTAG"""
        with self.assertRaises(InvalidSequenceError):
            parse_fasta_from_text(fasta_text)

    def test_validate_dna_sequence_valid(self):
        """Test DNA sequence validation with valid sequences."""
        valid_seq = validate_dna_sequence("ATCGATCG", "test_seq")
        self.assertEqual(valid_seq, "AUCGAUCG")  # T should be converted to U
        
        valid_seq_with_u = validate_dna_sequence("AUCGAUCG", "test_seq")
        self.assertEqual(valid_seq_with_u, "AUCGAUCG")

    def test_validate_dna_sequence_invalid_chars(self):
        """Test DNA sequence validation with invalid characters."""
        with self.assertRaises(InvalidSequenceError):
            validate_dna_sequence("ATCGXYZ", "test_seq")
        
        with self.assertRaises(InvalidSequenceError):
            validate_dna_sequence("ATCG123", "test_seq")

    def test_validate_dna_sequence_empty(self):
        """Test DNA sequence validation with empty sequences."""
        with self.assertRaises(InvalidSequenceError):
            validate_dna_sequence("", "test_seq")
        
        with self.assertRaises(InvalidSequenceError):
            validate_dna_sequence(None, "test_seq")

    def test_validate_dna_sequence_with_n(self):
        """Test DNA sequence validation with ambiguous nucleotides."""
        # Should remove N characters
        result = validate_dna_sequence("ATCGNATCG", "test_seq")
        self.assertEqual(result, "AUCGAUCG")

    def test_parse_fasta_file_valid(self):
        """Test parsing a valid FASTA file."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as f:
            f.write(self.valid_fasta_text)
            temp_path = f.name
        
        try:
            headers, seqs = parse_fasta_file(temp_path)
            self.assertEqual(len(headers), 3)
            self.assertEqual(len(seqs), 3)
            self.assertEqual(seqs[0], "AUGGCUAAGUAG")
        finally:
            os.unlink(temp_path)

    def test_parse_fasta_file_not_found(self):
        """Test parsing a non-existent FASTA file."""
        with self.assertRaises(FileNotFoundError):
            parse_fasta_file("/path/that/does/not/exist.fasta")

    def test_parse_mixed_valid_invalid_sequences(self):
        """Test parsing FASTA with mix of valid and invalid sequences."""
        headers, seqs = parse_fasta_from_text(self.mixed_valid_invalid)
        
        # Should have 2 valid sequences (gene1 and gene3)
        self.assertEqual(len(headers), 2)
        self.assertEqual(len(seqs), 2)
        self.assertEqual(headers[0], ">gene1")
        self.assertEqual(headers[1], ">gene3")

    def test_parse_fasta_case_insensitive(self):
        """Test that parsing handles both upper and lowercase."""
        fasta_text = """>gene1
atggctaagtag
>gene2
ATGTTTGCCTAG"""
        headers, seqs = parse_fasta_from_text(fasta_text)
        
        self.assertEqual(len(headers), 2)
        self.assertEqual(seqs[0], "AUGGCUAAGUAG")  # Should be uppercase and RNA
        self.assertEqual(seqs[1], "AUGUUUGCCUAG")


if __name__ == '__main__':
    unittest.main()