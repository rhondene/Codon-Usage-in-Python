"""
Core analysis functions for codon usage analysis
Adapted from the original codon usage analysis tools
"""

import pandas as pd
import numpy as np
import io
import os
import logging
from typing import Tuple, List, Optional


# Set up logging for better error reporting
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class CodonUsageError(Exception):
    """Custom exception for codon usage analysis errors."""
    pass


class InvalidSequenceError(CodonUsageError):
    """Raised when sequences contain invalid characters or format."""
    pass


class EmptyDataError(CodonUsageError):
    """Raised when no valid sequences are found for analysis."""
    pass


def validate_dna_sequence(sequence: str, sequence_id: str = "Unknown") -> str:
    """
    Validate and clean a DNA sequence.
    
    Args:
        sequence: DNA sequence string
        sequence_id: Identifier for the sequence (for error reporting)
        
    Returns:
        Cleaned and validated sequence
        
    Raises:
        InvalidSequenceError: If sequence contains invalid characters
    """
    if not sequence or not isinstance(sequence, str):
        raise InvalidSequenceError(f"Sequence {sequence_id}: Empty or invalid sequence provided")
    
    # Remove whitespace and convert to uppercase
    clean_seq = sequence.strip().upper()
    
    # Check for valid DNA characters (including T and U)
    valid_chars = set('ATCGUN')
    invalid_chars = set(clean_seq) - valid_chars
    
    if invalid_chars:
        raise InvalidSequenceError(
            f"Sequence {sequence_id}: Contains invalid characters: {', '.join(sorted(invalid_chars))}. "
            f"Only A, T, C, G, U, N are allowed."
        )
    
    # Convert T to U for RNA codon usage
    clean_seq = clean_seq.replace('T', 'U')
    
    # Warn about N characters
    if 'N' in clean_seq:
        logger.warning(f"Sequence {sequence_id}: Contains {clean_seq.count('N')} ambiguous nucleotides (N)")
        # Remove N characters for analysis
        clean_seq = clean_seq.replace('N', '')
    
    return clean_seq


def parse_fasta_from_text(fasta_text: str) -> Tuple[List[str], List[str]]:
    """
    Parse fasta text into a list of headers and sequences.
    
    Args:
        fasta_text: String containing FASTA format data
        
    Returns:
        Tuple of (headers, sequences)
        
    Raises:
        InvalidSequenceError: If FASTA format is invalid
        EmptyDataError: If no valid sequences found
    """
    if not fasta_text or not isinstance(fasta_text, str):
        raise EmptyDataError("No FASTA text provided")
    
    headers = []
    seqs = []
    current_seq = ''
    current_header = None

    lines = fasta_text.strip().split('\n')
    
    if not lines:
        raise EmptyDataError("Empty FASTA input")
    
    for line_num, line in enumerate(lines, 1):
        line = line.strip()
        if not line:  # Skip empty lines
            continue
            
        if line.startswith('>'):
            # Save previous sequence if exists
            if current_header is not None and current_seq:
                try:
                    validated_seq = validate_dna_sequence(current_seq, current_header)
                    headers.append(current_header)
                    seqs.append(validated_seq)
                except InvalidSequenceError as e:
                    logger.warning(f"Skipping sequence due to error: {e}")
            
            current_header = line.strip()
            current_seq = ''
        else:
            if current_header is None:
                raise InvalidSequenceError(f"Line {line_num}: Sequence data found before header")
            current_seq += line.strip()

    # Add the last sequence
    if current_header is not None and current_seq:
        try:
            validated_seq = validate_dna_sequence(current_seq, current_header)
            headers.append(current_header)
            seqs.append(validated_seq)
        except InvalidSequenceError as e:
            logger.warning(f"Skipping sequence due to error: {e}")

    if not headers or not seqs:
        raise EmptyDataError("No valid sequences found in FASTA input")

    logger.info(f"Successfully parsed {len(seqs)} sequences from FASTA text")
    return headers, seqs


def parse_fasta_file(file_path: str) -> Tuple[List[str], List[str]]:
    """
    Parse fasta file into a list of headers and sequences.
    
    Args:
        file_path: Path to FASTA file
        
    Returns:
        Tuple of (headers, sequences)
        
    Raises:
        FileNotFoundError: If file doesn't exist
        InvalidSequenceError: If FASTA format is invalid
        EmptyDataError: If no valid sequences found
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"FASTA file not found: {file_path}")
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            fasta_text = f.read()
    except UnicodeDecodeError:
        try:
            # Try with different encoding
            with open(file_path, 'r', encoding='latin-1') as f:
                fasta_text = f.read()
        except Exception as e:
            raise InvalidSequenceError(f"Unable to read file {file_path}: {e}")
    except Exception as e:
        raise InvalidSequenceError(f"Error reading file {file_path}: {e}")
    
    return parse_fasta_from_text(fasta_text)


# Codon to amino acid mapping
codon_to_aa = {
    "UUU":"Phe", "UUC":"Phe",
    "UCU":"Ser4", "UCC":"Ser4", "UCA":"Ser4", "UCG":"Ser4",
    "AGU":"Ser2", "AGC":"Ser2",
    "CUU":"Leu4", "CUC":"Leu4", "CUA":"Leu4", "CUG":"Leu4",
    "UUA":"Leu2", "UUG":"Leu2",

    "UAU":"Tyr", "UAC":"Tyr", "UAA":"STOP", "UAG":"STOP",
    "UGU":"Cys", "UGC":"Cys", "UGA":"STOP", "UGG":"Trp",
    "CGU":"Arg4", "CGC":"Arg4", "CGA":"Arg4", "CGG":"Arg4",
    "AGA":"Arg2", "AGG":"Arg2",
    "CCU":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",
    "CAU":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln",

    "AUU":"Ile", "AUC":"Ile", "AUA":"Ile", "AUG":"Met",
    "ACU":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",
    "AAU":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys",

    "GUU":"Val", "GUC":"Val", "GUA":"Val", "GUG":"Val",
    "GCU":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",
    "GAU":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu",
    "GGU":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly"
}


def compute_codon_frequencies(headers: List[str], seqs: List[str]) -> Tuple[pd.DataFrame, List[str]]:
    """
    Compute codon frequencies from sequences.
    
    Args:
        headers: List of sequence headers
        seqs: List of DNA/RNA sequences
        
    Returns:
        Tuple of (codon frequency DataFrame, list of skipped sequence IDs)
        
    Raises:
        EmptyDataError: If no valid sequences provided
        InvalidSequenceError: If sequences are invalid
    """
    if not seqs or not headers:
        raise EmptyDataError("No sequences provided for analysis")
    
    if len(headers) != len(seqs):
        raise InvalidSequenceError("Number of headers and sequences must match")
    
    codon_count = {codon: 0 for codon in codon_to_aa}
    skipped_sequences = []
    processed_count = 0

    for i, cds in enumerate(seqs):
        try:
            sequence_id = headers[i].split(' ')[0] if i < len(headers) else f"Sequence_{i}"
            
            # Validate sequence length
            if len(cds) == 0:
                logger.warning(f"Skipping empty sequence: {sequence_id}")
                skipped_sequences.append(sequence_id)
                continue
                
            if len(cds) % 3 != 0:
                logger.warning(f"Skipping sequence {sequence_id}: length {len(cds)} is not divisible by 3")
                skipped_sequences.append(sequence_id)
                continue

            # Ensure sequence is in RNA format
            cds_rna = cds.upper().replace('T', 'U')

            # Count codons
            codon_found = False
            for j in range(0, len(cds_rna), 3):
                codon = cds_rna[j:j+3]
                if len(codon) == 3:  # Ensure complete codon
                    if codon in codon_count:
                        codon_count[codon] += 1
                        codon_found = True
                    else:
                        logger.warning(f"Unknown codon '{codon}' in sequence {sequence_id} at position {j}")
            
            if codon_found:
                processed_count += 1
                
        except Exception as e:
            sequence_id = headers[i].split(' ')[0] if i < len(headers) else f"Sequence_{i}"
            logger.error(f"Error processing sequence {sequence_id}: {e}")
            skipped_sequences.append(sequence_id)

    if processed_count == 0:
        raise EmptyDataError("No valid sequences could be processed")
    
    logger.info(f"Processed {processed_count} sequences, skipped {len(skipped_sequences)}")

    df_codcount = pd.DataFrame(list(codon_count.items()))
    df_codcount.columns = ['Codon', 'Obs_Freq']
    df_codcount['Amino_Acid'] = [codon_to_aa[codon] for codon in df_codcount['Codon'].values]

    return df_codcount, skipped_sequences


def compute_rscu_weights(df_codcount: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate Relative Synonymous Codon Usage (RSCU).
    
    Args:
        df_codcount: DataFrame with codon frequencies
        
    Returns:
        DataFrame with RSCU values and relative adaptive weights
        
    Raises:
        EmptyDataError: If input DataFrame is empty
        InvalidSequenceError: If DataFrame has wrong format
    """
    if df_codcount.empty:
        raise EmptyDataError("No codon frequency data provided")
    
    required_columns = ['Codon', 'Obs_Freq', 'Amino_Acid']
    missing_columns = [col for col in required_columns if col not in df_codcount.columns]
    if missing_columns:
        raise InvalidSequenceError(f"Missing required columns: {missing_columns}")
    
    total_codons = df_codcount['Obs_Freq'].sum()
    if total_codons == 0:
        raise EmptyDataError("No codons found in the dataset")
    
    logger.info(f"Computing RSCU for {total_codons} total codons")
    
    aa_groups = df_codcount.groupby('Amino_Acid')
    df_list = []

    for a, group in aa_groups:
        d = group.copy()
        mean_freq = d['Obs_Freq'].mean()
        
        if mean_freq > 0:
            d['RSCU'] = d['Obs_Freq'] / mean_freq
            max_rscu = d['RSCU'].max()
            if max_rscu > 0:
                d['Relative_Adaptive_Weights'] = d['RSCU'] / max_rscu
            else:
                d['Relative_Adaptive_Weights'] = 0
            d['optimal'] = d['RSCU'] == max_rscu
        else:
            d['RSCU'] = 0
            d['Relative_Adaptive_Weights'] = 0
            d['optimal'] = False
        df_list.append(d)

    result_df = pd.concat(df_list, ignore_index=True)
    logger.info(f"RSCU computation completed for {len(result_df)} codons")
    return result_df


def compute_amino_acid_usage(df_rscu: pd.DataFrame) -> pd.DataFrame:
    """
    Compute amino acid usage statistics.
    
    Args:
        df_rscu: DataFrame with RSCU values
        
    Returns:
        DataFrame with amino acid usage statistics
        
    Raises:
        EmptyDataError: If input DataFrame is empty
        InvalidSequenceError: If DataFrame has wrong format
    """
    if df_rscu.empty:
        raise EmptyDataError("No RSCU data provided")
    
    required_columns = ['Codon', 'Obs_Freq', 'Amino_Acid']
    missing_columns = [col for col in required_columns if col not in df_rscu.columns]
    if missing_columns:
        raise InvalidSequenceError(f"Missing required columns: {missing_columns}")
    
    base_freq = {'U': 0.220, "A": 0.303, 'C': 0.217, 'G': 0.261}
    df_rscu = df_rscu[df_rscu['Amino_Acid'] != 'STOP']
    aa_usage = {}

    total_codons = df_rscu['Obs_Freq'].sum()
    if total_codons == 0:
        raise EmptyDataError("No valid codons found for amino acid analysis")

    logger.info(f"Computing amino acid usage for {len(df_rscu['Amino_Acid'].unique())} amino acids")

    for amino, group in df_rscu.groupby('Amino_Acid'):
        expected_aa_usage = 0
        for codon in group['Codon'].values:
            if len(codon) == 3:  # Ensure valid codon length
                try:
                    expected_aa_usage += base_freq[codon[0]] * base_freq[codon[1]] * base_freq[codon[2]]
                except KeyError as e:
                    logger.warning(f"Unknown nucleotide in codon {codon}: {e}")
                    continue

        expected_aa_usage = expected_aa_usage * 1.057  # correction factor
        obs_freq = group['Obs_Freq'].sum() / total_codons if total_codons > 0 else 0

        aa_usage[amino] = [
            expected_aa_usage * 100,
            obs_freq * 100,
            group['Obs_Freq'].sum(),
            len(group['Codon'].unique())
        ]

    aa_df = pd.DataFrame.from_dict(
        aa_usage,
        orient='index',
        columns=['Expected_Freq(%)', 'Obs_Freq(%)', 'Abs_Freq', 'Num_Codons']
    ).reset_index()
    aa_df.columns = ['Amino_acid', 'Expected_Freq(%)', 'Obs_Freq(%)', 'Abs_Freq', 'Num_Codons']

    logger.info(f"Amino acid usage computed for {len(aa_df)} amino acids")
    return aa_df


def compute_per_gene_rscu(headers: List[str], seqs: List[str]) -> Tuple[pd.DataFrame, List[str]]:
    """
    Compute RSCU for each gene individually.
    
    Args:
        headers: List of sequence headers
        seqs: List of DNA/RNA sequences
        
    Returns:
        Tuple of (per-gene RSCU DataFrame, list of skipped sequence IDs)
        
    Raises:
        EmptyDataError: If no valid sequences provided
    """
    if not seqs or not headers:
        raise EmptyDataError("No sequences provided for per-gene analysis")
    
    results = []
    skipped_sequences = []
    processed_count = 0

    logger.info(f"Computing per-gene RSCU for {len(seqs)} sequences")

    for i, cds in enumerate(seqs):
        try:
            sequence_id = headers[i].split(' ')[0] if i < len(headers) else f"Sequence_{i}"
            
            if len(cds) % 3 != 0:
                logger.warning(f"Skipping sequence {sequence_id}: length not divisible by 3")
                skipped_sequences.append(sequence_id)
                continue

            if len(cds) == 0:
                logger.warning(f"Skipping empty sequence: {sequence_id}")
                skipped_sequences.append(sequence_id)
                continue

            cds = cds.upper().replace('T', 'U')

            # Count codons for this sequence
            codon_count = {codon: 0 for codon in codon_to_aa}
            valid_codons = 0
            
            for j in range(0, len(cds), 3):
                codon = cds[j:j+3]
                if len(codon) == 3 and codon in codon_count:
                    codon_count[codon] += 1
                    valid_codons += 1

            if valid_codons == 0:
                logger.warning(f"No valid codons found in sequence {sequence_id}")
                skipped_sequences.append(sequence_id)
                continue

            # Create dataframe for this sequence
            df_seq = pd.DataFrame(list(codon_count.items()))
            df_seq.columns = ['Codon', 'Obs_Freq']
            df_seq['Amino_Acid'] = [codon_to_aa[codon] for codon in df_seq['Codon'].values]

            # Compute RSCU for this sequence
            rscu_seq = compute_rscu_weights(df_seq)
            rscu_seq['Gene_ID'] = sequence_id

            results.append(rscu_seq)
            processed_count += 1
            
        except Exception as e:
            sequence_id = headers[i].split(' ')[0] if i < len(headers) else f"Sequence_{i}"
            logger.error(f"Error processing sequence {sequence_id}: {e}")
            skipped_sequences.append(sequence_id)

    if not results:
        raise EmptyDataError("No sequences could be processed for per-gene RSCU analysis")
    
    logger.info(f"Per-gene RSCU computed for {processed_count} sequences")
    final_df = pd.concat(results, ignore_index=True)
    return final_df, skipped_sequences


def compute_codon_usage_per_1000(headers: List[str], seqs: List[str]) -> Tuple[pd.DataFrame, List[str]]:
    """
    Compute codon usage per 1000 codons.
    
    Args:
        headers: List of sequence headers
        seqs: List of DNA/RNA sequences
        
    Returns:
        Tuple of (codon usage DataFrame, list of skipped sequence IDs)
        
    Raises:
        EmptyDataError: If no valid sequences provided
    """
    df_codcount, skipped = compute_codon_frequencies(headers, seqs)

    total_codons = df_codcount['Obs_Freq'].sum()
    if total_codons > 0:
        df_codcount['Usage_per_1000'] = (df_codcount['Obs_Freq'] / total_codons) * 1000
        logger.info(f"Computed usage per 1000 for {total_codons} total codons")
    else:
        df_codcount['Usage_per_1000'] = 0
        logger.warning("No codons found, usage per 1000 set to 0")

    return df_codcount, skipped


def compute_relative_codon_frequencies(headers: List[str], seqs: List[str]) -> Tuple[pd.DataFrame, List[str]]:
    """
    Compute relative codon frequencies for each gene.
    
    Args:
        headers: List of sequence headers
        seqs: List of DNA/RNA sequences
        
    Returns:
        Tuple of (relative frequency DataFrame, list of skipped sequence IDs)
        
    Raises:
        EmptyDataError: If no valid sequences provided
    """
    if not seqs or not headers:
        raise EmptyDataError("No sequences provided for relative frequency analysis")
    
    results = []
    skipped_sequences = []
    processed_count = 0

    logger.info(f"Computing relative codon frequencies for {len(seqs)} sequences")

    for i, cds in enumerate(seqs):
        try:
            sequence_id = headers[i].split(' ')[0] if i < len(headers) else f"Sequence_{i}"
            
            if len(cds) % 3 != 0:
                logger.warning(f"Skipping sequence {sequence_id}: length not divisible by 3")
                skipped_sequences.append(sequence_id)
                continue

            if len(cds) == 0:
                logger.warning(f"Skipping empty sequence: {sequence_id}")
                skipped_sequences.append(sequence_id)
                continue

            cds = cds.upper().replace('T', 'U')

            # Count codons for this sequence
            codon_count = {codon: 0 for codon in codon_to_aa}
            for j in range(0, len(cds), 3):
                codon = cds[j:j+3]
                if len(codon) == 3 and codon in codon_count:
                    codon_count[codon] += 1

            total_codons_in_seq = sum(codon_count.values())
            
            if total_codons_in_seq == 0:
                logger.warning(f"No valid codons found in sequence {sequence_id}")
                skipped_sequences.append(sequence_id)
                continue

            # Create row for this sequence
            row = {'Gene_ID': sequence_id}
            for codon, count in codon_count.items():
                row[codon] = count / total_codons_in_seq if total_codons_in_seq > 0 else 0

            results.append(row)
            processed_count += 1
            
        except Exception as e:
            sequence_id = headers[i].split(' ')[0] if i < len(headers) else f"Sequence_{i}"
            logger.error(f"Error processing sequence {sequence_id}: {e}")
            skipped_sequences.append(sequence_id)

    if not results:
        raise EmptyDataError("No sequences could be processed for relative frequency analysis")
    
    logger.info(f"Relative frequencies computed for {processed_count} sequences")
    return pd.DataFrame(results), skipped_sequences