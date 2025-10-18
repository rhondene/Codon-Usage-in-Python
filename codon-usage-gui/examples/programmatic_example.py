"""
Example Python script showing how to use codon-usage-gui programmatically

This script demonstrates how biologists can use the package functions
directly in their own Python scripts for batch processing or integration
with other bioinformatics pipelines.
"""

# Import the analysis functions
from codon_usage_gui import (
    parse_fasta_file,
    compute_codon_frequencies,
    compute_rscu_weights,
    compute_amino_acid_usage,
    compute_per_gene_rscu
)
import pandas as pd


def analyze_codon_usage(fasta_file_path, output_prefix="analysis"):
    """
    Complete codon usage analysis pipeline.
    
    Args:
        fasta_file_path: Path to FASTA file
        output_prefix: Prefix for output files
    """
    print(f"Analyzing {fasta_file_path}...")
    
    # Step 1: Parse FASTA file
    try:
        headers, sequences = parse_fasta_file(fasta_file_path)
        print(f"Loaded {len(sequences)} sequences")
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return
    
    # Step 2: Compute codon frequencies
    try:
        df_codons, skipped = compute_codon_frequencies(headers, sequences)
        if skipped:
            print(f"Warning: Skipped {len(skipped)} sequences: {skipped[:3]}...")
    except Exception as e:
        print(f"Error computing codon frequencies: {e}")
        return
    
    # Step 3: Compute RSCU values
    try:
        df_rscu = compute_rscu_weights(df_codons)
        print("RSCU computation completed")
    except Exception as e:
        print(f"Error computing RSCU: {e}")
        return
    
    # Step 4: Compute amino acid usage
    try:
        df_aa = compute_amino_acid_usage(df_rscu)
        print("Amino acid usage analysis completed")
    except Exception as e:
        print(f"Error computing amino acid usage: {e}")
        return
    
    # Step 5: Per-gene analysis (optional for large datasets)
    try:
        if len(sequences) <= 100:  # Only for smaller datasets
            df_per_gene, skipped_per_gene = compute_per_gene_rscu(headers, sequences)
            print("Per-gene RSCU analysis completed")
        else:
            print("Skipping per-gene analysis (too many sequences)")
            df_per_gene = None
    except Exception as e:
        print(f"Error in per-gene analysis: {e}")
        df_per_gene = None
    
    # Step 6: Save results
    try:
        # Save codon frequencies and RSCU
        df_rscu.to_csv(f"{output_prefix}_rscu.csv", index=False)
        print(f"RSCU results saved to {output_prefix}_rscu.csv")
        
        # Save amino acid usage
        df_aa.to_csv(f"{output_prefix}_amino_acids.csv", index=False)
        print(f"Amino acid results saved to {output_prefix}_amino_acids.csv")
        
        # Save per-gene results if available
        if df_per_gene is not None:
            df_per_gene.to_csv(f"{output_prefix}_per_gene.csv", index=False)
            print(f"Per-gene results saved to {output_prefix}_per_gene.csv")
        
    except Exception as e:
        print(f"Error saving results: {e}")
    
    # Step 7: Print summary statistics
    print("\n=== ANALYSIS SUMMARY ===")
    print(f"Total sequences processed: {len(sequences) - len(skipped)}")
    print(f"Total codons analyzed: {df_codons['Obs_Freq'].sum()}")
    
    # Find most and least used codons
    used_codons = df_rscu[df_rscu['Obs_Freq'] > 0]
    if not used_codons.empty:
        most_used = used_codons.loc[used_codons['RSCU'].idxmax()]
        print(f"Most preferred codon: {most_used['Codon']} (RSCU = {most_used['RSCU']:.2f})")
        
        least_used = used_codons.loc[used_codons['RSCU'].idxmin()]
        print(f"Least preferred codon: {least_used['Codon']} (RSCU = {least_used['RSCU']:.2f})")
    
    # Codon bias strength (coefficient of variation of RSCU)
    rscu_std = df_rscu[df_rscu['Obs_Freq'] > 0]['RSCU'].std()
    rscu_mean = df_rscu[df_rscu['Obs_Freq'] > 0]['RSCU'].mean()
    bias_strength = rscu_std / rscu_mean if rscu_mean > 0 else 0
    print(f"Codon bias strength (CV): {bias_strength:.3f}")
    
    print("Analysis completed!")
    
    return df_rscu, df_aa, df_per_gene


def quick_codon_bias_check(fasta_file_path):
    """
    Quick check for codon bias in a dataset.
    
    Args:
        fasta_file_path: Path to FASTA file
        
    Returns:
        Dictionary with bias metrics
    """
    try:
        headers, sequences = parse_fasta_file(fasta_file_path)
        df_codons, skipped = compute_codon_frequencies(headers, sequences)
        df_rscu = compute_rscu_weights(df_codons)
        
        # Calculate bias metrics
        used_codons = df_rscu[df_rscu['Obs_Freq'] > 0]
        
        metrics = {
            'total_sequences': len(sequences),
            'processed_sequences': len(sequences) - len(skipped),
            'total_codons': df_codons['Obs_Freq'].sum(),
            'rscu_mean': used_codons['RSCU'].mean(),
            'rscu_std': used_codons['RSCU'].std(),
            'bias_strength': used_codons['RSCU'].std() / used_codons['RSCU'].mean(),
            'highly_biased_codons': len(used_codons[used_codons['RSCU'] > 1.5]),
            'avoided_codons': len(used_codons[used_codons['RSCU'] < 0.5])
        }
        
        return metrics
        
    except Exception as e:
        print(f"Error in bias check: {e}")
        return None


if __name__ == "__main__":
    # Example usage
    import os
    
    # Get the directory where this script is located
    script_dir = os.path.dirname(__file__)
    
    # Example 1: Full analysis
    sample_file = os.path.join(script_dir, "sample_sequences.fasta")
    if os.path.exists(sample_file):
        print("Running full analysis on sample sequences...")
        analyze_codon_usage(sample_file, "sample_analysis")
        print("\n" + "="*50 + "\n")
    
    # Example 2: Quick bias check
    ecoli_file = os.path.join(script_dir, "ecoli_sample.fasta")
    if os.path.exists(ecoli_file):
        print("Running quick bias check on E. coli sequences...")
        metrics = quick_codon_bias_check(ecoli_file)
        if metrics:
            print("Codon Bias Metrics:")
            for key, value in metrics.items():
                print(f"  {key}: {value}")
    
    print("\nTo use with your own data:")
    print("  python programmatic_example.py")
    print("Or import functions in your own scripts:")
    print("  from codon_usage_gui import parse_fasta_file, compute_rscu_weights")