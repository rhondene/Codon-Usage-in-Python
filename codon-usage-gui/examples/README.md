# Codon Usage Analysis Examples

This directory contains example FASTA files and usage guides for the codon-usage-gui package.

## Example Files

### 1. `sample_sequences.fasta`
A small set of 5 genes with different expression patterns:
- High expression genes (optimized codons)
- Medium expression genes 
- Low expression genes
- Housekeeping genes
- Tissue-specific genes

**Usage:**
```bash
# After installing the package
codon-usage-gui
# Then upload this file in the GUI
```

### 2. `ecoli_sample.fasta`
Example bacterial (E. coli) coding sequences showing typical bacterial codon usage patterns.

## Quick Start Guide 

### Step 1: Installation
```bash
pip install codon-usage-gui
```

### Step 2: Launch the GUI
```bash
codon-usage-gui
```

### Step 3: Upload Your Data
1. Click "Upload FASTA file" in the sidebar
2. Select your FASTA file containing coding sequences
3. Or copy/paste sequences in FASTA format

### Step 4: Choose Analysis Type
- **Transcriptome-wide RSCU**: Overall codon usage patterns across all genes
- **Per-gene RSCU**: Individual gene codon usage patterns
- **Amino Acid Usage**: Compare expected vs observed amino acid frequencies
- **Codon Usage per 1000**: Normalized codon frequencies
- **Relative Codon Frequencies**: Gene-by-gene relative usage

### Step 5: Interpret Results

#### RSCU Values (Relative Synonymous Codon Usage)
- **RSCU = 1**: Codon used at expected frequency
- **RSCU > 1**: Codon used more than expected (preferred)
- **RSCU < 1**: Codon used less than expected (avoided)
- **RSCU >> 1**: Highly preferred codon (often in highly expressed genes)

#### What to Look For:
1. **Highly expressed genes**: Should show strong codon bias (some codons with high RSCU)
2. **Lowly expressed genes**: Should show more uniform codon usage (RSCU closer to 1)
3. **GC-rich organisms**: Prefer codons ending in G or C
4. **AT-rich organisms**: Prefer codons ending in A or T

## Data Requirements

### Input Format
Your FASTA file should contain:
- Coding sequences (CDS) only
- Complete codons (sequence length divisible by 3)
- Standard nucleotide codes (A, T, C, G, U)

### Example FASTA Format
```
>gene_name_1
ATGGCAAGCGAATTTGCCGAAGCCCTGGACAAAGCG...
>gene_name_2 [optional description]
ATGGCTCTGGAAATTGCAGAGGCTCTGGATAAAACT...
```

### What NOT to Include
- ❌ Genomic sequences with introns
- ❌ UTR regions (5' or 3' untranslated regions)
- ❌ Partial sequences (incomplete codons)
- ❌ Sequences with ambiguous nucleotides (except N, which will be removed)

## Common Issues and Solutions

### Issue: "Skipped sequences (not multiple of 3)"
**Solution:** Your sequences contain incomplete codons. Check that:
- Sequences represent complete coding regions
- No extra nucleotides at start/end
- Sequences start with start codon (ATG) and end with stop codon

### Issue: "No valid sequences found"
**Solution:** Check your FASTA format:
- Each sequence must have a header starting with `>`
- Headers and sequences must be on separate lines
- Sequences should contain only A, T, C, G, U nucleotides

### Issue: Empty or zero results
**Solution:** 
- Ensure sequences are long enough (at least one complete codon)
- Check for invalid characters in sequences
- Verify sequences represent protein-coding genes

## Biological Interpretation Guide

### For Different Organism Types:

#### Prokaryotes (Bacteria)
- Expect strong codon bias in highly expressed genes
- Look for AT-rich or GC-rich bias depending on organism
- Ribosomal protein genes typically show strong bias

#### Eukaryotes
- Generally weaker codon bias than prokaryotes
- Highly expressed genes (histones, ribosomal proteins) show strongest bias
- Tissue-specific genes may show intermediate bias

#### Fast-growing organisms
- Strong selection for "optimal" codons
- Clear RSCU peaks for preferred codons
- Correlation between expression level and codon bias

#### Slowly-growing organisms
- Weaker codon bias overall
- More uniform RSCU values
- Less correlation between expression and bias

## Export and Further Analysis

All results can be exported as CSV files for further analysis in:
- Excel or Google Sheets
- R statistical software  
- Python pandas
- Other bioinformatics tools

The exported data includes:
- Codon frequencies and RSCU values
- Amino acid usage statistics
- Per-gene analysis results