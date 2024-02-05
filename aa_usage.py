#author: Rhondene Wint
# To compute Amino acid usage from FASTA of Coding sequences

"""If a particular amino acid is in some way adaptive, then it should occur more frequently than expected by chance. This can easily be tested by calculating the expected frequencies of amino acids and comparing to observed. The codons and observed frequencies of particular amino acids are given in the table.
- The frequencies of DNA bases in nature are 22.0% uracil, 30.3% adenine, 21.7% cytosine, and 26.1% guanine. The expected frequency of a particular codon can then be calculated by multiplying the frequencies of each DNA base comprising the codon. The expected frequency of the amino acid can then be calculated by adding the frequencies of each codon that codes for that amino acid.
- As an example, the RNA codons for tyrosine are UAU and UAC, so the random expectation for its frequency is (0.220)(0.303)(0.220) + (0.220)(0.303)(0.217) = 0.0292. Since 3 of the 64 codons are nonsense or stop codons, this frequency for each amino acid is multiplied by a correction factor of 1.057."""
import pandas as pd
import numpy as np
import argparse
import sys


def get_seqs(fasta: str) -> list:
    """Parse fasta file into a list of coding sequences, handling internal newlines."""
    seqs = []; headers=[]
    current_seq = ''  
    with open(fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                headers.append(line.strip())
                if current_seq:  
                    seqs.append(current_seq)
                    current_seq = ''  # Reset current sequence
                continue  # Skip the header line
            else:
                current_seq += line.strip()  
        if current_seq:  
            seqs.append(current_seq)
    return headers, seqs


## dictionary that maps codons to amino acids
"""break up 6-codon family into 2 and 4 fold (Ser (S), L(Leu), R (Arg))"""
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
    "GGU":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly"}


def get_cod_freq(headers:list, seqs:list)->pd.DataFrame:
    """ seqs: list of CDS
    """
    codon_count=dict() 
    codon_count = {codon: 0 for codon in codon_to_aa}
    for i,cds in enumerate(seqs):
        if len(cds)%3 !=0:
            ID = headers[i].split(' ')[0]
            print(f"WARNING: Skipping {ID} Length of CDS is not a multiple of 3.")
            continue
        cds = cds.upper().replace('T','U')
        ##count codons in cds
        for i in range(0, len(cds), 3):
            codon = cds[i:i+3]
            if codon in codon_count:
                codon_count[codon] += 1
    df_codcount=pd.DataFrame(list(codon_count.items()) )
    df_codcount.columns=['Codon', 'Obs_Freq']
    df_codcount['Amino_Acid'] = [codon_to_aa[codon] for codon in df_codcount['Codon'].values] ## add amino acid column
    
    return df_codcount

def compute_rscu_weights(df_codcount:pd.DataFrame)->pd.DataFrame:
    """ Caclculates Relative Synonymous codon usage (RSCU) """
    aa_groups = df_codcount.groupby('Amino_Acid')
    df_list = []
    for a, group in aa_groups:
        d = group.copy()
        d['RSCU'] = d['Obs_Freq'] / d['Obs_Freq'].mean()
        d['Relative_Adaptive_Weights'] = d['RSCU'] / d['RSCU'].max()
        d['optimal'] = d['RSCU'] == d['RSCU'].max()
        df_list.append(d)
    return pd.concat(df_list)

## to consolidate all codons of 6-box amino acids; you may or may not want to do this
def no_six(aa):
    return aa[:3]

def compute_aa_usage(species: str, df_rscu: pd.DataFrame)->pd.DataFrame:
    """Computes and saves amino acid usage."""
    base_freq = {'U': 0.220, "A": 0.303, 'C': 0.217, 'G': 0.261}  # Base frequency in nature
    df_rscu = df_rscu[df_rscu['Amino_Acid'] != 'STOP']
    aa_usage = {}
    
    for amino, group in df_rscu.groupby('Amino_Acid'):
        expected_aa_usage = 0
        for codon in group['Codon'].values:
            expected_aa_usage+= base_freq[codon[0]]*base_freq[codon[1]]*base_freq[codon[2]]
        expected_aa_usage = expected_aa_usage*1.057  ##correction factor
        obs_freq = group['Obs_Freq'].sum() / df_rscu['Obs_Freq'].sum()
        aa_usage[amino] = [expected_aa_usage*100, obs_freq * 100, group['Obs_Freq'].sum(), len(group['Codon'].unique())]

    aa_df = pd.DataFrame.from_dict(aa_usage, orient='index', columns=['Expected_Freq(%)', 'Obs_Freq(%)', 'Abs_Freq', 'Num_Codons']).reset_index()
    aa_df.columns = ['Amino_acid', 'Expected_Freq(%)', 'Obs_Freq(%)', 'Abs_Freq', 'Num_Codons']
    
    return aa_df
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Computes Observed and Expected Amino Acid Frequencies + RSCU from FASTA of Coding Sequences according to https://qubeshub.org/publications/979/serve/1/3067?el=1&download=1 .")
    parser.add_argument("-CDS", help="Path to the input FASTA file.", type=str, required=True, metavar='')
    parser.add_argument("-out", help="Prefix for the output files.",  type=str, required=True,metavar='')
    parser.add_argument("-keep_sixfold_aa", help=" True|False Whether to consolidate sixfold degenerate amino acids or break up into 4fold and 2fold. Default is False.", default=False, metavar='')
    args = parser.parse_args()

    headers,seqs = get_seqs(args.CDS)
    df_rscu = get_cod_freq(headers, seqs)  ##computes absolute codon frequencies
    rscu = compute_rscu_weights(df_rscu)  ##computes RSCU and adaptive weights
    aa_df = compute_aa_usage(args.out, rscu)  ##computes amino acid usage

    rscu.to_csv(f'{args.out}_rscu.csv', index=False, sep=',')
    aa_df.to_csv(f'{args.out}_aa_usage.csv', index=False, sep=',')