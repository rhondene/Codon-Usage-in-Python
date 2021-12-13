"""Compute Relative Synonymous Codon Usage (RSCU) of  a transcriptome(Sharp and Li 1986)"""

import pandas as pd
import numpy as np


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


## read in fasta genome file and parse into list

def get_seqs(species):
    seqs = []
    with open(species, 'r') as f:
        for line in f:
            if line.startswith('>') is False:
                seqs.append(line.rstrip())
    return seqs


def get_cod_freq(seqs):
    """ seqs: list of CDS
    	Returns a 59-dim dataframe of total absolute codon frequencies
    """
    
    codon_count=dict() 

    for codon in list(codon_to_aa.keys()):
        codon_count[codon]=0  ##dictionary to accumulate codon count
        
    for cds in seqs:
        cds = cds.upper().replace('T','U')
        codons = []
        ##make a list of codons
        for c in range(0,len(cds),3):
            if len(cds)%3 ==0:
                cod=cds[c:c+3]
                if 'N' not in cod:  ##ignore N and seqs not multiple of 3
                    codons.append(cod)
            else:
                continue

        for c in list(codon_count.keys()):
            codon_count[c]+= codons.count(c)
    
    df_codcount=pd.DataFrame(list(codon_count.items()) )
    df_codcount.columns=['Codon', 'Obs_Freq']
    df_codcount['Amino_Acid'] = [codon_to_aa[codon] for codon in df_codcount['Codon'].values] ## add amino acid column
    
    return df_codcount

def compute_rscu_weights(df_codcount):
    """ Caclculates Relative Synonymous codon usage (RSCU) wij = RSCUij/ RSCU i,max
    Input: 59-dim codon count dataframe
    Returns: 59-dim dataframe of RSCU values for each codon """
    aa_groups = df_codcount.groupby('Amino_Acid')
    aa =  df_codcount['Amino_Acid'].unique()  #make a list of all amino acids to iterate over
    df_list = []
    for a in aa:
        d=aa_groups.get_group(a)
        d['RSCU'] = d['Obs_Freq'].values/d['Obs_Freq'].mean() #obs/expected freq 
        d['Relative_Adaptive_Weights'] = d['RSCU'].values/d['RSCU'].max() 
        d['optimal'] = [True if rscu==d['RSCU'].max() else False for rscu in d['RSCU'].values] #marks optimal codon
        d['Species'] = species
        df_list.append(d)
    return pd.concat(df_list)


if __name__=='__main__':
    species= input('Enter fasta file address: ')
    out_name= input('Enter output file name: ')
    seqs = get_seqs(species)  ##formats fasta into list of sequences
    df_rscu = get_cod_freq(seqs)  ##computes absolute codon frequencies
    rscu = compute_rscu_weights(df_rscu)  ##computes RSCU and adaptive weights
    rscu.to_csv(f'{out_name}.csv', index=False)

   
                      
       
    
