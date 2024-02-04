#!/usr/bin/python

#Author: Rhondene Wint (rwint@ucmerced.edu, PhD candidate@UC-Merced,California USA)
#Python program that computes the fractional codon content of coding sequences

import argparse
import pandas as pd
import fix_fasta
from pandas.errors import SettingWithCopyWarning
import warnings
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

#dictionary of codons and their corresponding amino acids
codon_aa = {
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
	
##preprocess fasta file into csv

def preproc(fasta_file:str)->pd.DataFrame:
    """Converts fasta formated file to a dataframe"""
    #flybase fasta file has internal newline in the same seqeunce 
    fas=fix_fasta.fix_fasta(fasta_file)
    fasta_df=pd.DataFrame(list(fas))
    fasta_df.columns=['Header','Sequence']
    return fasta_df

def get_cod_freq(geneID:str, cds:str)->pd.DataFrame:
    """computes the fractional frequency of a codon in its CDS
    Input:
        gene: row from fastacsv table
        Output: row of normalized codon frequency """
    codon_count=dict() 
    
    #build dictionary to accumulate codon counts; omit  stop codons
    STOP_CODONS = ['UAA', 'UAG', 'UGA']
    codon_count = {codon: 0 for codon in codon_aa if codon not in STOP_CODONS }
    
    ##count codons in cds
    codons = []
    for c in range(0,len(cds),3):   
        cod=cds[c:c+3]
        if cod in STOP_CODONS:
            continue
        codon_count[cod]+=1
     
     # Calculate fractional frequencies.
    total_cod = sum(codon_count.values()) #total number of codons in the cds
    codon_count = {codon: count / total_cod for codon, count in codon_count.items()}
    
    # Create DataFrame from codon count.
    df_codcnt=pd.DataFrame(list(codon_count.items()) )
    df_codcnt.columns=['Codon', 'Fractional_Freq']
    df_codcnt=df_codcnt.set_index('Codon').T.reset_index(drop=True)
    df_codcnt['GeneID']=geneID
	#reorder columns
    cols2=[df_codcnt.columns[-1]]+sorted(df_codcnt.columns[:61])
    df_codcnt=df_codcnt[cols2]

    return df_codcnt

if __name__ =='__main__':
    #print a welcome message with my contacts
    msg = "Thanks for using CodonCount@2021. \n To send questions and suggestions, you may make a pull request at github OR  email rwint@ucmerced.edu"
    print(msg)
    parser = argparse.ArgumentParser(description="This program takes a fasta file of coding sequences and outputs a csv table of the fractional codon content for each gene. Written by Rhondene Wint (rwint@ucmerced.edu).")
    parser.add_argument('-CDS', help='Name or Path of fasta file with  coding sequences', type=str, required=True, metavar='')
    parser.add_argument('-out', help='Name or Path of destination folder for output file', type=str, metavar='') 
    args=parser.parse_args()
	
    transcriptome=preproc(args.CDS)
    comb=[]  #store the codon content for each gene
    for i in range(transcriptome.shape[0]):
        gene = transcriptome.iloc[i]
        geneID = gene['Header'].split(' ')[0][1:]
        cds = gene['Sequence'].upper().replace('T', 'U') 
        if len(cds) % 3 != 0:
            print(f"WARNING: Skipping {geneID}! Length of CDS is not a multiple of 3.")
            continue
        df_count = get_cod_freq(geneID, cds)
        comb.append(df_count)
    
    pd.concat(comb,axis=0).to_csv(args.out+'.csv',index=False)

    
		
	
	
   
	
    

    