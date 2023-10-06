#!/usr/bin/python

#Author: Rhondene Wint (rwint@ucmerced.edu, PhD candidate@UC-Merced,California USA)
#Python program that computes the fractional codon content of coding sequences


import argparse
import pandas as pd
import fix_fasta
from pandas.errors import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)


#codon_aa_table
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

def preproc(fasta_file):
    """Converts fasta formated file to a dataframe"""
    
    #flybase fasta file has internal newline in the same seqeunce 
    fas=fix_fasta.fix_fasta(fasta_file)
    fasta_df=pd.DataFrame(list(fas))
    fasta_df.columns=['Header','Sequence']

    return fasta_df


def get_cod_freq(gene):
    """computes the fractional frequency of a codon in its CDS
    Input:
    gene: row from fastacsv table"""
    header = gene.iloc[:,0].values[0].split(' ')
    geneID=header[0][1:]


    #get coding sequence
    cds = gene.iloc[:,1].values[0].upper().replace('T','U')
    codon_count=dict() 
    
    #build dictionary to accumulate codon counts; ignore with stop codons
    for codon in list(codon_aa.keys()):
        if codon not in [ "UAA","UAG", "UGA" ]:
            codon_count[codon]=0
   
    ##count codons in cds
    codons = []
    for c in range(0,len(cds),3):   #O(len cds)
        cod=cds[c:c+3]
        try:
            codon_count[cod]+=1
        except KeyError:
            continue
        
    #store the fractional freq of each codon in the codon dictionary
    total_cod=len(cds)/3 #total number of codons in the cds
    for c in list(codon_count.keys()):    #O(len codondict)
        codon_count[c]/=total_cod
    
    df_codcnt=pd.DataFrame(list(codon_count.items()) )
    df_codcnt.columns=['Codon', 'Fractional_Freq']
    df_codcnt=df_codcnt.set_index('Codon').T.reset_index(drop=True)
    
    df_codcnt['GeneID']=geneID
	#reorder columns
    cols2=[df_codcnt.columns[-1]]+sorted(df_codcnt.columns[:61])
    df_codcnt=df_codcnt[cols2]
    return df_codcnt	

if __name__=='__main__':
    #print a welcome message with my contacts
    msg = "Thanks for using CodonCount@2021. \n To send questions and suggestions, you may make a pull request at github OR  email rwint@ucmerced.edu"
    print(msg)
    parser = argparse.ArgumentParser(description="Computes the Fractional Codon Frequency for each coding sequence and outputs a csv table of results. Written by Rhondene Wint (rwint@ucmerced.edu).")
    parser.add_argument('filename', help='Name or path of input fasta file with coding sequences. E.g. CDS.fasta')
    parser.add_argument('out_name', help='name or path of output file. E.g. file_codcount')
    args=parser.parse_args()
	
    transcriptome=preproc(args.filename)
    comb=[]  #store the codon content for each gene
    for i in range(transcriptome.shape[0]):
        gene = transcriptome.iloc[[i]]
        df_count=get_cod_freq(gene)
        comb.append(df_count)
    pd.concat(comb,axis=0).to_csv(args.out_name+'.csv',index=False)

    
		
	
	
   
	
    

    