#! usr/bin/python
# fasta_to_csv
##make sequence table

import sys
import argparse
        
import pandas as pd


def make_csv(fasta,out_name):
    headers=[]
    seqs=[]
    with open(fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                headers.append(line.rstrip())
            else:
                seqs.append(line.rstrip())
        assert len(headers)==len(seqs),' Number of IDs and sequences do not match'
    
                
    df = pd.DataFrame(zip(headers,seqs))
    df.columns=['Header','Sequence']
    df.to_csv(out_name+'.csv',index=False)

if __name__=='__main__':

    parser = argparse.ArgumentParser(description="Converts a fasta file to a tabular comma separated value (csv) file. Written by Rhondene Wint (rwint@ucmerced.edu). Here is an example for running this sript: 'python fasta2csv.py myfile.fasta myfile.csv'",
                                     epilog='If an AssertionError was raised, try running my fix_fasta.py script (on my github) which outputs a corrected fasta file. Then re-run fasta2csv.py using the corrected fasta file')
    parser.add_argument('fasta_file', help='Input path of fasta file')
    parser.add_argument('out_name', help='Path for the output csv file')

    args = parser.parse_args()

    make_csv(args.fasta_file, args.out_name)
    
          
                  
