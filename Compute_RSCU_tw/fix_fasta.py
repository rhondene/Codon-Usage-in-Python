#!/usr/bin/python
# Author: Rhondene Wint
#Useful for correcting fasta file with newlines separating nucleotides of the same sequence.

import sys
import argparse



def fix_fasta(filename):
    headers = []
    sequences=[]
    file=[]
    with open(filename, 'r') as f:
        for line in f:
            file.append(line.strip())

        for i in range(len(file)):
            line = file[i]
            if line.startswith('>') is True:
                header =(line.rstrip())
                seq = ""
                j=i  #the index for the sequence
                while True:
                    j+=1  #the sequence is one position right of the header line
                    if j==len(file):
                        break
                    #when it encounters another header
                    if file[j].startswith('>') is True:
                        break  
                    else:
                        s = file[j]
                        seq+=s
                headers.append(header)
                sequences.append(seq)
            else:
                continue

    assert len(headers)==len(sequences), 'Number of headers does not equal number of sequences'
    return [headers,sequences]
	
	
#--------------------------
if __name__ =='__main__':
    
    parser = argparse.ArgumentParser(description="Corrects within-sequence newlines in fasta (.fasta, .fa, .fas) file. Example of using this script: \n\n 'python fix_fasta.py myfile.fasta myfile_fixed.fasta' ", epilog='Enjoy!')
    parser.add_argument('filename', help=' input path for the fasta file to be corrected. If fast file is in current directory then only the filename is needed', type=str)
    parser.add_argument('new_filename', help='output path for the corrected version of the fasta file', type=str)
    
    args = parser.parse_args()

    fasta=fix_fasta(args.filename)  
    
    ### re-rewrite the file as a proper fasta such that only newline between sequences and headers
    with open(args.new_filename, 'w') as f:
        for header,seq in fasta:
            f.write(header+'\n'+seq+'\n')

    

