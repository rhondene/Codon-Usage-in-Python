#!/usr/bin/python
# Useful for correcting fasta file with newlines separating nucleotides of the same sequence

import sys
import argparse



def fix_fasta(filename, new_filename):
    fasta = []
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
                    if file[j].startswith('>') is True:
                        break  #when it encounters another header
                    else:
                        s = file[j]
                        seq+=s
                fasta.append(header)
                fasta.append(seq)
            else:
                continue

    ### re-rewrite the file as a proper fasta such that only newline between sequences and headers
    with open(new_filename, 'w') as f:
        for line in fasta:
            f.write(line+'\n')
    return 
	
	
#--------------------------
if __name__ =='__main__':
    
    parser = argparse.ArgumentParser(description="Corrects within-sequence newlines in fasta (.fasta, .fa, .fas) file. Example of using this script: \n\n 'python fix_fasta.py myfile.fasta myfile_fixed.fasta' ", epilog='Enjoy!')
    parser.add_argument('filename', help=' input path for the fasta file to be corrected. If fast file is in current directory then only the filename is needed', type=str)
    parser.add_argument('new_filename', help='output path for the corrected version of the fasta file', type=str)
    
    args = parser.parse_args()

    fix_fasta(args.filename,args.new_filename)

    

