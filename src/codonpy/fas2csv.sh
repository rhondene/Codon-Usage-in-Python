#!/bin/bash

for file in $(ls *fixed*)
do
python ../fasta2csv.py $file ${file:0:5}
done
