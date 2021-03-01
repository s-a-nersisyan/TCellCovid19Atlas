#!/bin/bash

peptide=$1
echo $peptide > $peptide.pep
rm ../$peptide.txt
len=$(wc -l ../HLA_table.tsv | awk '{print $1}')
let ind=0
while read arg; do
    let ind++
    echo "$arg $ind $len"
    netMHCpan -BA -p $peptide.pep -a $arg >> $peptide.aff
done < "../HLA_table.tsv"
# cat ../HLA_table.tsv | parallel -j 48 "netMHCpan -BA -p $peptide.pep -a {} >> ../$peptide.txt"
python3 agregate_affinities.py $peptide
rm $peptide.aff
rm $peptide.pep
