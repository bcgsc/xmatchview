#!/bin/bash
#Rene Warren 2017
# COMMAND TO RUN:
# ./runXMV-conifer.sh FTL1_ss.fa_vs_FTL1_pa.fa.rep FTL1_ss.fa FTL1_pa.fa 200 10 1 FTL1_ss.txt FTL1_pa.txt

if [ $# -ne 8 ]; then
        echo "Usage: $(basename $0) <CROSS_MATCH .rep> <QUERY FASTA .fa> <REFERENCE FASTA .fa> <ALPHA TRANSPARENCY 0-255> <MISMATCH THRESHOLD> <SCALE> <QUERYfeatures.tsv> <REFERENCEfeatures.tsv>"
        exit 1
fi

# source PATH-TO-SOURCE (IF NEEDED)
python xmatchview-conifer.py -x $1 -s $3 -q $2 -a $4 -m $5 -b 10 -r 2 -l FTL1 -c $6 -f png -y $7 -e $8
