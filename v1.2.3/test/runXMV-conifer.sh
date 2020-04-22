#!/bin/bash
# Rene Warren 2017, 2019
# COMMAND TO RUN:
# ./runXMV-conifer.sh FTL1_pa.fa_vs_FTL1_ss.fa.rep FTL1_pa.fa FTL1_ss.fa 200 10 2 FTL1_pa.gff FTL1_ss.gff 
# ./runXMV-conifer.sh FTL1_paMF.fa_vs_FTL1_ssMF.fa.rep FTL1_paMF.fa FTL1_ssMF.fa 200 10 2 FTL1_pa.gff FTL1_ss.gff
# ./runXMV-conifer.sh FTL1_paMF.fa_vs_FTL1_ssMF.fa.paf FTL1_paMF.fa FTL1_ssMF.fa 200 99 2 FTL1_pa.gff FTL1_ss.gff 

if [ $# -ne 8 ]; then
        echo "Usage: $(basename $0) <CROSS_MATCH .rep> <QUERY FASTA .fa> <REFERENCE FASTA .fa> <ALPHA TRANSPARENCY 0-255> <MISMATCH THRESHOLD> <SCALE> <QUERYfeatures.tsv> <REFERENCEfeatures.tsv>"
        exit 1
fi

# source PATH-TO-SOURCE (IF NEEDED)
python3 ../xmatchview-conifer.py -x $1 -s $3 -q $2 -a $4 -m $5 -b 10 -l FTL1 -c $6 -f png -y $7 -e $8 -p ../../tarballs/fonts
