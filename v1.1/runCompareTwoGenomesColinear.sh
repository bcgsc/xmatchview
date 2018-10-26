#!/bin/bash
#RLW 2017
if [ $# -ne 8 ]; then
        echo "Usage: $(basename $0) <QUERY FASTA> <REFERENCE FASTA> <ALPHA TRANSPARENCY 0-255> <MISMATCH THRESHOLD> <SCALE> <QUERYfeatures.tsv> <REFERENCEfeatures.tsv> <PATH-TO-FONTS>"
        exit 1
fi

# source PATH-TO-SOURCE (IF NEEDED)
cross_match $2 $1 -minmatch 29 -minscore 59 -masklevel 101 > $1_vs_$2.rep
python xmatchview.py -x $1_vs_$2.rep -s $2 -q $1 -m $4 -b 10 -r 10 -c $5 -a $3 -f png -y $6 -e $7 -p $8
cross_match $1 $2 -minmatch 29 -minscore 59 -masklevel 101 > $2_vs_$1.rep
python xmatchview.py -x $2_vs_$1.rep -s $1 -q $2 -m $4 -b 10 -r 10 -c $5 -a $3 -f png -y $7 -e $6 -p $8
