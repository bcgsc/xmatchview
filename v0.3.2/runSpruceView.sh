#!/bin/bash
#RLW 2017
if [ $# -ne 8 ]; then
        echo "Usage: $(basename $0) <QUERY FASTA> <REFERENCE FASTA> <LABEL> <ALPHA TRANSPARENCY 0-255> <MISMATCH THRESHOLD> <QUERYfeatures.tsv> <REFERENCEfeatures.tsv> <PATH-TO-FONTS>"
        exit 1
fi

cross_match $2 $1 -minmatch 5 -minscore 10 -masklevel 101 > $1_vs_$2.rep
python xmatchview-conifer.py -x $1_vs_$2.rep -s $2 -q $1 -m $5 -b 10 -r 1 -c 2 -l $3 -a $4 -f png -y $6 -e $7 -p $8
cross_match $1 $2 -minmatch 5 -minscore 10 -masklevel 101 > $2_vs_$1.rep
python xmatchview-conifer.py -x $2_vs_$1.rep -s $1 -q $2 -m $5 -b 10 -r 1 -c 2 -l $3 -a $4 -f png -y $7 -e $6 -p $8
