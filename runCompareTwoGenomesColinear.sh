#!/bin/bash
#Rene Warren 2017
#example bash script for running cross_match and xmatchview 

if [ $# -ne 6 ]; then
        echo "Usage: $(basename $0) <QUERY FASTA> <REFERENCE FASTA> <ALPHA TRANSPARENCY 0-255> <MISMATCH THRESHOLD> <QUERYfeatures.tsv> <REFfeatures.tsv>"
        exit 1
fi

source /home/dmacmillan/anaconda2/envs/legacy/bin/activate
/home/pubseq/BioSw/phrap/current/cross_match $2 $1 -minmatch 29 -minscore 59 -masklevel 101 > $1_vs_$2.rep
/home/dmacmillan/anaconda2/envs/legacy/bin/python xmatchview.py -x $1_vs_$2.rep -s $2 -q $1 -m $4 -r 10 -l 1 -c 2 -a $3 -f png -y $5 -e $6
/home/pubseq/BioSw/phrap/current/cross_match $1 $2 -minmatch 29 -minscore 59 -masklevel 101 > $2_vs_$1.rep
/home/dmacmillan/anaconda2/envs/legacy/bin/python xmatchview.py -x $2_vs_$1.rep -s $1 -q $2 -m $4 -r 10 -l 1 -c 2 -a $3 -f png -y $6 -e $5
