#!/bin/bash
#RLW 2017
if [ $# -ne 7 ]; then
        echo "Usage: $(basename $0) <QUERY FASTA> <REFERENCE FASTA> <LABEL> <ALPHA TRANSPARENCY 0-255> <MISMATCH THRESHOLD> <QRY FEATURES> <REF FEATURES>"
        exit 1
fi

source /home/dmacmillan/anaconda2/envs/legacy/bin/activate
/home/pubseq/BioSw/phrap/current/cross_match $2 $1 -minmatch 29 -minscore 59 -masklevel 101 > $1_vs_$2.rep
/home/dmacmillan/anaconda2/envs/legacy/bin/python xmatchview-conifer.py -x $1_vs_$2.rep -s $2 -q $1 -m $5 -b 10 -r 1 -c 2 -l $3 -a $4 -f png -y $6 -e $7
/home/pubseq/BioSw/phrap/current/cross_match $1 $2 -minmatch 29 -minscore 59 -masklevel 101 > $2_vs_$1.rep
/home/dmacmillan/anaconda2/envs/legacy/bin/python xmatchview-conifer.py -x $2_vs_$1.rep -s $1 -q $2 -m $5 -b 10 -r 1 -c 2 -l $3 -a $4 -f png -y $7 -e $6
