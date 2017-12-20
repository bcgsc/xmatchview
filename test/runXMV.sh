#!/bin/bash
#Rene Warren 2017

source /home/dmacmillan/anaconda2/envs/legacy/bin/activate
/home/dmacmillan/anaconda2/envs/legacy/bin/python xmatchview.py -x FTL1_ss.fa_vs_FTL1_pa.fa.rep -s FTL1_pa.fa -q FTL1_ss.fa -m 10 -r 10 -l 1 -c 1 -a 200 -f png -y FTL1_ss.txt -e FTL1_pa.txt
