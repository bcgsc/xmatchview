#!/bin/bash
#RLW 2017,2019
if [ $# -ne 11 ]; then
        echo "Usage: $(basename $0)"
        echo " <QUERY FASTA>"
        echo " <REFERENCE/TARGET FASTA>"
        echo " <ALPHA TRANSPARENCY 0-255>"
        echo " <MISMATCH THRESHOLD>"
        echo " <BLOCK LENGTH (bp)>"
        echo " <LEAP LENGTH (bp)>"
        echo " <SCALE (1:n)>"
        echo " <QUERY features GFF .tsv>:"
        echo " <REFERENCE features GFF .tsv>"
        echo " <cross_match/minimap2>"
        echo " <PATH-TO-FONTS>"
        exit 1
fi

echo Running: $(basename $0) $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11}

# source PATH-TO-SOURCE (IF NEEDED)

echo "Make sure xmatchview.py, cross_match and minimap2 are in your PATH"

if [ ${10} == 'cross_match' ]; then

   # cross_match pipeline
   cross_match $1 $2 -minmatch 29 -minscore 59 -masklevel 101 > $1_vs_$2.rep
   xmatchview.py -x $1_vs_$2.rep -q $1 -s $2 -a $3 -m $4 -b $5 -r $6 -c $7 -f png -y $8 -e $9 -p ${11}

elif [ ${10} == 'minimap2' ]; then

   # minimap pipeline
   minimap2 $2 $1 -N200 -p0.0001 > $1_vs_$2.paf
   xmatchview.py -x $1_vs_$2.paf -q $1 -s $2 -a $3 -m $4 -b $5 -r $6 -c $7 -f png -y $8 -e $9 -p ${11}

else

   echo Unrecognizable option ${10}
   echo Make sure you specify: cross_match OR minimap2

fi

