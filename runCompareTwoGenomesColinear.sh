#!/bin/bash
#RLW 2017,2019
if [ $# -ne 11 ]; then
        echo "Usage: $(basename $0)"
        echo "  QUERY FASTA"
        echo "  REFERENCE/SUBJECT/TARGET FASTA"
        echo "  ALPHA TRANSPARENCY 0-255"
        echo "  MISMATCH THRESHOLD"
        echo "  BLOCK LENGTH (bp)"
        echo "  LEAP LENGTH (bp)"
        echo "  SCALE (1:n)"
        echo "  QUERY features GFF .tsv"
        echo "  REFERENCE features GFF .tsv"
        echo "  cross_match/minimap2"
        echo "  PATH-TO-FONTS"
        exit 1
fi

if [ ${11} == '.' ]; then
   if [ -z ${XM_FONTS+x} ]; then
      echo "WARN: No font path defined and ENV variable XM_FONTS not found"
   fi
else
   XM_FONTS=${11}
fi

echo Running: $(basename $0) $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${XM_FONTS}

FEATURE_OPTS=""
if [ $8 != '.' ]; then
   FEATURE_OPTS=" -y $8"
fi
if [ $9 != '.' ]; then
   FEATURE_OPTS=" -y $9"
fi

if [ ${10} == 'cross_match' ]; then
   if ! command -v cross_match &> /dev/null; then
      echo "ERROR: cross_match not found on path - for docker/singularity mount linux binary to /opt/cross_match/bin"
   fi
   # cross_match pipeline
   cross_match $1 $2 -minmatch 29 -minscore 59 -masklevel 101 > $1_vs_$2.rep
   xmatchview.py -x $1_vs_$2.rep -q $1 -s $2 -a $3 -m $4 -b $5 -r $6 -c $7 -f png $FEATURE_OPTS -p ${XM_FONTS}

elif [ ${10} == 'minimap2' ]; then
   if ! command -v minimap2 &> /dev/null; then
      echo "ERROR: minimap2 not found on path"
   fi
   # minimap pipeline
   minimap2 $2 $1 -N200 -p0.0001 > $1_vs_$2.paf
   xmatchview.py -x $1_vs_$2.paf -q $1 -s $2 -a $3 -m $4 -b $5 -r $6 -c $7 -f png $FEATURE_OPTS -p ${XM_FONTS}

else

   echo Unrecognizable option ${10}
   echo Make sure you specify: cross_match OR minimap2

fi

