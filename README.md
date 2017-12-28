# xmatchview
## Genome alignment visualization
## xmatchview v0.3 Rene L. Warren, 2005-2017
## email: rwarren [at] bcgsc [dot] ca
## Visit www.bcgsc.ca/bioinfo/software/xmatchview for additional information

### NAME
   <pre>
   xmatchview.py v0.3  November 2017
   XMatchView.py v0.2  March 2005/May 2005/January 2006
   </pre>

### SYNOPSIS
   xmatchview and xmatchview-conifer are imaging tools for comparing the synteny between DNA sequences. It allows users to align 2 DNA sequences in fasta format using cross_match and displays the alignment in a variety of image formats.
   xmatchview and xmatchview-conifer are written in python and run on linux and windows. They serve as visual tools for analyzing cross_match alignments. Cross_match (Green, P. (1994) http://www.phrap.org) uses an implementation of the Smith-Waterman algorithm for comparing DNA sequences that is sensitive.

### LICENSE PREAMBLE
   Copyright (c) 2005-2017 Rene Warren, Canada's Michael Smith Genome Science Centre.  All rights reserved.
   xmatchview is a utility for comparing, visually, two DNA/RNA sequences

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.

### IMPLEMENTATION
-------------

xmatchview and xmatchview-conifer are implemented in PYTHON and run on any OS where PYTHON is installed.


### COMMUNITY GUIDELINES
-------------

I encourage the community to contribute to the development of this software, by providing suggestions for improving the code and/or directly contributing to the open source code for these tools. Users and developers may report software issues, bug fix requests, comments, etc, at <https://github.com/warrenlr/xmatchview>


### INSTALL

Download the tar ball, gunzip and extract the files on your system using:

<pre>
gunzip xmatchview-0.3.tar.gz 
tar -xvf xmatchview-0.3.tar
</pre>

### DEPENDENCIES
-------------

<pre>
You will need to do the following before you can proceed:
      1) Download python from: http://www.python.org/  (The code was developed and tested on python2.3 or 2.4 but may work with newer versions of python (not tested))
      2) Download the Python Imaging Library (PIL) from: http://www.pythonware.com/products/pil/
      3) Either use PIL or true type fonts (ttf) and specify the full path to the font directory in xmatchview.py and xmatchview-conifer.py with the -p option. The fonts used by xmatchview and xmatchview-conifer include:

        arial.ttf
        arialbd.ttf
        arialbi.ttf
        -or-
        helvR14.pil
        helvR18.pil
        helvB18.pil
        helvBO18.pil
        helvB24.pil
        helvR24.pil
        helvB24.pil

        If fontpaths are not provided (not recommended), default fonts will be used to preserve code functionality. However, these systems fonts are very small and of limited utility. To download the above fonts, simply search the internet for "arial.ttf", "arialbd.ttf" and "arialbi.ttf". For convenience, they are included in the "fonts.tar" file in this directory. The "helv*.pil" fonts are distributed with PIL.

      4) Download cross_match for academic use, see http://www.phrap.org and http://www.phrap.org/consed/consed.html#howToGet
      5) Make sure cross_match is in your $PATH
</pre>

### USAGE 
---------------
<pre>
Usage: ['xmatchview.py'] v0.3.2
-x crossmatch file
-s reference genome fasta file
-q query contig/genome fasta file
-e reference features (eg. exons) coordinates tsv file (start end) - optional
-y query features (eg. exons) coordinates tsv file (start end) - optional
-m mismatch threshold (e.g. -m 10 allows representation of repeats having up to 10% mismatch
-r length (bp) of similarity block to display
-c scale (pixel to basepair scale, for displaying the image)
-l leap (bp) to evaluate repeat frequency (smaller numbers will increase the resolution, but will affect drastically the run time.  recommended -l=50)
-f output image file format (png, tiff, jpeg, or gif) NOTE: the png and tiff are better.
-a alpha value, from 0 (transparent) to 255 (solid, default)
-p full path to the directory with fonts on your system (please refer to the documentation for fonts used)

Usage: ['xmatchview-conifer.py'] v0.1.1
-x crossmatch file
-s reference genome fasta file
-q query contig/genome fasta file
-e reference features (eg. exons) coordinates tsv file (start end) - optional
-y query features (eg. exons) coordinates tsv file (start end) - optional
-m maximum mismatch threshold (e.g. -m 10 allows representation of repeats having up to 10% mismatch
-b minimum length (bp) of similarity block to display
-c scale (pixel to basepair scale, for displaying the image)
-r basepair length leap to evaluate repeat frequency (smaller numbers will increase the resolution, but will affect drastically the run time.  recommended -r=50)
-l label for the tree trunk (6 characters or less for best result)
-a alpha value, from 0 (transparent) to 255 (solid, default)
-f output image file format (png, tiff, jpeg, or gif) NOTE: the png and tiff are better.
-p full path to the directory with fonts on your system (please refer to the documentation for fonts used)

Note: Files for the -s and -q options must correspond to fasta files used to run cross_match

* A third column may be used to specify the color of a feature (default feature color is yellow or black, for xmatchview and xmatchview-conifer, respectively). Users may specify any of these color names: yellow, blue, cyan, green, lime, red, sarin, forest, dirtyred, dirtyyellow, grey, lightgrey, orange, beige, black, white.

</pre>

### RUNNING THE cross_match/xmatchview/xmatchview-conifer PIPELINES
------------- 
Refer to:
<pre>
./runCompareTwoGenomesColinear.sh 

Usage: runCompareTwoGenomesColinear.sh <QUERY FASTA> <REFERENCE FASTA> <ALPHA TRANSPARENCY 0-255> <MISMATCH THRESHOLD> <SCALE> <QUERYfeatures.tsv> <REFERENCEfeatures.tsv> <PATH TO FONTS>

and:

./runSpruceView.sh 

Usage: runSpruceView.sh <QUERY FASTA> <REFERENCE FASTA> <LABEL> <ALPHA TRANSPARENCY 0-255> <MISMATCH THRESHOLD> <QUERYfeatures.tsv> <REFERENCEfeatures.tsv> <PATH TO FONTS>
</pre>

### TEST xmatchview.py / xmatchview-conifer.py
-------------

At the unix prompt, once the package is installed, change directory to 
<pre>
cd ./test
</pre>
Once you have downloaded pyhon and PIL and changed the paths to fonts in the xmatchview.py and xmatchview-conifer.py
Execute:
<pre>
./runXMV-conifer.sh 
and
./runXMV.sh
</pre>
If all went well, images such as those provided in the test folder should be generated


### CITING xmatchview/xmatchview-conifer
-------------

Thank you for using, developing and promoting this free software.
If you use xmatchview/xmatchview-conifer for you research, please cite:

<pre>
Warren RL. 2017. Visualizing genome synteny with xmatchview. bioRxiv 238220; doi: https://doi.org/10.1101/238220
</pre>


### WHAT'S NEW in v0.3
------------------
<pre>
-Plot colinear blocks and sequence relationships with transparent color (alpha, supplied with -a).
-Plot the position of exons on the reference and query DNA segments (-e and -y arguments, optional).
-Plot the position of Ns in query and reference sequences.
-Bug fixes.
</pre>
---
Please find the v0.2 release in the corresponding subdirectory

Questions/Comments: Rene Warren : rwarren at bcgsc dot ca
