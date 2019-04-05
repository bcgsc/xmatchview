![Logo](https://github.com/warrenlr/xmatchview/blob/master/xmv-logo.png)

# xmatchview
## Genome alignment visualization
## xmatchview Rene L. Warren, 2005-2019
## email: rwarren [at] bcgsc [dot] ca
## Visit www.bcgsc.ca/bioinfo/software/xmatchview for additional information

### NAME
   <pre>
   xmatchview.py v1.1.1   December 2018
   xmatchview.py v1.1   October 2018
   xmatchview.py v1.0   January 2018 - Post JOSS review
   xmatchview.py v0.3.3 January 2018
   xmatchview.py v0.3   November 2017
   XMatchView.py v0.2   March 2005/May 2005/January 2006
   </pre>

### SYNOPSIS
   xmatchview and xmatchview-conifer are imaging tools for comparing the synteny between DNA sequences. It allows users to align 2 DNA sequences in fasta format using cross_match and displays the alignment in a variety of image formats.
   xmatchview and xmatchview-conifer are written in python and run on linux and windows. They serve as visual tools for analyzing cross_match alignments. Cross_match (Green, P. (1994) http://www.phrap.org) uses an implementation of the Smith-Waterman algorithm for comparing DNA sequences that is sensitive.

### LICENSE PREAMBLE
   Copyright (c) 2005-2018 Rene Warren, Canada's Michael Smith Genome Science Centre.  All rights reserved.
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

Download the tar file and extract the files on your system using:

<pre>
tar -xvf xmatchview_1-1-1.tar 
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
Usage: ['xmatchview.py'] v1.1.1
-x alignment file (cross_match .rep or Pairwise mApping Format .paf)
-s reference genome fasta file (sequence on a single line)
-q query contig/genome fasta file (sequence on a single line)
-e reference features (eg. exons) coordinates tsv file (start end *) - optional
-y query features (eg. exons) coordinates tsv file (start end *) - optional
-m mismatch threshold (e.g. -m 10 allows representation of repeats having up to 10% mismatch
-b length (bp) of similarity block to display
-c scale (pixel to basepair scale, for displaying the image)
-r leap length (bp) to evaluate repeat frequency (smaller numbers will increase the resolution, but will affect drastically the run time.  recommended -r=50)
-a alpha value, from 0 (transparent) to 255 (solid, default)
-f output image file format (png, tiff, jpeg, or gif) NOTE: the png and tiff are better.
-p full path to the directory with fonts on your system (please refer to the documentation for fonts used)


Usage: ['xmatchview-conifer.py'] v1.1.1
-x alignment file (cross_match .rep or Pairwise mApping Format .paf)
-s reference genome fasta file (sequence on a single line)
-q query contig/genome fasta file (sequence on a single line)
-e reference features (eg. exons) coordinates tsv file (start end *) - optional
-y query features (eg. exons) coordinates tsv file (start end *) - optional
-m maximum mismatch threshold (e.g. -m 10 allows representation of repeats having up to 10% mismatch
-b minimum length (bp) of similarity block to display
-c scale (pixel to basepair scale, for displaying the image)
-r leap length (bp) to evaluate repeat frequency (smaller numbers will increase the resolution, but will affect drastically the run time.  recommended -r=50)
-l label for the tree trunk (6 characters or less for best result)
-a alpha value, from 0 (transparent) to 255 (solid, default)
-f output image file format (png, tiff, jpeg, or gif) NOTE: the png and tiff are better.
-p full path to the directory with fonts on your system (please refer to the documentation for fonts used)

Note: Files for the -s and -q options must correspond to fasta files used to run cross_match.  Those files should each contain a single, non-justified, sequence (no line breaks).

fasta file format:

>yourSequence
AATAGCAGCTACGACGACGCAGCGCGACGTTTCATCAA...AATACAGACGCGACGACGCAGCATCATCGAGAC


* A third column may be used to specify the color of a feature (default feature color is yellow or black, for xmatchview and xmatchview-conifer, respectively). Users may specify any of these color names: yellow, blue, cyan, green, lime, red, sarin, forest, dirtyred, dirtyyellow, grey, lightgrey, orange, beige, black, white.

</pre>
Users can control whether to show the position of sequence features on the reference and query (*-e* and *-y* options), show co-linear blocks of a certain length (*-b* option) when their mismatch rates are below a threshold (*-m* option). The histogram is generated by moving a sliding window with a step length (*-r* recommended between 10-50). The color space in xmatchview is RGBA and the alpha channel is used for visualizing the relationship between co-linear blocks (*-a* option, transparent to solid, 0 to 255).



### RUNNING THE cross_match/xmatchview/xmatchview-conifer PIPELINES
------------- 
Demo shell scripts that pipeline cross_match and xmatchview/xmatchview-conifer are included with the distribution and provided for guidance (runCompareTwoGenomesColinear.sh and runSpruceView.sh, respectively).

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
To test your xmatchview install on your system, we provide a test folder where you can run xmatchview and xmatchview-conifer, using the shell script commands below. If all goes well, and you used the arial fonts provided, the image you generate should be identical to:
xmv-conifer_FTL1_ss.fa_vs_FTL1_pa.fa.rep_m10_b10_r1_c2_success.png
and
xmv_FTL1_ss.fa_vs_FTL1_pa.fa.rep_m10_r10_l1_c1_success.png 

Once you confirmed that it works as expected, you may explore the full range of parameters to test functionality (see USAGE above) 


At the unix prompt, once the package is installed, change directory to:
<pre>
cd ./test
</pre>
Once you have downloaded pyhon and PIL and changed the paths to fonts in the xmatchview.py and xmatchview-conifer.py
Execute:
<pre>
./runXMV-conifer.sh FTL1_ss.fa_vs_FTL1_pa.fa.rep FTL1_ss.fa FTL1_pa.fa 200 10 1 FTL1_ss.txt FTL1_pa.txt
and
./runXMV.sh FTL1_ss.fa_vs_FTL1_pa.fa.rep FTL1_ss.fa FTL1_pa.fa 200 10 1 FTL1_ss.txt FTL1_pa.txt
</pre>
If all went well, images such as those provided in the test folder should be generated


### CITING xmatchview/xmatchview-conifer
-------------

Thank you for using, developing and promoting this free software.
If you use xmatchview/xmatchview-conifer for you research, please cite:

<pre>
Warren, RL (2018). Visualizing genome synteny with xmatchview. Journal of Open Source Software, 3(21):497.
</pre>

[![DOI](http://joss.theoj.org/papers/10.21105/joss.00497/status.svg)](https://doi.org/10.21105/joss.00497)


### WHAT'S NEW in v1.1.1
------------------
<pre>
-Bug fixes (will now return an error when the alignment files are empty [instead of plotting empty graphs])
</pre>

### WHAT'S NEW in v1.1
------------------
<pre>
-Bug fixes (the forward synteny blocks were always printed, regardless of specified filters. This is fixed in both xmatchview and xmatchview-conifer v1.1)
</pre>

### WHAT'S NEW in v1.0
------------------
<pre>
-Published (JOSS peer-review) version
</pre>

### WHAT'S NEW in v0.3.3
------------------
<pre>
-Initial support for .paf (Pairwise mApping Format) alignment files.
</pre>

### WHAT'S NEW in v0.3.2
------------------
<pre>
-Made options consistent between xmatchview.py and xmatchview-conifer.py
-Included new option to specify font path (-p)
-Added option to specify feature colors in the tsv files provided with the -e and -y options. A third column may be used to specify the color of a feature (default feature color is yellow or black, for xmatchview and xmatchview-conifer, respectively). Users may specify any of these color names: yellow, blue, cyan, green, lime, red, sarin, forest, dirtyred, dirtyyellow, grey, lightgrey, orange, beige, black, white.
</pre>

### WHAT'S NEW in v0.3
------------------
<pre>
-Plot colinear blocks and sequence relationships with transparent color (alpha, supplied with -a)
-Plot the position of exons on the reference and query DNA segments (-e and -y arguments, optional)
-Plot the position of Ns in query and reference sequences
-Bug fixes
</pre>
---
Please find the v0.2 release in the corresponding subdirectory

Questions/Comments: Rene Warren : rwarren at bcgsc dot ca
