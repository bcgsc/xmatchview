# xmatchview
## Genome alignment visualization
## xmatchview v0.3 Rene L. Warren, 2005-2017
## email: rwarren [at] bcgsc [dot] ca
## Visit www.bcgsc.ca/bioinfo/software/xmatchview for additional information

### NAME
   xmatchview.py v0.3  November 2017
   XMatchView.py v0.2  March 2005/May 2005/January 2006

### SYNOPSIS
   Allows users to align 2 DNA sequences in fasta format using cross_match and displays the alignment in a variety of image formats
   xmatchview is written in python and runs on linux and windows. It serves as a visual tool for analyzing cross_match alignments. Cross_match (Green, P. (1994) http://www.phrap.org) uses an implementation of the Smith-Waterman algorithm for comparing DNA sequences that is sensitive.

### LICENSE
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


### INSTALL
<pre>
   A) you will need to do the following before you can proceed:
      1)Download python2.3 or 2.4 from: http://www.python.org/ and change the shebang line to reflect this
      2)Download the Python Imaging Library (PIL) from: http://www.pythonware.com/products/pil/
      3)Copy true type fonts from c:\WINDOWS\Fonts to a unix directory and change the line truetype= below to reflect the location of your ttf
      4)Change the sys.path.append line below to reflect the location of PIL
      5)Make sure cross_match is in your $PATH or change the line cross_match_exec= below
      6)Copy the image pbp.gif to the same directory where the XMatchView.py program resides, make a fake gif with that name or comment the whole "###Just for fun code block" below

   B) If you're running this program remotely, but on the GSC servers make sure you are running it on xhost01.bcgsc.ca
   C) You can run xmatchview.py on Windows XP, provided that you have installed python and PIL and that you changed the script line that specifies the location of the fonts to reflect their location in your windows computer. However, you won't be able to run crossmatch with it, unless you have obtained cross_match for windows.
</pre>

### RUN

A non-GUI version of xmatchview (xmatchview.py) exists and is included with this release.

<pre>
Usage: ['xmatchview.py'] v0.3
-x crossmatch file
-s reference genome fasta file
-q query contig/genome fasta file
-e reference exon coordinates tsv file (start end) - optional
-y query exon coordinates tsv file (start end) - optional
-m mismatch threshold (e.g. -m 10 allows representation of repeats having up to 10% mismatch
-r length of similarity block to display
-c scale (for displaying the image)
-l leap to evaluate repeat frequency (smaller numbers will increase the resolution, but will affect drastically the run time.  recommended -l=50)
-f file format (bmp, jpeg, png, ps, gif, pdf, tiff) NOTE: the png, ps, tiff and bmp are much better.
-a alpha value, from 0 (transparent) to 255 (solid, default)
-p transform bacterial ORF into protein (i.e. plot alignment between ORF products? 1/0)
* Files for the -s and -q options must correspond to fasta files used to run cross_match
</pre>

### What's new in v0.3
------------------

Plot colinear blocks and sequence relationships with transparent color (alpha, supplied with -a)
Plot the position of exons (yellow half rectangle) on the reference and query DNA segments (-e and -y arguments, optional)
Plot the position of Ns (red ticks) in query and reference sequences 
Bug fixes

---
Please find the v0.2 release in the corresponding subdirectory

Questions/Comments: Rene Warren : rwarren at bcgsc dot ca
