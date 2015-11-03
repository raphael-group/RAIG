Copyright 2014 Brown University, Providence, RI.

                         All Rights Reserved

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose other than its incorporation into a
commercial product is hereby granted without fee, provided that the
above copyright notice appear in all copies and that both that
copyright notice and this permission notice appear in supporting
documentation, and that the name of Brown University not be used in
advertising or publicity pertaining to distribution of the software
without specific, written prior permission.

BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
http://compbio.cs.brown.edu/software

README for RAIG (Recurrent Aberrations from Interval Graph)

Software that reconstructs a cancer genome as a rearrangement of segments,
or intervals, from the reference genome using paired end sequencing data.

If you use this software in your research, please cite:

H. Wu, I. Hajirasouliha, B.J. Raphael (2014)
Detecting independent and recurrent copy number aberrations using interval 
graphs. Bioinformatics (Proceedings of 22nd Annual International 
Conference on Intelligent Systems in Molecular Biology (ISMB)).
Vol. 30 ISMB2014, pages i195-i203.

Corresponding Authors:
Hsin-Ta Wu (bournewu@cs.brown.edu)
Ben Raphael (braphael@cs.brown.edu)

Version: 1.02
Version Date: Nov 3, 2015

SUMMARY ==========================================================

This README file lists and annotates the files of the software 
for the RAIG algorithm.

CONTENTS ========================================================

(i) Software (in src/ subdirectory ):
* Source code in src

(ii) Genome information (in required_files/ subdirectory ):
* Required files for RAIG, including gene database, cytoband database, and probe/marker information.

(iiii) Input (in input/ subdirectory ):
* Example input files: example.input.gbm.seg.txt
