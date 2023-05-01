#!/bin/bash

#An example of using annotate peaks to create a history of relative location of motifs to peak. This uses the files with all motifs, 
#but you can use this for a specific motif by extracting from the discovery results or creating your own motif file as described on Homer's website
#The output can be graphed in excel
annotatepeaks.pl ./peakcalls/TS_01.bam_summits.bed Fence_Lizard -m ./TS_01.bam_summits/1homerMotifs.all.motifs -size 1000 -hist 10 > TS_01_annotatepeaks_hist.txt