###Expression Analysis Folder README file

This file contains the expression analysis scripts by Katerina Holan and Devin Molnau.

Devin Molnau used R version 3.3.1 and Katerina used 3.2.1

Expression\_Analysis\_Markdown.Rmd is an Rmarkdown file that contains their combined expression analysis for figures 1-4 of the study.

Asparagus.RSEM_genecounts.txt is the raw counts data from the paper and is used as an input for the expression analysis. 

DEG_counts.csv is the counts data of just the experimentally expressed genes that we found in or expression analysis.

The Df\_gene\_names directory contains 3 text files with the names of the differentially expressed genes from the three pairwise comparisons, pulled from supplementary table 1.

The figures directory contains all figures generated during analysis.

TMM\_normalized\_FPKM\_matrix.txt is a file used during heatmap clustering.

SRR1639681\_counting.genes.results and SRR1642915\_counting.genes.results are the 2 samples that were run in the assembly pipeline that are used for generating the FPKM heatmap.
