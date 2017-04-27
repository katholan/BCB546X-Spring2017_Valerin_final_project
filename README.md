# BCB546X Final project
##Team Valerin: Valeria Velásquez, Katerina Holan, and Devin Molnau
#### Research article: Harkess, A., Mercati, F., Shan, H. Y., Sunseri, F., Falavigna, A., & Leebens-Mack, J. (2015). Sex-biased gene expression in dioecious garden asparagus (Asparagus officinalis). New Phytologist, 207(3), 883–892. http://doi.org/10.1111/nph.13389

This Github repository contains all the analysis that were made for the final project of the BCB546X course. The repository consists of 4 major directories, the final presentation, and this README file. The structure of directories are as follows: 

#### Paper_supplementary _data folder:
This directory contains a README file, the original paper, and the supplementary tables and materials.

####Trimming_Analysis folder: 

* The first part of the transcriptome profiling analysis, which is comprised of the quality assessment of the raw and trimmed reads using [FASTQC](https://wiki.hpcc.msu.edu/display/Bioinfo/FastQC+Tutorial), the trimming using [cutadapt](https://github.com/marcelm/cutadapt) and the [MULTIQC](http://multiqc.info/) for summarizing the results.

* This folder also contains one README.md file with the script for the pipeline and the output FASTQC and MULTIQC files for visualization. 

#### Alignment_Counting folder: 

* This folder contains the second step in the transcriptome profiling analysis. The alignment and expression quantification using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [RSEM](https://github.com/deweylab/RSEM). 

* The folder contains one README.md file with the script and the output counting files, summaries of the RSEM statistical models in pdf.

*Note: Since the transcriptome profiling is highly intensive at the computational level, it was done for only 2 samples as an example*
  
#### Expression_Analysis folder: 

* This folder contains all the analysis for the differential expression, using R, specially the [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html) package. The counting data for all the samples was retrieved by contacting the paper's authors.
* The folder contains one README file with the description of all the md notebooks that were made for filtering the genes and making each group of figures, the input files, and all the md notebooks for making each figure.

##### Figures folder:

* This subfolder of the Expression Analysis contains the figures our group generated by expression analysis. 



*Valeria Velásquez was in charge of doing the transcriptome profiling analysis, and also took part in the edgeR gene filtering procedure to get the candidate genes for the differential expression analysis*

*Devin Molnau was in charge of generating figures 1,2, and 3 in the expression analysis using the R package edgeR.*

*Katerina Holan was in charge of generating figure 4 and formatting the figures. Additionally, she made the comparison between the candidate gene list obtained by the authors and by our analysis as a quality control measurement for our pipeline.*
