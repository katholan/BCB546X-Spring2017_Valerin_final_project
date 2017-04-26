# BCB546X Final project
##Team Valerin: Valeria Velásquez, Katerina Holan, and Devin Molnau
#### Research article: Harkess, A., Mercati, F., Shan, H. Y., Sunseri, F., Falavigna, A., & Leebens-Mack, J. (2015). Sex-biased gene expression in dioecious garden asparagus (Asparagus officinalis). New Phytologist, 207(3), 883–892. http://doi.org/10.1111/nph.13389

This Github repository contains all the analysis that were made for the final project of the BCB546X course. The structure of directories is split in 4: 

#### First, we have original paper and the supplementary tables in the Paper_supplementary _data folder and the final presentation.

#### Second, the Trimming_Analysis folder: 

* It is the first part of the transcriptome profiling analysis, which comprisses the qualitity assesment of the raw and trimmed reads using [FASTQC](https://wiki.hpcc.msu.edu/display/Bioinfo/FastQC+Tutorial), the trimming using [cutadapt](https://github.com/marcelm/cutadapt) and the [MULTIQC](http://multiqc.info/) for summarizing the results.

* The folder  contains one README.md file with the script for the pipeline and the output FASTQC and MULTIQC files for visualization. 

#### Third, the Alignment_Counting folder: 

* This folder contains the second step in the transcriptome profiling analysis. The alignment and expression quantification using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [RSEM](https://github.com/deweylab/RSEM). 

* The folder contains one README.md file with the script and the output counting files, summaries of the RSEM statistical models in pdf.

*Note: Since the transcriptome profiling is highly intensive at the computationl level, it was done for only 2 samples as an example*
  
#### Expression_Analysis folder: 

* This folder contains all the analysis for the differential expression, using R, specially the [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html) package. The counting data for all the samples was retrieved by contacting the paper's authors.
* The folder contains one README file with the description of all the md notebooks that were made for filtering the genes and making each group of figures, the input files, and all the md notebooks for making each figure.

#### Figures folder

* This folder contains the figures our group made for our analysis and the presentation.
	* 



*Valeria Velásquez was in charge of doing the transcriptome profiling analysis, and also took part in the edgeR gene filtering procedure to get the candidate genes for the differential expression analysis*

*Devin Molnau also participated in the edgeR gene filtering, including the testing of the experimental design for the pairwise comparisons, and also took part in making the figures 1, 2 and 3*

*Katerina Holan was involved in formatting the figures. Additionally, she made the comparison between the candidate genelist obtained by the authors and by our analysis as a quality control measurment for our pipeline. She also made the figure 4 by taking the normalized TMM matrix from the edgeR output*
