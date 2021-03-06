# Alignment and expression estimation with RSEM

*Author: Valeria Velasquez Zapata*

This folder contains a serie the results of the attempt of read aligment and quantification of expression, which corresponds to the step after the read cleaning procedure, which is located in the FASTQC_files folder. For this analysis we will take the SRR1642915 sample adn then the SRR1639681 sample. 

The output of this pipeline is the counting data that is being used in the differential expression analysis in edgeR.
 
### Using bowtie2 to align clean reads to the transcriptome

To start we should get the transcriptome that the authors assembled *de novo*
For that purpose we can use the wget function

	$ wget http://datadryad.org/bitstream/handle/10255/dryad.90932/Asparagus_officinalis.trinity.fa.bz2?sequence=1

Then we can decompress the file and change its name

	$ mv Asparagus_officinalis.trinity.fa.bz2\?sequence\=1 asparagus_transcriptome.fa.bz2
	$ bzip2 -dk  asparagus_transcriptome.fa.bz2

Once bowtie2 is installed we can build the indexes for the alignment

	$ nohup bowtie2-build asparagus_transcriptome.fa asparagus_transcriptome &

### Using TOPHAT to map the reads

There are several modes to run Tophat, we sill try 2 of them: very sensitive and very fast (only with the SRR1642915 sample): 

	$ nohup tophat2 -o very_sensitive --b2-very-sensitive asparagus_transcriptome SRR1642915_Adapt_Quali_9bp_noShortReads_trimmed.fastq &
 
	$ nohup tophat2 -o very_fast --b2-very-sensitive asparagus_transcriptome SRR1642915_Adapt_Quali_9bp_noShortReads_trimmed.fastq &

### RSEM analysis for counting

For more detailed information about RSEM you can visit [https://github.com/bli25ucb/RSEM_tutorial](https://github.com/bli25ucb/RSEM_tutorial)

To install RSEM 

    $ wget https://github.com/deweylab/RSEM.git
    $ make install
 
#### Preparing reference

The first command that we will use will take the annotation file (supplementary information table S3 [http://onlinelibrary.wiley.com/doi/10.1111/nph.13389/full](http://onlinelibrary.wiley.com/doi/10.1111/nph.13389/full)

We formatted the file, by taking just the first two columns and adding "Aspof_" prefix to each entry. Once this was dome we are ready for building the reference for the counting analysis. We should build this reference based on the trancriptome file, since that was the way how the authors did the analysis (this step can also be done with the genome), by using the `--transcript-to-gene-map` function

    $ rsem-prepare-reference --transcript-to-gene-map ~/BCB546_final_project/data/rsem/asparagus_annotation.txt --bowtie2 --num-threads 10 ~/BCB546_final_project/data/bowtie_analysis/asparagus_transcriptome.fa asparagus_ref

Then we can do the counting, having the reference where the mapping and counting will be done:

    $ nohup rsem-calculate-expression -p 10 --bowtie2 ~/BCB546_final_project/data/trimmed_files/SRR1642915_Adapt_Quali_9bp_noShortReads_trimmed.fastq ~/BCB546_final_project/data/rsem/asparagus_ref SRR1642915_counting &

After finishing we can inpect the nohup.out file
    
    $ tail nohup.out
	17000000 alignment lines are loaded!
	18000000 alignment lines are loaded!
	19000000 alignment lines are loaded!
	20000000 alignment lines are loaded!
	21000000 alignment lines are loaded!
	Bam output file is generated!
	Time Used for EM.cpp : 0 h 04 m 05 s
	rm -rf SRR1642915_counting.temp

That way we confirm that the analysis have finished as expected.

We take the SRR1642915_counting.genes.results file as out counting file, especifically the expected_count column that indicates the counts per million (CPM).

    gene_id	transcript_id(s)length	effective_length	expected_count  TPM FPKM
    Aspof_comp100003_c0 	Aspof_comp100003_c0_seq1	728.00  540.63  0.000.000.00
    Aspof_comp100012_c0 	Aspof_comp100012_c0_seq1	1840.00 1652.63 0.000.000.00
    Aspof_comp100014_c0 	Aspof_comp100014_c0_seq1	619.00  431.63  0.000.000.00
    Aspof_comp100017_c0 	Aspof_comp100017_c0_seq1	569.00  381.63  0.000.000.00
    Aspof_comp100019_c0 	Aspof_comp100019_c0_seq1	569.00  381.63  0.000.000.00
    Aspof_comp100023_c0 	Aspof_comp100023_c0_seq1	374.00  186.63  0.000.000.00
   
If we want to get the statistics of the alignment process we can call 

    $ rsem-plot-model SRR1642915_counting SRR1642915_models.pdf

That will produce a pdf file with the lenght of the reads, a comparison of the observed quality and the phred score and a percentage of mapped reads

Only about 2% of reads were mapped. In addition, 46,7% of the reads are repeated, as we could see in the multiqc analysis. After looking at these results we realized that that was a paired-end data, so that can have an influence on the final results, as we treated it as single-end. We also made a different trimming for the file, which also could influence the results

Now we will take the SRR1639681 sample and do the alignment and expression quantification, and we will treat it as paired-end data (and make sure that the trimming was done with only the parameters that the authors report)

	$ nohup rsem-calculate-expression -p 10 --paired-end --bowtie2 ~/BCB546_final_project/data/trimmed_files/SRR1639681_1_trimmed.fastq ~/BCB546_final_project/data/trimmed_files/SRR1639681_2_trimmed.fastq ~/BCB546_final_project/data/rsem/asparagus_ref SRR1639681_counting &

If we inspect the SRR1639681_counting.genes.results file, we can see:

	$ less SRR1639681_counting.genes.results
	gene_id transcript_id(s)        length  effective_length        expected_count  TPM     FPKM
	Aspof_comp100003_c0     Aspof_comp100003_c0_seq1        728.00  536.60  1.00    1.78    1.49
	Aspof_comp100012_c0     Aspof_comp100012_c0_seq1        1840.00 1648.59 0.00    0.00    0.00
	Aspof_comp100014_c0     Aspof_comp100014_c0_seq1        619.00  427.62  0.00    0.00    0.00
	Aspof_comp100017_c0     Aspof_comp100017_c0_seq1        569.00  377.65  0.00    0.00    0.00
	Aspof_comp100019_c0     Aspof_comp100019_c0_seq1        569.00  377.65  0.00    0.00    0.00
	Aspof_comp100023_c0     Aspof_comp100023_c0_seq1        374.00  183.83  0.00    0.00    0.00
	
Therefore, if we do all this process per each sample and then we compile all the results per gene and per specimen we will end up with the file that was made as imput for the edgeR analysis of differential expression.

If we want to get the statistics of the alignment process we can call 

    $ rsem-plot-model SRR1639681_counting SRR1639681_models.pdf

This analysis showed that with this sample we got the same number of reads, that is mainly because the initial quality of the reads was very good so the trimming procedure really did not make a big difference to the file. 86% of the reads were mapped to the transcriptome which corresponds with the data that the authors report.