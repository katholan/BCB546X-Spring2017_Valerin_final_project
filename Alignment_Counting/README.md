# Alignment and expression estimation with RSEM

*Author: Valeria Velasquez Zapata*

This folder contains the results of the attempt of read aligment and quantification of expression, which corresponds to the step after the read cleaning procedure, which is located in the FASTQC_files folder. For this analysis we will take the SRR1642915 sample. 

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

There are several modes to run Tophat, we sill try 2 of them: very sensitive and very fast: 

	$ nohup tophat2 -o very_sensitive --b2-very-sensitive asparagus_transcriptome SRR1642915_Adapt_Quali_9bp_noShortReads_trimmed.fastq &
 
	$ nohup tophat2 -o very_fast --b2-very-sensitive asparagus_transcriptome SRR1642915_Adapt_Quali_9bp_noShortReads_trimmed.fastq &

### RSEM analysis for counting

For more detailed information about RSEM you can visit [https://github.com/bli25ucb/RSEM_tutorial](https://github.com/bli25ucb/RSEM_tutorial)

To install RSEM 

    $ wget https://github.com/deweylab/RSEM.git
    $ make install
 
#### Preparing reference

The first command that we will use will take the annotation file (supplementary information table S3 [http://onlinelibrary.wiley.com/doi/10.1111/nph.13389/full](http://onlinelibrary.wiley.com/doi/10.1111/nph.13389/full)

We formatted the file, by taking just the first two columns and adding Annotation file. Having the annotation file from the supplementary information and the  Aspof_ prefix to each entry. Once this was dome we are ready for building the reference for the counting analysis. We should build this reference based on the trancriptome file, since there is no genome available for this specie yet (by using the `--transcript-to-gene-map` function)

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
   
Therefore, if we do all this process per each sample and then we compile all the results per gene and per specimen we will end up with the file that was made as imput for the edgeR analysis of differential expression.
