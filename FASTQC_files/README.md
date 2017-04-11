
#FASTQC_files for inspecting raw RNA sequencing data

*Author: Valeria Velasquez Zapata*

This folder contains a serie of FASTQC analyses that were made in order to show an example of the quality assesing and trimming that the authirs of the paper have done before obtaining the counting data. 
### Getting data

The data is contained in the NCBI database, so, in order to get it we can download an configure toolkit from ncbi

    $ wget "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz"

    $ tar -xzf sratoolkit.current-centos_linux64.tar.gz

After downloading and decompressing it we can get the raw sequences by calling fastq-dump on the SRA experiment number, as we can see in this web [page](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=259909):

Asparagus officinalis line 2 XX female postmeiotic bud(SRR1642915)


    $ ./sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump SRR1642915

Asparagus officinalis line 8A XX female spear tip(SRR1639282)

    $ /sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump SRR1639282

### Getting and running FASTQC

To download FASTQC we use:

    $ curl -O http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.10.1.zip
    $ unzip fastqc_v0.10.1.zip
    $ cd FastQC
    $ chmod +x fastqc
    $ cd ..

To run [FASTQC](http://wiki.bits.vib.be/index.php/Linux_command_line#Automating_FASTQC_analyses)

    $ ./FastQC/fastqc SRR1642915.fastq
    $ ./FastQC/fastqc SRR1639282.fastq

After looking at the results it was possible to see that both dataset need to be trimmed. In order to achieve that we are using CUTADAPT. For the SRR1642915 datset we will perform a quality, kmer and adapter trim, since those were the specific issues that we found after the FASTQC analysis. Then, for the SRR1639282 dataset we will perform a quality and kmer trimming.

### Installing CUTADAPT

To install CUTADAPT we can use:

    $ pip install --user --upgrade cutadapt

For help:

    $ cutadapt -h

We will use a qphred threshold of 20 to trim the data:

    $ cutadapt -q 20 -o SRR1639282_trimmed.fastq SRR1639282.fastq

Then we can trim out the kmers, which are localized at the beginning of the sequences (about 10b when we see the FASTQC results)

    $ cutadapt -u 10 -o SRR1639282_Quali_10bp_trimmed.fastq SRR1639282_trimmed.fastq

Each of those files can be analyzed by FASTQC angain, in order to see the read quality improvement:

    $ ./FastQC/fastqc SRR1639282_trimmed.fastq
    $ ./FastQC/fastqc SRR1639282_Quali_10bp_trimmed.fastq

After looking at the results we could see that the reads were sucessfully cleaned and ready for mapping to the transcriptome data.

For SRR1642915 there are overrepresented sequences that correspond to adapters, so it is necessary to trim them out

    $ cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATG -o SRR1642915_AdapTrimmed.fastq SRR1642915.fastq
    $ cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGC -o SRR1642915_AdapTrimmed2.fastq SRR1642915_AdapTrimmed.fastq 
  
Then we can trim them by quality and remove the kmers
  
    $ cutadapt -q 20 -o SRR1642915_Adapt_Quali_trimmed.fastq SRR1642915_AdapTrimmed2.fastq 
    $ cutadapt -u 9 -o SRR1642915_Adapt_Quali_9bp_trimmed.fastq SRR1642915_Adapt_Quali_9bp_trimmed.fastq

Finally, we re run the FASTQC analysis and we can track the changes of each trimming procedure:

    $ ./fastqc SRR1642915_AdapTrimmed.fastq
    $ ./fastqc SRR1642915_AdapTrimmed2.fastq
    $ ./fastqc SRR1642915_Adapt_Quali_trimmed.fastq
    $ ./fastqc SRR1642915_Adapt_Quali_9bp_trimmed.fastq

### Using MUlTIQC to summarize the results

For installing MULTIQC we do:
    
    $ pip install multiqc

Then for runnig it:

    $ multiqc .

Then, it is earier to compare and look at the results. We can also look at each result separately
