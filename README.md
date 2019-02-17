**_Advanced Data Analysis - Introduction to NGS data analysis_**<br>
*Department of Animal and Plant Sciences, University of Sheffield*

# Aligning Illumina RNA-seq data
#### Victor Soria-Carrasco

The aim of this practical is to learn how to align Illumina RNA-seq data to a reference genome. We will be using a dataset of expression data of 8 individuals of *Heliconius melpomene*. For each individual, two different wing regions have been sequenced.

## Table of contents

1. Prepare reference genome

2. Align RNA-seq reads to reference

3. Visualise alignments

4. Assess mapping quality

5. Assemble transcripts

---
---

## Initial set up
First of all, this tutorial must be run using an interactive session in ShARC. You will also submit jobs to ShARC. For that, you should log in into ShARC with `ssh`, and then request an interactive session with `qrsh`. Your shell prompt should show `sharc-nodeXXX` (XXX being a number between 001 and 172) and not `@sharc-login1` nor `@sharc-login2`.

For this particular tutorial, we are going to create and work on a directory called `align` in your /fastdata/$USER directory:
```bash
cd /fastdata/$USER
mkdir align
cd align
```
Remember you can check your current working directory anytime with the command `pwd`.
It should show something like:
```bash
pwd
```
>`´/fastdata/myuser/align`´<br>

---
---

## Programmes

There are a number of different programs that can be used to align RNA-seq reads to a reference genome. There are two types of aligners: Splice-unaware and Splice-aware. Splice-unaware aligners are not aware of exon/intron junctions and are therefore only appropriate for mapping to a transcriptome. Splice-aware aligners map reads over exon/intro junctions and are appropriate for aligning reads to a genome reference and the analysis of gene expression levels. This is an important point to consider when choosing an aligner and why it is necessary to use different programs for RNA- versus DNA-seq data.

Potential aligners include [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml), [STAR](https://github.com/alexdobin/STAR) and [Tophat2](https://ccb.jhu.edu/software/tophat/index.shtml). We will be using Tophat2 in this practical as it is a standard, widely used approach for Illumina data. However, HISAT2 is the next generation of spliced aligner from the same research group that developed TopHat and has considerable advantages in terms of speed and memory requirements. When choosing an aligner appropriate to your study, it is useful to refer to the numerous studies that compare performance statistics across programs (eg [PeerJ review](https://peerj.com/preprints/27283/), [PLoS One review](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190152). The blog [RNA-Seq Blog](https://www.rna-seqblog.com/review-of-rna-seq-data-analysis-tools/) is a valuable resource.  

## 1. Prepare reference genome

Prepare the reference genome before mapping RNA-seq reads with [Bowtie2](http://bowtie-bio.sourceforge.net/tutorial.shtml)

* **Tophat2 index**

        >bowtie2-build reference_genome.fa reference_genome_db

This command indexes the reference genome fasta file and outputs four files: reference_genome_db.1.ebwt, reference_genome_db.2.ebwt, reference_genome_db.rev.1.ebwt, and reference_genome_db.rev.2.ebwt. These files constitute the index.

## a. PRACTICAL ACTIVITY

Make an index of the reference genome.

* First load an interactive session on ShARC.

        qrsh

* Make a new folder in your align folder in fastdata

        mkdir /fastdata/$USER/align/ref

* Copy the reference genome to this folder

        cp /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/NGS_data/Heliconius_melpomene.Hmel1.dna.toplevel.fa.gz /fastdata/$USER/align/ref

* We need to load the bashrc file containing the path to bowtie2.
        
        source /usr/local/extras/Genomics/.bashrc
        
* Index the genome

        bowtie2-build Heliconius_melpomene.Hmel1.dna.toplevel.fa.gz Hmel_db
   
## 2. Align RNA-seq reads to reference

Map RNA-seq reads to a reference genome using [Tophat2](https://ccb.jhu.edu/software/tophat/index.shtml).
 
 * **Map reads to genome with Tophat2**
 
        >tophat2 reference_genome_db sample1_1.fq sample1_2.fq
    
There are many different mapping parameters you can specify, see [here](https://ccb.jhu.edu/software/tophat/manual.shtml). While it is often suffient to run Tophat2 with default settings, there are a number of parameters that should be considered:
    
**--mate-inner-dist** This is the expected (mean) inner distance between mate pairs. For, example, for paired end runs with fragments selected at 300bp, where each end is 50bp, you should set -r to be 200. The default is 50bp. You can obtain this information from the sequencing centre.
    
**--solexa-quals or --solexa1.3-quals** You need to specify the quality score system of your data. --solex-quals tells Tophat to use the Solexa scale for quality values in FASTQ files. --solexa1.3-quals specifies that quality scores are encoded in Phred-scaled base-64. Use this option for FASTQ files from pipeline 1.3 or later.

**--num-threads** Use this many threads to align reads. The default is 1.

**--library-type** The default is unstranded (fr-unstranded). Can specify fr-firststrand or fr-secondstrand.

**--GTF** Supply TopHat with a set of gene model annotations and/or known transcripts, as a GTF 2.2 or GFF3 formatted file. If this option is provided, TopHat will first extract the transcript sequences and use Bowtie to align reads to this virtual transcriptome first. Only the reads that do not fully map to the transcriptome will then be mapped on the genome. The reads that did map on the transcriptome will be converted to genomic mappings (spliced as needed) and merged with the novel mappings and junctions in the final tophat output.

**-p**  Use this many threads to align reads. The default is 1.

The following two parameters should be specified if it is necessary to obtain very accurate and specific mapping data. A pair of reads that align with the expected relative mate orientation and range of distances between mates is said to align “concordantly”. If both mates have unique alignments, but the alignments do not match paired-end expectations (i.e. the mates aren’t in the expected relative orientation, or aren’t within the expected distance range, or both), the pair is said to align “discordantly”. Discordant alignments may be of particular interest, for instance, when seeking structural variants.

**--no-discordant** For paired reads, report only concordant mappings.

**--no-mixed** For paired reads, only report read alignments if both reads in a pair can be mapped (by default, if TopHat cannot find a concordant or discordant alignment for both reads in a pair, it will find and report alignments for each read separately; this option disables that behavior).

Tophat outputs a BAM file with the alignment information for each read.

## b. PRACTICAL ACTIVITY

Map RNA-seq reads to the reference genome. This practical will illustrate the effect of different parameters on the stringency of mapping. In particular, we will focus on the effect of **no-mixed**. In some cases, such as for SNP calling, we may want to ensure that mapping is particularly stringent and that we have high confidence a read is correctly aligning. We can ensure higher confidence in our mapping by specifying no-mixed, which means that a read is only mapped if both pairs align to the reference. You will compare mapping statistics between two Tophat runs, one with no-mixed and one without. This activity will comprise part of the course assessment.

You should submit these commands as jobs to Sharc. Remember to specify different output folders for each command so that you don't overwrite any folders (-o). Our data is paired end, strand-specific and has quality scores in phred33 format. We will focus on one individual (60A). 

* First, make an output folder.

        mkdir /fastdata/$USER/align/60A
 
* Copy the two fastq files into your output folder. You must be in interactive mode to do this.

        cp /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/NGS_data/Trimmed_files/60A.trimA_1.fastq.gz /fastdata/$USER/align/60A
        cp /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/NGS_data/Trimmed_files/60A.trimA_2.fastq.gz /fastdata/$USER/align/60A
        
* Make an executable script where you can specify the job requirements. 
        
        emacs Tophat_60A.sh
        
* Tophat2 normally requires around 10Gb of memory and across 10 threads will take around 10 hours to finish. Specify this information in the executable script. You can use the UNIX command pwd to get the full path of a folder. Remember to specify the .bashrc file which includes the path to Tophat2.

        #$ -l h_rt=10:00:00
        #$ -l rmem=10G
        #$ -pe smp 10
        #$ -wd /fastdata/$USER/align/60A
        
        source /usr/local/extras/Genomics/.bashrc
        
* Next, add the following commands to map reads to the reference without specifying no-mixed. You need to replace PATH with paths to the relevant folders.

        tophat2 -p 10\
        /fastdata/$USER/align/ref/Hmel_db\
        --library-type fr-firststrand\
        -o /fastdata/$USER/align/60A\
        /fastdata/$USER/align/60A/60A.trimA_1.fastq.gz /fastdata/$USER/align/60A/60A.trimA_2.fastq.gz 
        
* Finally, submit your job
        
        qsub Tophat_60A.sh
        
* Next, repeat this process to run Tophat again but now specifying --no-mixed. You need to create a new working directory and executable script. It is essential that you specify the new output folder and working directory to ensure files aren't overwritten.

        mkdir /fastdata/$USER/align/60A_nomixed
        emacs Tophat_60A_nomixed.sh
        
        #$ -l h_rt=10:00:00
        #$ -l rmem=10G
        #$ -pe smp 10
        #$ -wd /fastdata/$USER/align/60A_nomixed
        
        source /usr/local/extras/Genomics/.bashrc

        tophat2 -p 10\
        PATH/Hmel_db\
        --library-type fr-firststrand\
        --no-mixed\
        -o /fastdata/$USER/align/60A_nomixed\
        /fastdata/$USER/align/60A_nomixed/60A.trimA_1.fastq.gz /fastdata/$USER/align/60A_nomixed/60A.trimA_2.fastq.gz 

## 3. Visualise alignments

* **Interactive Genome Viewer**

## b. PRACTICAL ACTIVITY

## 4. Assess mapping quality

Whilst IGV allows us to examine specific regions of the genome, it not easy to summarise this information across the genome. It is useful to measure the general performance of the aligner for a number of reasons including i) choosing the best performing aligner, ii) optimising mapping parameters, or iii) identifying a subset of high quality reads. There are a number of tools we can use to calculate quality statistics of mapping quality.

* **Tophat Output**

Tophat produces an output log file with various statistics including i) reads processed, ii) number of reads mapped iii) number of pairs mapped etc

* **BAM format**

Information about read mapping is specified in the BAM file. BAM and SAM formats are designed to contain the same information. The [SAM format](https://en.wikipedia.org/wiki/SAM_(file_format)) is more human readable, and easier to process by conventional text based processing programs. The BAM format provides binary versions of most of the same data, and is designed to compress reasonably well.

For example, you might want to identify reads with high quality alignments. This may be because you want to extract a set of reads for downstream analyses. This information is encoded in the BAM/SAM file as [flags](https://samtools.github.io/hts-specs/SAMv1.pdf). The following table gives an overview of the mandatory fields in the SAM format:

![alt text](https://github.com/alielw/APS-NGS-day2-AM/blob/master/SAM%20fields.jpg)

The fields that will be most useful are the FLAGS and MAPQ. The Broad Institute have a useful website to understand the [FLAGS](https://broadinstitute.github.io/picard/explain-flags.html). In addition, the MAPQ field indicates the mapping quality of the read.

You can visualise a BAM file using [Samtools](http://www.htslib.org/doc/samtools-1.0.html).

        >samtools view file.bam

You can also get general statistics about the BAM files eg the number of mapped reads using [Samtools](http://www.htslib.org/doc/samtools-1.0.html).

        > samtools flagstat file.bam

## c. PRACTICAL ACTIVITY

Count how many reads have all of the following features; are paired, first in the pair, on the reverse strand and also have a read mapped in proper pair. This information is encoded in the FLAG section of the SAM/BAM file. We can do this using a bam file we already generated using Tophat (96I.bam). We can do this in interactive mode.

* Make a new folder and transfer the bam files into this folder. 

        mkdir /fastdata/$USER/align/Quality_assessment
        cp /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/NGS_data/Tophat_output/96I.bam /fastdata/$USER/align/Quality_assessment

* Then, view the first few lines of the bam file with samtools (samtools view BAM), Unix pipe (|) and Unix head. Compare this output with the table above. Identify which column contains the FLAG field.

        samtools view 96I.bam | head

* Next, using samtools (samtools view BAM), with Unix pipe (|) and cut (cut -f["Column number"]), extract the FLAG field.

* Then, using the [Broad Institute website](https://broadinstitute.github.io/picard/explain-flags.html) identify the FLAG which specifies that a read is all of the following:
        
        paired
        is mapped in proper pair
        first in pair
        on the reverse strand

* Finally, using samtools (samtools view BAM), Unix pipe (|) and cut (cut -f["Column number"]), Unix pipe (|) and grep (grep -c "FLAG"), count how many reads have this FLAG.

* Ask a demonstrator for the correct answer!

## 5. Assemble transcripts

We can use [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/) to assemble gene transcripts. Cufflinks is part of the  “classic” RNA-Seq workflow, which includes read mapping with TopHat followed by assembly with Cufflinks. As Tophat is being superceded by [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml), Cufflinks is being replaced by [StringTie](https://ccb.jhu.edu/software/stringtie/). Although often you may have a set of annotated genes in the reference genome, and therefore a reference gtf, this may be incomplete and some genes may not be annotated. Cufflinks will identify potential new transcripts.

* **Generate gtf file for each sample**
A GTF file contains information about assembled transcripts. We can use cufflinks to generate a set of transcripts for each sample. 

        >cufflinks sample.bam 
    
As with Tophat, there are many different mapping parameters you can specify, see [here](http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html). While it is often suffient to run Cufflinks with default settings, there are a number of parameters that should be considered:

**–GTF-guide** Tells Cufflinks to use the supplied reference annotation a GFF file to guide RABT assembly. Reference transcripts will be tiled with faux-reads to provide additional information in assembly. Output will include all reference transcripts as well as any novel genes and isoforms that are assembled.

**–library-type** Library type

Cufflinks produces three output files: transcripts.gtf, isoforms.fpkm_tracking, genes.fpkm_tracking

**transcripts.gtf**
This GTF file contains Cufflinks’ assembled isoforms. The first 7 columns are standard GTF, and the last column contains attributes, some of which are also standardized (“gene_id”, and “transcript_id”). There one GTF record per row, and each record represents either a transcript or an exon within a transcript. Further description of the format can be found [here](https://github.com/alielw/APS-NGS-day2-AM/blob/master/GTF%20format.jpg)

* **Generate list of gtf files**

Generates a text file with a list (one per line) of GTF files to merge together into a single GTF file.

        >ls */transcripts.gtf > gtf_list.txt

* **Merge gtf files**

Cufflinks includes a script called [cuffmerge](http://cole-trapnell-lab.github.io/cufflinks/cuffmerge/index.html) that you can use to merge together several Cufflinks assemblies along with the reference set of already annotated transcripts. We must merge transcripts to get a complete reference set of transcripts across all individuals.

        >cuffmerge -g reference.gtf - s gtf_list.txt

## d. PRACTICAL ACTIVITY

Cufflinks takes a couple of hours to run, so we have already generated gtf files for you for each sample. In this practical, we will merge all the gtf files, together with the reference set of annotated transcripts, to generate a complete set of transcripts.

* First, copy the folder containing the gtf files to your fastdata folder

        cp -r /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/NGS_data/Cufflinks_output /fastdata/$USER/align/

* Copy the reference gtf to the new folder. This contains all the annotated transcripts in the genome.

        cp -r /usr/local/extras/Genomics/workshops/NGS_AdvSta_2019/NGS_data/Heliconius_melpomene.Hmel1.42.gtf /fastdata/$USER/align/Cufflinks_output

* Generate list of gtf files

        ls */transcripts.gtf > gtf_list.txt

* Merge gtf files with cuffmerge

        cuffmerge -g Heliconius_melpomene.Hmel1.42.gtf - s gtf_list.txt

        
