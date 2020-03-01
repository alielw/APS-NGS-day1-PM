**_Advanced Data Analysis - Introduction to NGS data analysis_**<br>
*Department of Animal and Plant Sciences, University of Sheffield*
APS Advanced Stats Delivery https://github.com/visoca/AdvDataAna-introNGS

# Aligning Illumina RNA-seq data
#### Alison Wright, Nicola Nadeau, Helen Hipperson

The aim of this practical is to learn how to align Illumina RNA-seq data to a reference genome and assemble transcripts. We will be using a dataset of expression data of multiple individuals of *Heliconius melpomene*.

## Table of contents

1. Prepare reference genome

2. Align RNA-seq reads to reference

3. Visualise alignments

4. Assess mapping quality

5. Assemble transcripts

---

## Initial set up
First of all, this tutorial must be run using an interactive session in ShARC. You will also submit jobs to ShARC. For that, you should log in into ShARC with `ssh`, and then request an interactive session with `qrsh`. Your shell prompt should show `sharc-nodeXXX` (XXX being a number between 001 and 172) and not `@sharc-login1` nor `@sharc-login2`.

For this particular tutorial, we are going to create and work on a directory called `align` in your /fastdata/$USER directory. Remember you need to replace $USER with your username.

        cd /fastdata/$USER
        mkdir align
        cd align

Remember you can check your current working directory anytime with the command `pwd`.
It should show something like:

        pwd
        /fastdata/$USER/align

The data we will be using today is the folder below. You need to be in interactive mode `qrsh` to access this data.

	/usr/local/extras/Genomics/workshops/NGS_AdvSta_2020/NGS_data

---

## Programmes

There are a number of different programs that can be used to align RNA-seq reads to a reference genome. There are two types of aligners: Splice-unaware and Splice-aware. Splice-unaware aligners are not aware of exon/intron junctions and are therefore only appropriate for mapping to a transcriptome. Splice-aware aligners map reads over exon/intro junctions and are appropriate for aligning reads to a genome reference. This is an important point to consider when choosing an aligner and why it is necessary to use different programs for RNA- versus DNA-seq data.

Potential aligners include [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml), [STAR](https://github.com/alexdobin/STAR) and [Tophat2](https://ccb.jhu.edu/software/tophat/index.shtml). We will be using HISAT2. This is the next generation of spliced aligner from the same research group that developed TopHat and has considerable advantages in terms of speed and memory requirements. When choosing an aligner appropriate to your study, it is useful to refer to the numerous studies that compare performance statistics across programs (eg [PeerJ review](https://peerj.com/preprints/27283/), [PLoS One review](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190152)). The blog [RNA-Seq Blog](https://www.rna-seqblog.com/review-of-rna-seq-data-analysis-tools/) is a valuable resource.  

---

## 1. Prepare reference genome

Prepare the reference genome before mapping RNA-seq reads with [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml)

* **HISAT2 index**

        >hisat2-build -f reference_genome.fa reference_genome_db

This command indexes the reference genome fasta file and outputs six files: .1.ht2  .2.ht2  .3.ht2  .4.ht2  .5.ht2  .6.ht2 ... These files constitute the index.

## a. PRACTICAL ACTIVITY

Make an index of the reference genome.

* First load an interactive session on ShARC.

        qrsh

* Make a new folder in your align folder in fastdata. Remember to change $USER to your username.

        mkdir /fastdata/$USER/align/ref
        cd /fastdata/$USER/align/ref

* Copy the reference genome to this new folder

        cp /usr/local/extras/Genomics/workshops/NGS_AdvSta_2020/NGS_data/Reference/Hmel2.fa /fastdata/$USER/align/ref
        
* Index the genome

        hisat2-build -f Hmel2.fa Hmel2

This commands takes around five minutes to finish. While this is running, move onto section two to familiarise yourself with HISAT2 and prepare your command to map reads to the reference genome. When `hisat2-build` has finished you can then run HISAT2.
   
## 2. Align RNA-seq reads to reference

Map RNA-seq reads to a reference genome using [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml).
 
 * **Map reads to genome with HISAT2**
 
        >hisat2 reference_genome forward.fastq reverse.fastq
		
There are many different mapping parameters you can specify, see [here](https://ccb.jhu.edu/software/hisat2/manual.shtml). While it is often suffient to run HISAT2 with default settings, there are a number of parameters that should be considered:

**-q** Reads are FASTQ files. FASTQ files usually have extension .fq or .fastq. FASTQ is the default format.

**--dta/--downstream-transcriptome-assembly** Report alignments tailored for transcript assemblers including StringTie. With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computationa and memory usage. You must specify this option when running StringTie after HISAT2 to quantify gene expression.

**-S path** File to write SAM alignments to. By default, alignments are written to the "standard out" or "stdout" filehandle (i.e. the console).

**-p**  Use this many threads to align reads. The default is 1.

**--phred33** Input qualities are ASCII chars equal to the Phred quality plus 33. This is also called the "Phred+33" encoding, which is used by the very latest Illumina pipelines.

**--met-file path** Write hisat2 metrics to file <path>. Having alignment metric can be useful for debugging certain problems, especially performance issues. See also: --met. Default: metrics disabled.

The following two parameters should be specified if it is necessary to obtain very accurate and specific mapping data. A pair of reads that align with the expected relative mate orientation and range of distances between mates is said to align “concordantly”. If both mates have unique alignments, but the alignments do not match paired-end expectations (i.e. the mates aren’t in the expected relative orientation, or aren’t within the expected distance range, or both), the pair is said to align “discordantly”. Discordant alignments may be of particular interest, for instance, when seeking structural variants.

**--no-discordant** For paired reads, report only concordant mappings.

**--no-mixed** For paired reads, only report read alignments if both reads in a pair can be mapped. By default, when hisat2 cannot find a concordant or discordant alignment for a pair, it then tries to find alignments for the individual mates. This option disables that behavior.

HISAT2 outputs a SAM file with the alignment information for each read.

## b. PRACTICAL ACTIVITY

Map RNA-seq reads to the reference genome. This practical will illustrate the effect of different parameters on the stringency of mapping. In particular, we will focus on the effect of `**no-mixed**`. In some cases, such as for SNP calling, we may want to ensure that mapping is particularly stringent and that we have high confidence a read is correctly aligned. We can ensure higher confidence in our mapping by specifying `no-mixed`, which means that a read is only reported if both pairs align to the reference. You can compare mapping statistics between two Tophat runs, one with `no-mixed` and one without.

You should submit these commands as jobs to ShARC. Our data is paired end, strand-specific and has quality scores in phred33 format. We will focus on one individual (60A). 

* First, make some output folders.

        mkdir /fastdata/$USER/align/Trimmed_data
        
        mkdir /fastdata/$USER/align/HISAT2
        mkdir /fastdata/$USER/align/HISAT2/60A
 
* Copy two fastq files into your output folder. These are subsampled versions of the files you generated yourself to speed up computational time in this practical. You must be in interactive mode to do this.

        cp /usr/local/extras/Genomics/workshops/NGS_AdvSta_2020/NGS_data/Trimmed_files/60A_1.fq.gz /fastdata/$USER/align/Trimmed_data
        cp /usr/local/extras/Genomics/workshops/NGS_AdvSta_2020/NGS_data/Trimmed_files/60A_2.fq.gz /fastdata/$USER/align/Trimmed_data

* Make an executable script where you can specify the job requirements. 
        
        cd /fastdata/$USER/align/HISAT2/60A
        nano HISAT2_60A.sh
        
* Specify the requirements for the hpc and give the path to the Genomics Software Repository.

        #!/bin/bash
        #$ -l h_rt=00:15:00
        #$ -l rmem=5G
        #$ -pe smp 3
        #$ -wd /fastdata/$USER/align/HISAT2/60A
        
        source /usr/local/extras/Genomics/.bashrc

* Next, add the following commands to map reads to the reference genome.

        hisat2 \
		/fastdata/$USER/align/ref/Hmel2 \
		-1 /fastdata/$USER/align/Trimmed_data/60A_1.fq.gz \
		-2 /fastdata/$USER/align/Trimmed_data/60A_2.fq.gz \
		-p 3 \
		-q \
		--dta \
		-S /fastdata/$USER/align/HISAT2/60A/60A.sam \
		--met-file /fastdata/$USER/align/HISAT2/60A/60A.stats
        
* Finally, submit your job. Only submit this if `hisat2-build` has finished and you have an indexed reference genome. This should take around 10 minutes to run so move onto the next step.
        
        qsub HISAT2_60A.sh
        
* Next, repeat this process to run HISAT2 again but now specifying --no-mixed. You need to create a new working directory and executable script. It is essential that you specify the new output folder and working directory to ensure files aren't overwritten. You need to replace `$USER` with your username.

        mkdir /fastdata/$USER/align/HISAT2/60A_nomixed
        
        cd /fastdata/$USER/align/HISAT2/60A_nomixed
        
        nano HISAT2_60A_nomixed.sh
        
        #!/bin/bash
        #$ -l h_rt=00:15:00
        #$ -l rmem=5G
        #$ -pe smp 3
        #$ -wd /fastdata/$USER/align/HISAT2/60A_nomixed
        
        source /usr/local/extras/Genomics/.bashrc
        
        hisat2 \
		/fastdata/$USER/align/ref/Hmel2 \
		-1 /fastdata/$USER/align/Trimmed_data/60A_1.fq.gz \
		-2 /fastdata/$USER/align/Trimmed_data/60A_2.fq.gz \
		-p 3 \
		-q \
		--dta \
		--no-mixed \
		-S /fastdata/$USER/align/HISAT2/60A/60A.nomixed.sam \
		--met-file /fastdata/$USER/align/HISAT2/60A/60A.nomixed.stats

* Run the script. This will take ~10 minutes to run so move onto the next step.
        
       qsub HISAT2_60A_nomixed.sh
        
---

## 3. Visualise alignments

* **Interactive Genome Viewer**

We can view read alignments to the reference genome with [IGV](http://software.broadinstitute.org/software/igv/) on the Desktop. 

You need to first sort the SAM file with [Samtools](http://www.htslib.org/doc/samtools-1.0.html). By default Samtools sorts files by coordinate. It is easier to handle BAM files as they are compressed so at the same time we convert to a BAM format.

        samtools view -b sample.sam | samtools sort > sample.sorted.bam

You need to then index the BAM file with [Samtools](http://www.htslib.org/doc/samtools-1.0.html).

	samtools index sample.sorted.bam

Then you have to transfer the SAM file, index and reference genome to the Desktop.

Load the reference genome into IGV

        Genomes/Load Genome from File

Load the SAM file into IGV

        File/Load from File

## c. PRACTICAL ACTIVITY

We can view read alignments to the reference genome with [IGV](http://software.broadinstitute.org/software/igv/) on the Desktop.
	
* You need to sort the SAM file you generated (60A.sam) with [Samtools](http://www.htslib.org/doc/samtools-1.0.html)

		samtools view -b 60A.sam | samtools sort > 60A.sorted.bam

* You need to index the new BAM file.

		samtools index 60A.sorted.bam

* Lets extract reads aligning to one scaffold (Hmel200115). We can do this with [Samtools](http://www.htslib.org/doc/samtools-1.0.html)
	
	        samtools view -b -h 60A.sorted.bam "Hmel200115" > 60A_Hmel200115.bam

* You need to index the new BAM file with [Samtools](http://www.htslib.org/doc/samtools-1.0.html)
	
	        samtools index 60A_Hmel200115.bam
	
* copy the BAM file, index and reference genome to your desktop.
	
	        /fastdata/$USER/align/Quality_assessment/96I_Hmel200115.bam
	
	        /fastdata/$USER/align/Quality_assessment/96I_Hmel200115.bam.bai
	
	        /fastdata/$USER/align/ref/Hmel2.fa
                
* Download [IGV](http://software.broadinstitute.org/software/igv/download) on Windows to your Desktop.

* Unzip downloaded file. 

* Click on the Windows icon on the bottom left of the screen. Type `cmd` into the search code. Click on `x86 Nature Tools Command Prompt for VS`. A terminal will open.

* Drag and drop the unzipped file called igv in the IGV folder onto the terminal window. Press enter. There may be an error message asking for permissions. Just click cancel if this happens. IGV should now load.
	
* Load the reference genome into IGV
	
	        Genomes/Load Genome from File
	
* Load the BAM file into IGV

	        File/Load from File

* Lets look at the scaffold `Hmel200115`. According to the gff file, there is one annotated coding gene on this scaffold with three exons.

		grep "Hmel20011" /usr/local/extras/Genomics/workshops/NGS_AdvSta_2020/NGS_data/Reference/Hmel2.gff
        
        	scaffold        feature start   stop
        	Hmel200115	gene	1	1683	
        	Hmel200115	exon	1	415
        	Hmel200115	exon	496	1280
        	Hmel200115	exon	1564	1683

* Using the coordinates above, locate the scaffold & exons in the IGV viewer. 

		Search Hmel200115 in the search box in IGV viewer. You can also use the scaffold dropdown.

* Do reads map to these exons? Can you see the intronic regions? Is there any evidence for additional exons that have not been annotated?
        
---

## 4. Assess mapping quality

Whilst IGV allows us to examine specific regions of the genome, it not easy to summarise this information across the genome. It is useful to measure the general performance of the aligner for a number of reasons including i) choosing the best performing aligner, ii) optimising mapping parameters, or iii) identifying a subset of high quality reads. There are a number of tools we can use to calculate quality statistics of mapping quality.

* **Tophat Output**

Tophat produces an output log file (align_summary.txt) with various statistics including i) reads processed, ii) number of reads mapped iii) number of pairs mapped. Eg:

        Left reads:
                Input     :  18214710
                Mapped   :   7665434 (42.1% of input)
                of these:    366828 (4.8%) have multiple alignments (4783 have >20)
        Right reads:
                Input     :  18214710
                Mapped   :   7678151 (42.2% of input)
                of these:    368728 (4.8%) have multiple alignments (4830 have >20)
        42.1% overall read mapping rate.

        Aligned pairs:   6851101
        of these:        337896 (4.9%) have multiple alignments
                         27231 (0.4%) are discordant alignments
        37.5% concordant pair alignment rate.

* **BAM format**

Information about read mapping is specified in the BAM file. BAM and SAM formats are designed to contain the same information. The [SAM format](https://en.wikipedia.org/wiki/SAM_(file_format)) is more human readable, and easier to process by conventional text based processing programs. The BAM format provides binary versions of most of the same data, and is designed to compress reasonably well.

For example, you might want to identify reads with high quality alignments. This may be because you want to extract a set of reads for downstream analyses. This information is encoded in the BAM/SAM file as [flags](https://samtools.github.io/hts-specs/SAMv1.pdf). The following table gives an overview of the mandatory fields in the SAM format:

![alt text](https://github.com/alielw/APS-NGS-day1-PM/blob/master/SAM%20fields.jpg)

The fields that will be most useful are the FLAGS and MAPQ. The Broad Institute have a useful website to understand the [FLAGS](https://broadinstitute.github.io/picard/explain-flags.html). In addition, the MAPQ field indicates the mapping quality of the read.

You can visualise a BAM file using [Samtools](http://www.htslib.org/doc/samtools-1.0.html).

        >samtools view file.bam

You can also get general statistics about the BAM files eg the number of mapped reads using [Samtools](http://www.htslib.org/doc/samtools-1.0.html).

        >samtools flagstat file.bam
        
Eg:
        
        16989295 + 0 in total (QC-passed reads + QC-failed reads)
        1645710 + 0 secondary
        0 + 0 supplementary
        0 + 0 duplicates
        16989295 + 0 mapped (100.00% : N/A)
        15343585 + 0 paired in sequencing
        7665434 + 0 read1
        7678151 + 0 read2
        13310548 + 0 properly paired (86.75% : N/A)
        13702202 + 0 with itself and mate mapped
        1641383 + 0 singletons (10.70% : N/A)
        70752 + 0 with mate mapped to a different chr
        39370 + 0 with mate mapped to a different chr (mapQ>=5)

## d. PRACTICAL ACTIVITY

Count how many reads have all of the following features; are paired, first in the pair, on the reverse strand and also have a read mapped in proper pair. This information is encoded in the FLAG section of the SAM/BAM file. We can do this using a bam file we already generated using Tophat (96I.bam). We can do this in interactive mode.

* Lets move into the folder you created earlier containing the BAM file. 
     
        cd /fastdata/$USER/align/Quality_assessment

* Then, view the first few lines of the BAM file with samtools (`samtools view BAM`), Unix pipe (`|`) and Unix `head`. Compare this output with the table above. Identify which column contains the FLAG field.

        samtools view 96I.bam | head

* Next, print the FLAG field for the first few lines using samtools (`samtools view BAM`), Unix pipe (`|`) and cut (`cut -f[Column number]`), Unix pipe (`|`) and Unix `head`. 

* Then, using the [Broad Institute website](https://broadinstitute.github.io/picard/explain-flags.html) identify the FLAG which specifies that a read is all of the following. The FLAG will be an integer. 
        
        paired
        is mapped in proper pair
        first in pair
        on the reverse strand

* Finally, count how many reads have this FLAG using samtools (`samtools view BAM`), Unix pipe (`|`) and cut (`cut -f[Column number]`), Unix pipe (`|`) and grep (`grep -c "FLAG"`).

* Ask a demonstrator for the correct answer!

---

## 5. Assemble transcripts

We can use [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/) to assemble gene transcripts. Cufflinks is part of the  “classic” RNA-Seq workflow, which includes read mapping with TopHat followed by assembly with Cufflinks. As Tophat is being superceded by [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml), Cufflinks is being replaced by [StringTie](https://ccb.jhu.edu/software/stringtie/). Although often you may have a set of annotated genes in the reference genome, and therefore a reference gtf, this may be incomplete and some genes may not be annotated. Cufflinks will identify potential new transcripts.

* **Generate gtf file for each sample**
A GTF file contains information about assembled transcripts. We can use cufflinks to generate a set of transcripts for each sample. 

        >cufflinks sample.bam 
    
As with Tophat, there are many different mapping parameters you can specify, see [here](http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html). While it is often sufficient to run Cufflinks with default settings, there are a number of parameters that should be considered:

**–GTF-guide** Tells Cufflinks to use the supplied reference annotation a GFF/GTF file to guide RABT assembly. Reference transcripts will be tiled with faux-reads to provide additional information in assembly. Output will include all reference transcripts as well as any novel genes and isoforms that are assembled.

**–library-type** Library type

Cufflinks produces four output files: transcripts.gtf, isoforms.fpkm_tracking, genes.fpkm_tracking, skipped.gtf

**transcripts.gtf**
This GTF file contains Cufflinks’ assembled isoforms. The first 7 columns are standard GTF, and the last column contains attributes, some of which are also standardized (“gene_id”, and “transcript_id”). There one GTF record per row, and each record represents either a transcript or an exon within a transcript. Further description of the format can be found [here](https://github.com/alielw/APS-NGS-day1-PM/blob/master/GTF%20format.jpg)

* **Generate list of gtf files**

Generates a text file with a list (one per line) of GTF files to merge together into a single GTF file.

        >ls */transcripts.gtf > gtf_list.txt

* **Merge gtf files**

Cufflinks includes a script called [cuffmerge](http://cole-trapnell-lab.github.io/cufflinks/cuffmerge/index.html) that you can use to merge together several Cufflinks assemblies. We must merge transcripts to get a complete reference set of transcripts across all individuals.

        >cuffmerge -g reference.gtf - s gtf_list.txt

## e. PRACTICAL ACTIVITY

Cufflinks takes a couple of hours to run, so we have already generated gtf files for each sample. In this practical, we will merge all the gtf files to generate a complete set of transcripts. We will also include the reference gff as Cufflinks will use it to attach gene names and other metadata to the merged catalog.

* First, copy the folder containing the gtf files to your fastdata folder

        cp -r /usr/local/extras/Genomics/workshops/NGS_AdvSta_2020/NGS_data/Cufflinks_output /fastdata/$USER/align/

* Copy the reference gff to the new folder. This contains all the annotated transcripts in the genome.

        cp /usr/local/extras/Genomics/workshops/NGS_AdvSta_2020/NGS_data/Reference/Hmel2.gff /fastdata/$USER/align/Cufflinks_output

* View list of gtf files

        cd /fastdata/$USER/align/Cufflinks_output

        for i in */*transcripts.gtf; do echo $i; done

* Generate list of all gtf files

         for i in */*transcripts.gtf; do echo $i; done > gtf_list.txt

* Merge gtf files with cuffmerge. Cuffmerge produces a GTF file, merged.gtf that merges together the input assemblies.

        cuffmerge -g Hmel2.gff -s /fastdata/$USER/align/ref/Hmel2.fa gtf_list.txt 

* Count how many transcripts there are. This should be equal to the number of lines minus the four summary lines at the end of the file. You can count line numbers with Unix (`wc -l`).
        
* Correct a bug in cuffmerge. We noitced that there are a few reads that aren't stranded which then generate an error. These need to be removed from the gtf file with the following:

        grep -v -P "\t\.\t\.\t\." merged_asm/merged.gtf > merged_asm/merged_allstranded.gtf

---
Return to main course page: https://github.com/visoca/AdvDataAna-introNGS
