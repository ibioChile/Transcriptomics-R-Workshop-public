# Session 1 | RNA-Seq data preprocessing

Before starting with these intructions, please watch this video: [Session1-RNAseq_Introduction](https://drive.google.com/file/d/1ogZf5UNYp2iL3x4-b3RRnU5Vmk7qJVBQ/view). 
Also, don't forget to read and follow the [Requirements](https://github.com/ibioChile/Transcriptomics-R-Workshop-public/blob/master/Session1-RNAseq_Data_Preprocessing/1_Requirements.md) section.

### **Watch [this video](https://drive.google.com/file/d/104scpjSD8ZeEIcqTHo0JNqueY5RGxJdH/view) for a more detailed explanation of the following pipeline**

We must activate the conda **env1** environment that we prepared with the Requirements instructions. This way, all required programs and packages will be accessible for the following pipeline.

    conda activate env1

## 1. Quality Control.

### 1.1 Quality control check.

We should always check reads quality, even when the authors/sequencing company state that data was already checked. For this step, we will use **fastqc**. 
> The modifier **-t** indicates the number of processors to be used by the program... modify it accordingly to your machine specs.

    cd test1/
    mkdir raw_fastqc
    fastqc -t 8 reduced/*.gz -o raw_fastqc
    
Within **raw_fastqc** folder, we will find results as different **.html** files, which can be opened with any web browser.

If you are using the WSL system (Linux in Windows), in order to be able to use your favorite Windows web explorer to see this results you must copy the **.html** files to any folder that could be seen by Windows (for example, the root of the disk C).

First we create a folder to store the data:

    mkdir /mnt/c/RNAseq-tutorial
    
Now we copy the data to the created folder:

    cp -r raw_fastqc /mnt/c/RNAseq-tutorial
    
Then, you can use your **Windows file explorer** and reach the **.html** files in C:\RNAseq-tutorial\raw_fastqc and open it by double-clic.

These links can give you more information about how to interpret the output of **fastqc**:  
[fastqc tutorial and faq with good/bad images examples](https://rtsf.natsci.msu.edu/genomics/tech-notes/fastqc-tutorial-and-faq/)  
[fastqc documentation with good/bad reports examples](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

### 1.2 Trimming adapters and regions of poor quality.

For this step, we will use the program Trimmomatic. First, we copy a file containing adapters sequence.

    cp ~/test1/Downloads/Trimmomatic-0.39/adapters/TruSeq3-PE.fa .

At this point, we have two ways to proceed with the adapter trimming process:

#### Option 1) To run one Trimmomatic command for each pair of sequence files

The "standard" command is as follows using the first sequence pairs as example:

    java -jar ~/test1/Downloads/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 SRR5440784_1_1M.fastq.gz SRR5440784_2_1M.fastq.gz -baseout trimmed/SRR5440784_1M_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:10:30 MINLEN:50 AVGQUAL:25;

So, if you want, you can copy and paste 30 times that command, **replacing the names of the files with the names of the fastq.gz files** and **replacing the baseout statement by the corresponding one**.

#### Option 2) To write a bash script to process all files

We create a file containing a script to automatically run Trimmomatic on each read file pair. You can use **nano** editor to create this file:

    nano doTrimming.sh

Then, paste the following command, **exit and save**:
> The modifier **-threads** indicates the number of processors to be used by the program... don't forget to review this parameter before saving the file.

    mkdir -p trimmed/
    for f in reduced/*1_1M.fastq.gz; do
		base=${f##*/};
		java -jar ~/test1/Downloads/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 $f ${f%_*_*}_2_1M.fastq.gz -baseout trimmed/${base%_*_*}_1M_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:10:30 MINLEN:50 AVGQUAL:25;
    done

To execute the file:

    bash doTrimming.sh

This script will create a folder called **trimmed**, where compressed trimmed (paired and unpaired) files are created. During the next steps, we will only use paired trimmed files.

### 1.3 Check results of quality trimming.
Once the process is finished we can check the results using fastqc again.

We repeat the command used in **1.1** to check the quality of trimmed reads:
> The modifier **-t** indicates the number of processors to be used by the program... don't forget to review this parameter before running the command.

    mkdir trimmed_fastqc
    fastqc -t 8 trimmed/*P.fastq.gz -o trimmed_fastqc

Remember that, **as mentioned in section 1.1**, if you are using the WSL system (Linux in Windows), in order to be able to use your favorite Windows web explorer to see this results you must copy the **.html** files to any folder that could be seen by Windows (for example, the root of the disk C).

## 2. Mapping and reads quantification

Now that we have prepared our reads, we can align the reads to an existing reference genome of Arabidopsis. The current version of the reference genome is TAIR10. Here we will use [HiSat2](http://daehwankimlab.github.io/hisat2/) to align these reads. HiSat2 is the descendent of TopHat, one of the first widely-used aligners, but alternative mappers could be used, such as STAR, salmon, kalisto, etc. There are often numerous mapping parameters that we can specify, but usually the default mapping parameters are fine. However, library type (paired-end vs single-end) and library strandness (stranded vs unstranded) require some different settings when mapping and counting, so they are two important pieces of information to know about samples. This Arabidopsis data comprises unstranded, paired-end reads so we will specify that where necessary. HiSat2 can output a mapping summary file that tells what proportion of reads mapped to the reference genome. As we’re only using a subset of 1,000,000 reads per sample, aligning should just take a few minutes or so. To run the full samples from this dataset would take longer.

### 2.1. HiSat2: index

The HiSat2 program requires an index file for every genome in order to do the mapping. Some pre-computed indexes can be downloaded from the program's website, but for plants and many other species, we have to prepare the index files using HiSat2. We will build the index files for the Arabidopsis genome, using the function **hisat2-build**. This will generate a set of 6 files with suffixes .1.ht2, .2.ht2, .3.ht2, .4.ht2, .5.ht2, .6.ht2, .7.ht2.  Athaliana is the prefix assigned to output.
> The modifier **-p** indicates the number of processors to be used by the program... don't forget to review this parameter before running the program.

    cd databases/
    gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
    hisat2-build -p 8 Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Athaliana
    cd ..

### 2.2 HiSat2: mapping

At this point, we have two ways to proceed with the sequence mapping process:

#### Option 1) To run a set of commands for each pair of sequence files

The "standard" command is as follows using the first trimmed sequence pairs as example:

    hisat2 -p 8 --dta -x databases/Athaliana -1 SRR5440784_1M_trimmed_1P.fastq.gz -2 SRR5440784_1M_trimmed_2P.fastq.gz -S mapping/SRR5440784_1M_trimmed_Athaliana.sam 2> mapping/SRR5440784_1M_trimmed.align.stat.txt
    samtools sort -@ 8 -o mapping/SRR5440784_1M_trimmed_Athaliana.bam mapping/SRR5440784_1M_trimmed_Athaliana.sam
    rm mapping/SRR5440784_1M_trimmed_Athaliana.sam

> The modifier **-p** indicates the number of processors to be used by the program... don't forget to review this parameter before saving the file.

So, if you want, you can copy and paste 30 times these three commands, **replacing the names of the files with the names of the fastq.gz files** and **replacing the baseout statement by the corresponding one**.

#### Option 2) To write a bash script to process all files

We will write an script with the command below to align each set of reads from shoot tissues to the genome of Arabidopsis. HiSat2 generates a SAM file with mapped reads, Sequence Alignment/Map (SAM) format is a generic alignment text file that describes the genome coordinates were the reads were aligned, the number of possible match sites in the genome and other related information. Before we move to next step, we need to convert output SAM files to sorted BAM files. Therefore, we add an extra line to our script where we use Samtools to convert SAM to BAM files.

    nano doHisat2.sh
   
Paste the following text, **save and exit**:
> The modifier **-p** indicates the number of processors to be used by the program... don't forget to review this parameter before saving the file.

    mkdir -p mapping
    for fn in trimmed/*_1M_trimmed_1P.fastq.gz;
    do	
    	base=${fn##*/};
    	echo "Processing sample ${base%_*}"
    	hisat2 -p 8 --dta -x databases/Athaliana -1 $fn -2 ${fn%_*}_2P.fastq.gz -S mapping/${base%_*}_Athaliana.sam 2> mapping/${base%_*}.align.stat.txt
    	samtools sort -@ 8 -o mapping/${base%_*}_Athaliana.bam mapping/${base%_*}_Athaliana.sam
    	rm mapping/${base%_*}_Athaliana.sam
    done

Run the script:

    bash doHisat2.sh
    
The ".align.stat.txt" files contains the stats of the alignment.

### 2.3 Rsubread: quantification

To assign the reads to a gene, we use the R software (R Core Team 2015). In this guide we will use the Rsubread package that facilitates the RNA-seq read data analyses, proportionating several metrics: quality assessment of sequence reads, read alignment, read summarization, exon-exon junction among others (Liao, Smyth, & Shi, 2013). 
**This package is available in Unix/Linux/Mac/Windows, but in Windows sometimes it's not easy to install.**

> If you are using the WSL system (Linux in Windows), must read the note at the end of this tutorial.
    
First, open Rstudio. Then, make sure that the working directory is "WorkDir". To install Rsubread package, run the following instructions:

	>setwd("~/test1/mapping") # or the path where you store the bam files.

	>BiocManager::install("Rsubread")
	>library(Rsubread)
	
 Then, it is necessary to run the feature counts function of this library for all SAM files. First, select the SAM files:

	>bam.list <- dir(pattern=".bam")
	
Since there are different types of library preparation and run configurations, the user should use the appropriate one. In this example, we use the non-stranded (strandSpecific=0) option. The following command runs the program using 6 threads (modifier **-nthreads**) and store the output in "fc0" object. 

	>fc0 <- featureCounts(bam.list, annot.ext= "../databases/Arabidopsis_thaliana.TAIR10.34.gtf",isGTFAnnotationFile=T, allowMultiOverlap=T, isPairedEnd=T, nthreads=6, strandSpecific=0)
	
Then, we save the counts and associated stats in tab delimited files:

	>write.table(fc0$counts, "../fc0.counts.txt", sep="\t", col.names=NA, quote=F)
	>write.table(fc0$stat, "../fc0.stat.txt", sep="\t", row.names=F, quote=F)

Now, with **fc0.counts.txt** file you are ready to start the [RNAseq analysis pipeline in R](https://github.com/ibioChile/Transcriptomics-R-Workshop-public/tree/master/Session2-Temporal_Analysis).

### 2.4 Rsubread: quantification in Windows using Linux with WSL system

If you are using the WSL system (Linux in Windows), you can run R in the bash terminal to complete the Rsubread steps and then, with the results, you can move to Windows Rstudio.

The steps are:

Install R on a conda environment

    conda deactivate
    conda create -n env2 -c conda-forge r-base=4.1.3
    conda activate env2
 
Enter to the R environmment

    R

Install required packages

    >install.packages("BiocManager")
    >BiocManager::install("Rsubread")
    
Follow "2.3 Rsubread: quantification" steps

    >library(Rsubread)
    >setwd("~/test1/mapping") # or the path where you store the bam files.
    
    >bam.list <- dir(pattern=".bam")
    >head(bam.list) # a list of filenames
    >str(bam.list) # 30 files as expected
    
    # In the following command, the modifier **-nthreads** indicates the number of processors to be used by the program 
    >fc0 <- featureCounts(bam.list, annot.ext= "../databases/Arabidopsis_thaliana.TAIR10.34.gtf",isGTFAnnotationFile=T, allowMultiOverlap=T, isPairedEnd=T, nthreads=8, strandSpecific=0)
    
    >write.table(fc0$counts, "../fc0.counts.txt", sep="\t", col.names=NA, quote=F)
    >write.table(fc0$stat, "../fc0.stat.txt", sep="\t", row.names=F, quote=F)

    # Once all done, you can quit R command line interface with this command
    >q()
    ## Respond "n" to the following question as this time it is not necessary: Save workspace image? [y/n/c]
    
Now, you can copy the **fc0.counts.txt** to your Windows filesystem and use it with Windows RStudio to follow the [RNAseq analysis pipeline in R](https://github.com/ibioChile/Transcriptomics-R-Workshop-public/tree/master/Session2-Temporal_Analysis).

First we create a folder to store the data (in case you don't created it on the previous steps):

    mkdir /mnt/c/RNAseq-tutorial
    
Now we copy the data:

    cp -r fc0.counts.txt /mnt/c/RNAseq-tutorial/

Then, you can use your **Windows file explorer** and reach the **fc0.counts.txt** in C:\RNAseq-tutorial\


