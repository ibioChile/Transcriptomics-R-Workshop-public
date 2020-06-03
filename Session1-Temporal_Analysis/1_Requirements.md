# Session 1 | Requirements

Before starting with these intructions, please watch this video: [Session1-RNAseq_Introduction](https://drive.google.com/open?id=1ogZf5UNYp2iL3x4-b3RRnU5Vmk7qJVBQ).

For this session you will need a **Unix compatible system**, ie. MAC or Linux operating system to work with **Unix Shell** (terminal/command line). Windows10 users can activate the Windows Subsystem for Linux feature to enable instalation of Ubuntu inside Windows10.

Please read the [Workshop Setup](https://github.com/ibioChile/Transcriptomics-R-Workshop/blob/master/Workshop%20Setup.md) document for more info.

More info about how to use **Unix Shell** (command line) [can be found here](https://github.com/ibioChile/Transcriptomics-R-Workshop/blob/master/Command%20line.md).

### **Watch [this video](https://drive.google.com/file/d/1GZbBg3DW284LVvTSucif0fhKskjvMr4X/view?usp=sharing) for a more detailed explanation of the following steps**

## 1. Programs download and instalation

In a bash terminal, create a working folder where you will store documents and files used during this session.

    mkdir test1
    cd test1

Within test1, create a folder to store downloaded programs.

    cd test1
    mkdir Downloads
    cd Downloads

### 1.1 Downloading and working with Trimmomatic (sequencing trimming by quality and adapters).

    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
    unzip Trimmomatic-0.39.zip
    
> An alternative program to download files in command line is "curl", available on Linux and MacOS.  
> example:  
> curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip  

> If wget, unzip, tree or any other system application is not installed on your system:
>> In MacOS (with Homebrew installed) use "brew install \<package\>"  
>> In Ubuntu use "sudo apt install \<package\>"  

### 1.2 Anaconda/Miniconda installation

[Conda](https://conda.io/en/latest) is an open source package management system that quickly installs, runs and updates programs and their dependencies. Once installed, the user creates an 'environment' where any program of interest can be installed. It is possible to create multiple environments, which allow us to install programs whose dependencies are uncompatible. 

>If you already have installed Anaconda or Minicoda, you can proceed with the creation of a working **environment**, which we will use for this session (Section 1.3).

Download links of Anaconda can be found [here](https://www.anaconda.com/distribution/), where the user must select *macOS* or *Linux* operative system accordingly. 

> **Don't download the Windows version!**. If you are a MS_Windows user, remember that our recommendation is to use the Windows Subsystem for Linux (WSL). In that case, Anaconda/Miniconda must be installed inside the Linux system running under Windows (**bash terminal**).

A lightweight alternative to Anaconda is to use Miniconda. This package will install only essential linux resources. Download links can be found [here](https://docs.conda.io/en/latest/miniconda.html)

* **We will use Python 3.7**

Once in the site, click on "Download" and select operative system. Right click on download link and "Copy link".

Now, in your terminal, inside Download folder, run the following command to download file:

    wget <paste the Anaconda link>
    or
    wget <paste the Miniconda link>

Example for Linux:

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

Example for MacOS

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

To install the program, run the following command and follow instructions on screen: 

    bash Miniconda3-latest-Linux-x86_64.sh

During the first step of this installation, the program will show you the 'End User License Agreement'. You must press 'Enter' to scroll down the text, until you have the option of accepting license terms. Type 'yes' to start installation. Press 'Enter' again to confirm the program location. 

Once the installation finishes, you must answer 'yes' to the last message, which will let you inside conda's 'environment1'.  

Now you **must exit the system and enter again**. This means closing the terminal and open it again. Preceeding the default command line, you will see the word **(base)**. Now you are inside the "base" system of anaconda.

The following command prevents Anaconda's "base" system to automatically start, which is recommended. 

    conda config --set auto_activate_base false

Finally, **exit and enter the system again** to update changes on the environment.

### 1.3 Installing programs using Anaconda

Now, let's create an Anaconda environment (conda) and call it 'env1'. In the same command, we will install the programs that will be used during this session on our "env1" environment.
- Install **fastqc**, which will allow us to check reads quality.
- Different **mapping** y **quantifying** programs (in this case, we picked 'Hisat2' for mapping).
#
    conda create -n env1 -c conda-forge -c bioconda fastqc hisat2 samtools

**In addition**, if you wish, you can also install and compare other "mapping" tools (this will take more time and use more space in your computer):

    conda install -n env1 -c conda-forge -c bioconda bowtie2 minimap last star salmon kalisto

>To obtain more information about differences on mapping tools, we recommed to check this forum's discussion:
>[https://bioinformatics.stackexchange.com/questions/4507/better-aligner-than-bowtie2](https://bioinformatics.stackexchange.com/questions/4507/better-aligner-than-bowtie2)

Finally, we **activate** Anaconda's environment:

    conda activate env1

Once we finish our work on this environment, we can exit it by either closing the terminal or using the following command:

    conda deactivate

## 2. Download reads

For this session, we will use data published in the publication ["Temporal transcriptional logic of dynamic regulatory
networks underlying nitrogen signaling and use
in plants"](https://www.pnas.org/content/pnas/115/25/6494.full.pdf). For now, we will only use RNA sequencing data from samples of Arabidopsis thaliana's shoot taken over time: 0, 5, 10, 15, 20, 30, 45, 60, 90 and 120 min. This data was sequenced on Illumina HiSEq 2500 v4 platform (100 bp, PE).

The original files (122 Gb of compressed data, 390 Gb of incompressed data) and a reduced version (4.3 Gb of compressed data) can be found here:

Run these commands to download a *reduced* version of these files (needed for this session):

    cd test1
    mkdir reduced && cd reduced
    wget http://genius.bio.puc.cl/genius/workshop/reduced_reads_session1.tar
    wget http://genius.bio.puc.cl/genius/workshop/reduced_reads_session1.tar.md5
    
    Linux:  
    md5sum reduced_reads_session1.tar  
    
    OSX:  
    md5 reduced_reads_session1.tar  
    
The previous command will produce a hash that must be identical to the hash inside the .md5 file.  
    
Use the following  command to see the .md5 content and compare the hash with the output of the previous command.  

    cat reduced_reads_session1.tar.md5
    
If the hash is identical, then proceed. Otherwise, download the reduced_reads_session1.tar again.  

    tar xf reduced_reads_session1.tar
    cd ..

Run these commands **if you want** to download the *original* version of these files (alternative):  

    cd test1
    mkdir raw && cd raw
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR544/FOLDER/*

This table present the meaning of each sample in terms of sample time and replicate number.

| FOLDER               | Replicate | Time |
|----------------------|-----------|------|
| ```006/SRR5440786``` | 1         | 0    |
| ```005/SRR5440785``` | 2         | 0    |
| ```004/SRR5440784``` | 3         | 0    |
| ```004/SRR5440814``` | 1         | 5    |
| ```002/SRR5440832``` | 2         | 5    |
| ```003/SRR5440823``` | 3         | 5    |
| ```005/SRR5440815``` | 1         | 10   |
| ```003/SRR5440833``` | 2         | 10   |
| ```004/SRR5440824``` | 3         | 10   |
| ```004/SRR5440834``` | 1         | 15   |
| ```005/SRR5440825``` | 2         | 15   |
| ```006/SRR5440816``` | 3         | 15   |
| ```005/SRR5440835``` | 1         | 20   |
| ```006/SRR5440826``` | 2         | 20   |
| ```007/SRR5440817``` | 3         | 20   |
| ```006/SRR5440836``` | 1         | 30   |
| ```007/SRR5440827``` | 2         | 30   |
| ```008/SRR5440818``` | 3         | 30   |
| ```007/SRR5440837``` | 1         | 45   |
| ```008/SRR5440828``` | 2         | 45   |
| ```009/SRR5440819``` | 3         | 45   |
| ```008/SRR5440838``` | 1         | 60   |
| ```009/SRR5440829``` | 2         | 60   |
| ```000/SRR5440820``` | 3         | 60   |
| ```009/SRR5440839``` | 1         | 90   |
| ```000/SRR5440830``` | 2         | 90   |
| ```001/SRR5440821``` | 3         | 90   |
| ```000/SRR5440840``` | 1         | 120  |
| ```001/SRR5440831``` | 2         | 120  |
| ```002/SRR5440822``` | 3         | 120  |

    cd ..

Furthermore, we need to download from TAIR, *Arabidopsis thaliana*'s genome and other files:

    mkdir databases && cd databases
    wget ftp://ftp.ensemblgenomes.org/pub/release-34/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
    wget ftp://ftp.ensemblgenomes.org/pub/release-34/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.34.gtf.gz
    gunzip *.gz

With the previous steps we are ready to follow the [Data pre-processing instructions](https://github.com/ibioChile/Transcriptomics-R-Workshop/blob/master/Session1-Temporal_Analysis/2_Pipeline-Data_preprocessing.md).

<END>
