# Homework 2: RNA-Seq data analysis

*Written by Pamela Camejo*  
*Updated in July 2023 - Jonathan Maldonado*

This homework tests some of the knowledge acquired during the three "RNA-seq data analysis" sessions of this workshop and it is a **requirement** to receive your certificate of completion. 
* The assignment is due Sunday, July 30 (17hrs).
* **Code with figures must be prepared using RMarkdown format (.Rmd) in RStudio** except for Part 1. You can download and use any of the published .Rmd files as template (sessions 2 to 4).
* After completion, **send Us an email with all requested files and your RMarkdown (Part 2-5)**.
* Replies to the questions of the homework must be included inside the RMarkdown file, following the corresponding command.

You should not use any artificial intelligence tool to assist in the development of the homework, otherwise you will lose the possibility to obtain the course certificate.

Just in case: [info about RMarkdown](https://rmarkdown.rstudio.com/lesson-1.html) and a [useful Cheatsheet](https://www.rstudio.com/wp-content/uploads/2015/02/rmarkdown-cheatsheet.pdf).

In this assignment, we will work with RNA-sequencing data obtained from two strains of *Saccharomyces cerevisiae* (*S. cerevisiae*) when exposed to a micotoxin called patulin. You can find more information of this data in the article [Distinct Transcriptional Changes in Response to Patulin Underlie Toxin Biosorption Differences in Saccharomyces cerevisiae](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6669508/). This data was sequenced on Illumina HiSeq 4000 (150 bp, PE). 

The goal of this homework is to analyze this data, starting from reads pre-processing, mapping to the genome of *S. cerevisiae* and identification of DEGs, expression patterns and their functions. 

The original sequencing files can be downloaded [here](https://www.ebi.ac.uk/ena/browser/view/PRJNA540854), but for Part 1, you will work with a reduced version of these files (1M reads).

## Part 1: Reads pre-processing.

The first part of this pipeline can be performed using command line. Follow the steps in [Session 1, Part 1 | Requirements](https://github.com/ibioChile/Transcriptomics-R-Workshop-public/blob/master/Session1-RNAseq_Data_Preprocessing/1_Requirements.md) and in [Session 1, Part 2 : RNA-Seq data preprocessing](https://github.com/ibioChile/Transcriptomics-R-Workshop-public/blob/master/Session1-RNAseq_Data_Preprocessing/2_Pipeline-Data_preprocessing.md) to perform this analysis.

Prepare a simple text file to record any used command and functions including your comments. Name it ```commands_SCR64.txt```.

1.1	Use ```wget```, ```curl -O``` or any other way to download the reduced version of the original files from this link and extract files:

	http://genius.bio.puc.cl/workshop/reduced_homework2.tar
	http://genius.bio.puc.cl/workshop/reduced_homework2.tar.md5
	
> Extract files using ```tar xvf file.tar```

1.2	Download the genome and .gtf files of *S. cerevisiae* from the following link. Unpack using ```gunzip file.gz```.

	ftp://ftp.ensembl.org/pub/release-84/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
	ftp://ftp.ensembl.org/pub/release-84/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.84.gtf.gz

1.3	Check reads quality using *fastqc*. Store results in the folder ```raw_fastqc```.

1.4	Use Trimommatic to trim adapters and regions of poor quality. Set minimum lenght of 100bp in Trimmomatic. Store trimmed files in folder ```trimmed```.

1.5	Check the quality of the trimmed reads using *fastqc*. Store results in the folder ```trimmed_fastqc```.

1.6	Use HiSat2 to map trimmed reads to the genome of *S. cereviseae*. Don't forget to index the genome first.

1.7	Use Rsubread to generate a table with counts per gene.

1.8	Save this table, name it ```reduced_counts_SCR64.tsv```.

1.9	Send Us. both files: ```commands_SCR64.txt``` and ```reduced_counts_SCR64.tsv```

## Part 2: Data analysis 

Checkout pipelines in [Session 1, Part 2 : RNA-Seq data preprocessing](https://github.com/ibioChile/Transcriptomics-R-Workshop-public/blob/master/Session1-RNAseq_Data_Preprocessing/2_Pipeline-Data_preprocessing.md) and [Session 3, Part 2: Gene expression analysis, multivariate experiments](https://github.com/ibioChile/Transcriptomics-R-Workshop-public/blob/master/Session3-Treatment_and_Multivariate/3_Pipeline-Differential-expression-multivariate.md) to perform this analysis. 

Download ```fc0.counts.txt``` and ```metadata_Sc_Patulin.txt``` from the Homework2 folder (github).   
TIP: Open each github file, look at the "Raw" button, right clik and select "save as".

2.1	Import the downloaded counts and metadata table. Sort samples in ```metadata``` table, first by Strain and then by Treatment. Format ```counts``` headers, so they have the same order than samples in ```metadata``` table.
> How many genes and samples are found in this table?

2.2	Create an edgeR object called ```dgList``` with this data. Use the parameters ```counts``` and ```genes``` to create this object. Add metadata to the object.

2.3	Normalize data using the ```calcNormFactors``` function. Normalize using two different methods: "TMM" and "upperquartile". Store normalized objects on ```dgList_TMM``` and ```dgList_UPQ```.

2.4	Create boxplots of:
- Raw data.
- Normalization with "UPQ"
- Normalization with "TMM"

> Make sure you use log2(cpm), replacing zero counts by ones (use ```cpm``` function with ```log2=TRUE``` and ```prior.count = 1```) when plotting data. 
> Use ```par(mfrow=c(1,3))``` to plot the three graphs next to each other.
> Remember to add titles to each graph.

> What's the difference between the output of these methods?

Carry out the rest of the pipeline with data normalized using TMM.
	
2.5	Use function ```filterByExpr``` to filter out genes with  < 10 counts. **No need to group data**.

2.6	Create a density plot comparing raw and normalized/filtered data. Make sure you use cpm with log2 and 0=1. **No need to color lines**.

2.7	Create a PCA plot with normalized/filtered data. Color samples by Treatment and use a different shape for each strain.

> Which of the two factors, Treatment or Strain, has a greater influence on sample clustering? Print your answer.

## Part 3: DGEs analysis

Checkout pipeline in [Session 3, Part 2: Gene expression analysis, multivariate experiments](https://github.com/ibioChile/Transcriptomics-R-Workshop-public/blob/master/Session3-Treatment_and_Multivariate/3_Pipeline-Differential-expression-multivariate.md) to perform this analysis. 

3.1	First, create a ```model.matrix``` using metadata information. In this matrix, we are going to evaluate both individual factors (Strain and Treatment) as the interaction of both variables. In order to accomplish that, the formula used for the design should be: 
	
	~ Strain + Treatment + Strain:Treatment
	
> Also, you could have used ```~ Strain*Treatment``` for the model. Look for more information about formulas for model matrices [here](https://genomicsclass.github.io/book/pages/expressing_design_formula.html).
	
3.2	Estimate common, trended and gene specific dispersion.

3.3	Fit a negative binomial generalized log-linear model with the function ```glmFit```.

> Check coefficients of this model. Some of them will be related to the individual variables, and some to the interaction of both factors.
  
3.4	Create 3 lists of DEGs. Adjust p-values with 'fdr' method and use as threshold p-value = 0.01:

- ```strain_DEG_sig```: List of genes differentially expressed in each strain.
- ```treat_DEG_sig```: List of genes differentially expressed by treatment.
- ```int_DEG_sig```: List of genes differentially expressed by strain and treatment.

> How many DEGs are found in each list? Print your answer.

3.5	Create a table with cpm counts (log2 and 0=1) of each of these gene sets. Name this tables:```strain_DEG_sig_exp```,```treat_DEG_sig_exp``` and ```int_DEG_sig_exp```.

3.6	Create a heatmap of gene expression for each of these gene sets with the function ```Heatmap``` from the *ComplexHeatmap* package. Include in the heatmap a top bar colored by Treatment and Tissue.

> Don't forget to scale the data before creating the heatmap.

## Part 4: DGEs clusters

Checkout pipeline in [Session 2, Part 1 : Gene expression analysis](https://github.com/ibioChile/Transcriptomics-R-Workshop-public/blob/master/Session2-Temporal_Analysis/1_Pipeline-Gene_expression_analysis.md) to perform this analysis. For this section we are only going to work with the table ```int_DEG_sig_exp```. 

4.1	We will use hierarchical clustering to partition genes into different clusters. There’s two steps to this clustering procedure:
- Calculate a “distance” metric between each pair of genes.
- Cluster the genes hierarchically.

In R, we can use the dist() function for the first step and hclust() for the second step. Create a gene dendogram using these functions as shown here:

	dend <- hclust(dist(scaledMatrix))
	
> Make sure you use a scaled version of int_DEG_sig_exp. The same that you used before to create the Heatmap.

4.2	Plot dendogram using ```plot(dend)```. Add a line cutting the tree at h=5.

4.3	Create a set of clusters named ```geneClust``` with the function ```cutree``` at height = 5.

> How many clusters were generated? Print your answer.

4.4	Obtain the mean expression (core) per sample of each cluster. Calculate cores for each cluster.

> Use a scaled version of int_DEG_sig_exp to calculate cores.

4.5	Create a molten table of cores and Sample Names. Follow similar steps to the ones shown in Session 2, 7.4.2.
> Instead of Time use Sample_Name.
> 
> The variable "Time" is a number but the variable Sample_Name is a character so, the dmolten data table must be adjusted with this command to be sure that expression values are used as "numbers".   
> ```dmolten$value <- as.numeric(dmolten$value)```
> 
> Before plotting, sort your samples order with this code: ```dmolten$Sample_Name <- factor(dmolten$Sample_Name ,levels = dgList_TMM$metadata$Sample_Name)```
> 

4.6	Make a ```geom_tile``` plot with the mean expression of each core (cluster vs. Sample Name). Fill tiles by expression. This will create a heatmap where every row represents the mean expression of each cluster across samples. 
> If ```geom_tile``` don't work in your R session, try loading the library ```farver```.
> Add x and y-axis titles. 
> Use ```scale_fill_continuous(type = "viridis") ``` to change color.

## Part 5: DGEs function

Checkout pipeline in [Session 3, Part 1 : Gene expression analysis, control/treatment experiment](https://github.com/ibioChile/Transcriptomics-R-Workshop-public/blob/master/Session3-Treatment_and_Multivariate/2_Pipeline-Differential-expression-control-treatment.md) to perform this analysis. For this section we are only going to work with the table ```int_DEG_sig_exp```. 

5.1	Download gene ontology terms of S. cerevisiae from this link http://current.geneontology.org/annotations/sgd.gaf.gz and unzip it.

5.2	Import it to your R session.

> Use ```read.delim``` as shown in Session 2.1 adding the command ```quote = ""```

5.3	This table has a different format than the one used during Session 2. We will need to extract gene IDs from column 11 of this file. Use ```sapply``` and ```strsplit``` with this column to return the first string (Gene ID) (HINT: For strsplit, use "|" as split and fixed=TRUE). Store this vector containing gene IDs as ```gene_id```.

5.4 Create a data frame called ```goframeData```. The first column of this data frame should be ```gene_id``` and the other ones, columns 3, 5 and 7 from ".gaf" file.

5.5	Build a customizable database for ViSEAGO.

5.6	Create a topGOdata object considering the following information: 
> Use as “background universe” the total list of genes being expressed. You can obtain this information from ```dgList_TMM$genes```.
> Use “biological processes (BP)” category.
> Use as ```selection``` the list of genes in ```int_DEG_sig_exp``` list.
> Use ```nodeSize = 5```

5.7	Use these parameters to perform an enrichment test over DEGs. Run "topGO" enrichment test with parameters ```algorithm ="classic",  statistic = "fisher"```. Display results as a table and save it as "GO_intersect.BP.txt". Read written table and store information as ```BP_sResults_table```.
> What is the most significantly enriched (lowest p-value) gene ontology (GO) term for these genes. Print your answer. Include GO ID and term.

5.8	Plot gene ontology (GO) using ```showSigOfNodes```. Save plot using ```printGraph```. ```showSigOfNodes``` will need the ```Rgraphviz``` package installed.
