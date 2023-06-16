# Homework 1: R Basics and Visualization

*Written by Pamela Camejo*  
*Updated in May 2022 - Nate Johnson*

This homework tests some of the knowledge acquired during 'R Basics' and 'R Visualization' online classes and it is a **requirement** to continue with the following 'RNA-seq analysis' sessions. We will have one session on Friday May 20 2022 to discuss and help with any challenges. The assignment is due **Monday May 23 2022 at 5:00pm** and must be uploaded by that time. The assignment can be submitted via "dropbox requests" using the following link:  [link to upload folder](https://www.dropbox.com/request/ksKczHquA6ZIJgCOfADQ)

| Field | Info |
| ------- | ---- |
| **Help session** | Friday June 23 2022 at 11:00am |
| **Due date** | Saturday June 24 2022 at 10:00am |
| **Upload link** | [link to upload folder](https://www.dropbox.com/request/dgFJl6uRYQDt9v9qx6ox) |
| **Upload format** | Save as an R script from R studio in the following name **Homework1_[Nombre].[Apellido].R** |

The goal of this homework is to start working with data resulting from RNA-sequencing. We will import data, process it, create functions, load libraries, use some popular functions (some of them different from the ones used in online classes) and create and modify plots. We will work with data from the 'pasilla' package. This package provides per-gene read counts computed for selected genes from RNA-seq data that were presented in the article ["Conservation of an RNA regulatory map between Drosophila and mammals"](https://www.ncbi.nlm.nih.gov/pubmed/20921232). 

* Start an R Studio session, create a new R Script file and save it as **Homework1_[FirstName].[LastName].R**. Write an script to follow the next instructions:

1. **Download the dataset**, you need to install the pasilla library from Bioconductor then load this library. **HINT**: Use the function of BiocManager to install this package. You will find the file directory using the following command. 

        >counts = system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
        >counts

Import counts table as 'pasillaCountsTable'. Make sure that column names are equal to the sample names and row names to gene IDs.

> This can also be found in the github folder fo the assignment where you can directly import this table directly without installing the 'pasilla' package. See read.delim().

2. Write a line in R to **print the answer** to each of this questions.

 2.1. How many rows (genes) in the table?
 2.2. How many columns (samples) in the table?

3. **Create metadata table**. Create a data frame called 'pasillaMetadataTable' containing the following sample information:

<center> 
        
| | condition   | libType |
|:-|:--------------:|-------------:|
| **untreated1** |  untreated  | single-end  |
| **untreated2** |  untreated  | single-end  |
| **untreated3** |  untreated  | paired-end  |
| **untreated4** |  untreated  | paired-end  |
| **treated1** |  treated  | single-end  |
| **treated2** |  treated  | paired-end  |
| **treated3** |  treated  | paired-end  |

 </center>

**HINT**: Use the function *data.frame* to generate this table. Row.names of this table (first column in the table shown here) correspond to the sample names, which are 'pasillaCountTable' column names.

4. What is the class of each column in metadata table? Write a code line to **print this answer**. **HINT**: Check *lapply* function to perform this task in one line.

5. **Select samples** in metadata with 'paired-end' sequencing. Name this new table 'paired'.

6. **Create a new table** from 'pasillaCountsTable' with 'paired-end' sequenced samples. Name this new table 'countTable_paired'.

7. **Create a function 'median_of_ratios'** to normalize counts using the ['median of ratios' method](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html). In this method, each count is divided by a sample-specific size factor. The following diagram describes the steps followed to perform this normalization:

<pre><code>	┌─────────────────────────┐
	│    Step 1: Calculate    │
	│row-wise geometric mean. │
	└───────────┬─────────────┘
	            │              
	┌───────────▼─────────────┐
	│ Step 2: Calculate ratio │
	│of each count to the row │
	│     geometric mean.     │
	└───────────┬─────────────┘
	            │              
	┌───────────▼─────────────┐
	│    Step 3: Calculate    │
	│column-wise median value.│
	│This is the normalization│
	│ factor for each sample  │
	│     (size factor).      │
	└───────────┬─────────────┘
	            │              
	┌───────────▼─────────────┐
	│    Step 4: Calculate    │
	│  normalized values by   │
	│ dividing each count by  │
	│ the sample size factor. │
	└─────────────────────────┘</code></pre>

Following this diagram, the 'median_of_ratios' function that you will design here should take as input any table with counts (formatted like 'pasillaCountsTable') and perform the following calculations:

7.1 Create a vector equal to the geometric mean across all samples. This vector cotains one value per each gene and should have the same number of rows that input table. Use this formula to calculate the geometric mean:                                             

        exp(mean(log(x)))   # Where x is a vector of numbers
 
  **HINT**: Use *rowMeans* to calculate the average of rows in a matrix or data.frame.
        
7.2 Divide each row of counts by the geometric mean. This new table is known as 'table of ratios'.
 
7.3 Calculate the median of each column of 'table of ratios' (ignore NA in median calculation). **HINT**: Check function *colMedians*, which is part of the package *matrixStats*. The input of this function should be a matrix and not a data.frame.
  
7.4 Finally, calculate the normalized count values using the normalization factor (median of columns). This is performed by dividing each count value in a given sample by that sample’s normalization factor. **HINT**: Check the function *apply* to divide each column by a different number, then transpose the resulting matrix.
  
7.5 Print the normalized table inside the function.

8. **Use this function** to normalize 'countTable_paired' and store results in 'norm_count'.

9. **Filter** rows with at least 5 counts in at least half of the samples. **HINT**: Use the function *rowSums* to sum logical TRUEs in table.

10. How many genes were removed? Write a code to **print this answer**.

11. **Log normalize filtered table**. Add 1 count to each value before applying log2. Store this normalized and filtered log table as 'log_filt_norm_table'.

12. Before plotting, we need to **convert the filtered table into a molten data frame**. For that, we will use the function *melt*, which is part of the package *reshape2*. After installing this package and loading the library, use *melt* with 'log_filt_norm_table'. Store melted table as 'melted_log_table'. Change colnames of 'melted_log_table' to: Gene, Sample and Expression.

13. Add an extra column with the corresponding treatment variable to 'melted_log_table'. **HINT**: Use *match* with the column "Sample" from 'melted_log_table' and row.names of 'pasillaMetadataTable' to create an INDEX with the matching data. Use this index to select data from 'pasillaMetadataTable$condition'. Then, you can use the function *cbind* to bind this column to 'melted_log_table'. When using 'cbind' name the new column "Treatment". This is an example of how to perform this step:

        >to_plot <- cbind(melted_log_table, "Treatment" = pasillaMetadataTable$condition[INDEX])
 
We will use this table to create some plots.

14. Install *ggplot2* package and load library.

15. **Create density plot** to visualize the distribution of expression values (x axis). Additionally:
* Color by sample and facet by treatment.
* Add title "Density plot of expression values". The title should be centered in the plot.
* Change x-axis title to 'log2(expression)' and remove y-axis label.
* Change the size of x and y axis titles to 10.
* Save plot as 'density_plot.pdf', adjust width and height to improve visualization.

16. **Create boxplots** of expression for each sample. Additionally:
* Color by treatment.
* Change colors to 'red' for treated and 'blue' for untreated.
* Add title "Boxplot of expression values". The title should be centered in the plot.
* Change x-axis title to 'Sample' and y-axis title to 'log2(expression)'.
* Use theme_bw().
* Save plot as 'boxplot.pdf', adjust width and height to improve visualization.
