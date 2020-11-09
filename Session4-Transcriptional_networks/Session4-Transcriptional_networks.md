Session 4 - Transcriptional networks
================
tm
5/28/2020

The rapid increase in the availability of transcriptomics data
represents an opportunity for biologists to integrating and interpreting
these data sets. The opportunity is to use this information to generate
testable hypothesis to understand molecular mechanisms controlling gene
expression and biological processes. Asuccessful strategy to generate
tractable hypotheses from transcriptomics data has been to buildnetwork
graphs based on patterns of gene co-expression. This guide includes
basic instructions for the operation of widely used open source
platforms R and Cytoscape. Even though the data we used in this example
was obtained from Arabidopsis thaliana, the workflow developed in this
guide can be easily adapted to work with RNA-seq data from any organism.
The instructions are based in this book chapter:
<https://pubmed.ncbi.nlm.nih.gov/29525965/>

### Importing and formatting data.

First set your working directory:

``` r
setwd("~/Documentos/pcamejo/")
```

Start importing counts table and metadata associated to the samples
(previously downloaded from
[Data](https://github.com/ibioChile/Transcriptomics-R-Workshop-public/tree/master/Session3-Treatment_and_Multivariate/Data)
folder).

``` r
counts <- read.table("fc0.original.counts.session2-2.txt",sep="\t",header=T,row.names = 1)

kable(head(counts))
```

|           | SRR5440784.trimmed.Athaliana.bam | SRR5440785.trimmed.Athaliana.bam | SRR5440786.trimmed.Athaliana.bam | SRR5440787.trimmed.Athaliana.bam | SRR5440788.trimmed.Athaliana.bam | SRR5440789.trimmed.Athaliana.bam | SRR5440790.trimmed.Athaliana.bam | SRR5440791.trimmed.Athaliana.bam | SRR5440792.trimmed.Athaliana.bam | SRR5440793.trimmed.Athaliana.bam | SRR5440794.trimmed.Athaliana.bam | SRR5440795.trimmed.Athaliana.bam | SRR5440796.trimmed.Athaliana.bam | SRR5440797.trimmed.Athaliana.bam | SRR5440798.trimmed.Athaliana.bam | SRR5440799.trimmed.Athaliana.bam | SRR5440800.trimmed.Athaliana.bam | SRR5440801.trimmed.Athaliana.bam | SRR5440802.trimmed.Athaliana.bam | SRR5440803.trimmed.Athaliana.bam | SRR5440804.trimmed.Athaliana.bam | SRR5440805.trimmed.Athaliana.bam | SRR5440806.trimmed.Athaliana.bam | SRR5440807.trimmed.Athaliana.bam | SRR5440808.trimmed.Athaliana.bam | SRR5440809.trimmed.Athaliana.bam | SRR5440810.trimmed.Athaliana.bam | SRR5440811.trimmed.Athaliana.bam | SRR5440812.trimmed.Athaliana.bam | SRR5440813.trimmed.Athaliana.bam | SRR5440814.trimmed.Athaliana.bam | SRR5440815.trimmed.Athaliana.bam | SRR5440816.trimmed.Athaliana.bam | SRR5440817.trimmed.Athaliana.bam | SRR5440818.trimmed.Athaliana.bam | SRR5440819.trimmed.Athaliana.bam | SRR5440820.trimmed.Athaliana.bam | SRR5440821.trimmed.Athaliana.bam | SRR5440822.trimmed.Athaliana.bam | SRR5440823.trimmed.Athaliana.bam | SRR5440824.trimmed.Athaliana.bam | SRR5440825.trimmed.Athaliana.bam | SRR5440826.trimmed.Athaliana.bam | SRR5440827.trimmed.Athaliana.bam | SRR5440828.trimmed.Athaliana.bam | SRR5440829.trimmed.Athaliana.bam | SRR5440830.trimmed.Athaliana.bam | SRR5440831.trimmed.Athaliana.bam | SRR5440832.trimmed.Athaliana.bam | SRR5440833.trimmed.Athaliana.bam | SRR5440834.trimmed.Athaliana.bam | SRR5440835.trimmed.Athaliana.bam | SRR5440836.trimmed.Athaliana.bam | SRR5440837.trimmed.Athaliana.bam | SRR5440838.trimmed.Athaliana.bam | SRR5440839.trimmed.Athaliana.bam | SRR5440840.trimmed.Athaliana.bam | SRR5440841.trimmed.Athaliana.bam | SRR5440842.trimmed.Athaliana.bam | SRR5440843.trimmed.Athaliana.bam | SRR5440844.trimmed.Athaliana.bam | SRR5440845.trimmed.Athaliana.bam | SRR5440846.trimmed.Athaliana.bam | SRR5440847.trimmed.Athaliana.bam | SRR5440848.trimmed.Athaliana.bam | SRR5440849.trimmed.Athaliana.bam | SRR5440850.trimmed.Athaliana.bam | SRR5440851.trimmed.Athaliana.bam | SRR5440852.trimmed.Athaliana.bam | SRR5440853.trimmed.Athaliana.bam | SRR5440854.trimmed.Athaliana.bam | SRR5440855.trimmed.Athaliana.bam | SRR5440856.trimmed.Athaliana.bam | SRR5440857.trimmed.Athaliana.bam | SRR5440858.trimmed.Athaliana.bam | SRR5440859.trimmed.Athaliana.bam | SRR5440860.trimmed.Athaliana.bam | SRR5440861.trimmed.Athaliana.bam | SRR5440862.trimmed.Athaliana.bam | SRR5440863.trimmed.Athaliana.bam | SRR5440864.trimmed.Athaliana.bam | SRR5440865.trimmed.Athaliana.bam | SRR5440866.trimmed.Athaliana.bam | SRR5440867.trimmed.Athaliana.bam | SRR5440868.trimmed.Athaliana.bam | SRR5440869.trimmed.Athaliana.bam | SRR5440870.trimmed.Athaliana.bam | SRR5440871.trimmed.Athaliana.bam | SRR5440872.trimmed.Athaliana.bam | SRR5440873.trimmed.Athaliana.bam | SRR5440874.trimmed.Athaliana.bam | SRR5440875.trimmed.Athaliana.bam | SRR5440876.trimmed.Athaliana.bam | SRR5440877.trimmed.Athaliana.bam | SRR5440878.trimmed.Athaliana.bam | SRR5440879.trimmed.Athaliana.bam | SRR5440880.trimmed.Athaliana.bam | SRR5440881.trimmed.Athaliana.bam | SRR5440882.trimmed.Athaliana.bam | SRR5440883.trimmed.Athaliana.bam | SRR5440884.trimmed.Athaliana.bam | SRR5440885.trimmed.Athaliana.bam | SRR5440886.trimmed.Athaliana.bam | SRR5440887.trimmed.Athaliana.bam | SRR5440888.trimmed.Athaliana.bam | SRR5440889.trimmed.Athaliana.bam | SRR5440890.trimmed.Athaliana.bam | SRR5440891.trimmed.Athaliana.bam | SRR5440892.trimmed.Athaliana.bam | SRR5440893.trimmed.Athaliana.bam | SRR5440894.trimmed.Athaliana.bam | SRR5440895.trimmed.Athaliana.bam | SRR5440896.trimmed.Athaliana.bam | SRR5440897.trimmed.Athaliana.bam |
| --------- | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: | -------------------------------: |
| AT1G01010 |                             1040 |                             1512 |                             1162 |                              373 |                              587 |                              157 |                              588 |                              339 |                              508 |                              458 |                              615 |                              568 |                              258 |                              355 |                              412 |                              555 |                              261 |                              643 |                              699 |                                0 |                              692 |                              384 |                              373 |                              268 |                              445 |                              321 |                              548 |                              556 |                              423 |                              475 |                              517 |                              379 |                              440 |                              279 |                              384 |                              362 |                              421 |                              221 |                              265 |                              339 |                              406 |                              597 |                              394 |                              367 |                              263 |                              307 |                              379 |                              283 |                              369 |                              296 |                              588 |                              258 |                              383 |                              265 |                              400 |                              199 |                              487 |                             2697 |                             2746 |                             1650 |                              784 |                              333 |                              595 |                             1278 |                             1084 |                              960 |                              901 |                             1296 |                              842 |                              702 |                              791 |                              384 |                              698 |                             1220 |                             1064 |                              390 |                             1033 |                              979 |                              927 |                              599 |                              585 |                              961 |                              944 |                              922 |                              937 |                              843 |                              903 |                             1081 |                              471 |                              216 |                              641 |                             1015 |                              428 |                              739 |                              667 |                             1393 |                              868 |                              508 |                              258 |                              867 |                              614 |                              596 |                              631 |                              779 |                              869 |                             1234 |                              459 |                              271 |                              653 |                              609 |                              525 |                              475 |                              586 |                             1851 |
| AT1G01020 |                              376 |                              722 |                              679 |                              215 |                              180 |                               64 |                              162 |                              157 |                              249 |                              253 |                              293 |                              226 |                              217 |                              246 |                              172 |                              337 |                              153 |                              267 |                              283 |                                0 |                              240 |                              154 |                              191 |                              195 |                              246 |                              209 |                              186 |                              235 |                              242 |                              220 |                              206 |                              226 |                              213 |                              291 |                              244 |                              224 |                              209 |                              327 |                              260 |                              260 |                              197 |                              368 |                              217 |                              168 |                              171 |                              178 |                              340 |                              277 |                              225 |                              206 |                              282 |                              103 |                              222 |                              243 |                              267 |                              201 |                              370 |                             2151 |                             1381 |                             1226 |                              649 |                               76 |                               75 |                              283 |                              245 |                              233 |                              281 |                              230 |                              130 |                              835 |                               78 |                               29 |                              263 |                              212 |                              281 |                               93 |                              311 |                              209 |                              651 |                               66 |                               41 |                              332 |                              366 |                              277 |                              200 |                              368 |                              335 |                              484 |                               41 |                               24 |                              228 |                              270 |                              179 |                              260 |                              191 |                              666 |                              640 |                               53 |                               21 |                              371 |                              319 |                              425 |                              325 |                              275 |                              435 |                              406 |                               49 |                               25 |                              268 |                              205 |                              299 |                              358 |                              181 |                             1137 |
| AT1G03987 |                                1 |                                0 |                                4 |                                1 |                                1 |                                0 |                                1 |                                1 |                                1 |                                1 |                                3 |                                2 |                                1 |                                1 |                                0 |                                1 |                                1 |                                0 |                                0 |                                0 |                                1 |                                0 |                                0 |                                0 |                                0 |                                2 |                                0 |                                0 |                                0 |                                0 |                                0 |                                1 |                                0 |                                1 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                1 |                                1 |                                0 |                                1 |                                0 |                                0 |                                0 |                                0 |                                1 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                1 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                1 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                1 |                                1 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                1 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                1 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                1 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |                                0 |
| AT1G01030 |                              116 |                              195 |                              179 |                               73 |                               41 |                               15 |                               55 |                               55 |                               60 |                               54 |                               70 |                               81 |                               58 |                               47 |                               49 |                               78 |                               49 |                               56 |                               69 |                                0 |                               69 |                               53 |                               56 |                               59 |                               53 |                               58 |                               55 |                               76 |                               66 |                               62 |                               51 |                               60 |                               55 |                               48 |                               55 |                               37 |                               60 |                               65 |                               45 |                               77 |                               56 |                               81 |                               32 |                               42 |                               32 |                               33 |                               58 |                               49 |                               85 |                               60 |                               78 |                               50 |                               50 |                               37 |                               43 |                               23 |                               61 |                              173 |                               82 |                               94 |                               47 |                               22 |                               30 |                               26 |                               24 |                               17 |                               26 |                               25 |                               12 |                               59 |                               34 |                               23 |                               18 |                               20 |                               25 |                                3 |                               14 |                                9 |                               53 |                               34 |                               24 |                               21 |                               35 |                               24 |                               17 |                               23 |                               17 |                               30 |                               33 |                               11 |                               12 |                               18 |                               12 |                               19 |                               13 |                               38 |                               70 |                               33 |                               30 |                               12 |                               18 |                               21 |                               16 |                               21 |                               23 |                               24 |                               20 |                               22 |                               25 |                               17 |                               18 |                               22 |                               15 |                               73 |
| AT1G01040 |                             3931 |                             6002 |                             5513 |                             1482 |                             1784 |                              532 |                             1622 |                             1461 |                             1800 |                             1739 |                             2211 |                             1855 |                             1585 |                             1659 |                             1446 |                             2195 |                             1214 |                             1805 |                             2389 |                                0 |                             2295 |                             1247 |                             1415 |                             1438 |                             1746 |                             1519 |                             1605 |                             1865 |                             1908 |                             1744 |                             1418 |                             1503 |                             1436 |                             1415 |                             1787 |                             1571 |                             2006 |                             1340 |                             1413 |                             1702 |                             1382 |                             2453 |                             1369 |                             1388 |                             1365 |                             1257 |                             1622 |                             1458 |                             1737 |                             1384 |                             2240 |                              983 |                             1627 |                             1346 |                             1833 |                             1174 |                             2119 |                             4606 |                             3944 |                             3058 |                             1555 |                             1185 |                             1157 |                             1320 |                             1025 |                             1165 |                              980 |                             1062 |                              962 |                             1633 |                             1567 |                             1055 |                              850 |                             1149 |                             1361 |                              644 |                             1369 |                             1222 |                             1634 |                             1528 |                             1290 |                             1138 |                             1202 |                             1219 |                             1165 |                             1223 |                             1364 |                             1524 |                             1061 |                              570 |                              847 |                              884 |                              723 |                             1091 |                              786 |                             1741 |                             1336 |                             1287 |                             1024 |                             1123 |                             1085 |                              883 |                             1127 |                              964 |                             1083 |                             1251 |                             1054 |                              727 |                              875 |                              860 |                              893 |                              629 |                              780 |                             2564 |
| AT1G03993 |                              327 |                              564 |                              513 |                              102 |                              140 |                               32 |                               99 |                               99 |                              151 |                              156 |                              209 |                              142 |                              109 |                              131 |                              116 |                              160 |                               94 |                              169 |                              194 |                                0 |                              188 |                               93 |                              110 |                               82 |                              159 |                              104 |                              154 |                              168 |                              159 |                              143 |                               94 |                               94 |                               81 |                               38 |                              106 |                              109 |                              155 |                               65 |                               82 |                              111 |                              116 |                              182 |                               75 |                               65 |                              122 |                               65 |                              101 |                               98 |                              105 |                               67 |                              151 |                               89 |                              108 |                               73 |                              143 |                               80 |                              195 |                              144 |                              176 |                              105 |                               78 |                               24 |                               50 |                               69 |                               42 |                               93 |                               42 |                               55 |                               94 |                              117 |                               64 |                               66 |                               33 |                              102 |                               74 |                               39 |                               92 |                               92 |                               75 |                               65 |                               63 |                               39 |                               30 |                               69 |                               94 |                               67 |                               69 |                               61 |                               49 |                               20 |                               27 |                               48 |                               44 |                               54 |                               62 |                               57 |                               78 |                               38 |                               90 |                               37 |                               28 |                               31 |                               33 |                               56 |                               41 |                               45 |                               50 |                               27 |                               30 |                               22 |                               49 |                               25 |                               48 |                               89 |

Normalization of gene reads counts from RNA-seq data.

Reads counts for each gene can be used to calculate correlation of gene
expression . In this example, read counts should be normalized to the
number of reads that were effectively mapped. In this example, we use
median normalization of the data. This method can be applied using the
EBSeq R package. EBSeq generates a normalized matrix using different
methods including median normalization. The normalization will be
applied to all the mapped reads selected in previous sections.

The following commands can perform the mentioned normalization. First
install EBSeq:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EBSeq")
```

    ## 
    ## The downloaded binary packages are in
    ##  /var/folders/mq/q_8_31_564376clx2y7s1g4c0000gn/T//RtmpFMWAak/downloaded_packages

``` r
library(EBSeq)
```

Then normalize the data in “counts” file executing EBSeq and store the
normalized data in “NormData” object:

``` r
NormData <- GetNormalizedMat(counts, MedianNorm(counts))
```

The resulting gene expression matrix “NormData” contains unique row
identifiers and row counts obtained from different experiments from the
previous steo in each column. After normalization, it is recommended to
delete very low counts data or sum a unit to all data in order to avoid
values equal to zero. Then is useful to generate a logarithmic matrix of
the data to standardize the variance. All previous recommendations can
be done by simply:

``` r
NormData.log <- log2(NormData+1)
```

The “NormData.log” object stores normalized counts in logarithmic scale
from all libraries.

\#\#Calculating correlation of gene expression for each pair of genes

We are now getting close to the goal of this guide which is to generate
a gene co-expression network. In this step, we will determine
correlation of every pair of DEGs across the complete data set. First,
the user should select DEGs obtained by DESeq2 from the “NormData.log”
table generated previously and extract their normalized counts, using
the following
command:

``` r
DEGs.deseq2<-read.table("regulated_DESEQ2_log2FC1_padj0.01.txt",header=T,row.names=1)

select<-DEGs.deseq2[abs(DEGs.deseq2$log2FoldChange)>1,]

# For this example, we will use a stricter criterion and we will select the genes that change at least 4 times with respect to the control. The number of genes selected will influence the calculation time of the correlation.

Norm.interest <- NormData.log[rownames(select),]

dim(Norm.interest)
```

    ## [1] 2536  114

The “Norm.interest” object contains the normalized counts for DEGs. In
order to build the gene co-expression network, the user needs to
calculate correlation and correlation significance between each pair of
DEGs in the whole experimental data set. There are several tools to
calculate correlation, here we use the R package “psych” (
(<https://www.rdocumentation.org/packages/psych/>). Psych is a
general-purpose toolbox originally developed for psychometrics analyses,
but with useful tools for data analyses including correlation analysis.

The user can download the package from the Comprehensive R Archive
Network (CRAN) using the following commands:

``` r
#If it is not installed, install with:
#install.packages("mnormt")
install.packages("https://cran.r-project.org/src/contrib/Archive/psych/psych_1.7.5.tar.gz", repos=NULL, type="source")

library("psych")
```

Then, set the instructions for the correlation calculation for every
DEGs pairs. Between the available methods select Pearson correlation to
analyze the normalized data. The correlation results will be stored in
“Norm.interest.corr”
object:

``` r
Norm.interest.corr <- corr.test( t(Norm.interest), method="pearson", ci=F)
```

Among the many results of this function, there are two triangular
matrices with the needed data. One matrix is symmetric and contains the
correlation values. The other matrix contains the p-values in the upper
part and the adjusted p-values in the lower part. To generate a table
comprising the data organized properly to visualize the network, the
user should execute the following commands:

``` r
Norm.interest.corr$p[lower.tri( Norm.interest.corr$p,diag=TRUE)]=NA
Pval.adj <- as.data.frame(as.table(Norm.interest.corr$p))
#The "Pval.adj" object contains all p-values from the lower part of the matrix.

Norm.interest.corr$r [lower.tri( Norm.interest.corr$r,diag=TRUE)]=NA
Correlation <- as.data.frame(as.table(Norm.interest.corr$r))

#The "Correlation" object contains all correlations from the matrix.

Cor.table <- na.exclude(cbind( Correlation, Pval.adj))[,c(1,2,3,6)]

#The "Cor.table" object contains all correlations and p-values selected from the matrix. The following command adds the names to the columns in "Cor.table":

colnames(Cor.table) <- c("gene1","gene2","cor","p.adj")

kable(head(Cor.table))
```

|      | gene1     | gene2     |         cor | p.adj |
| ---- | :-------- | :-------- | ----------: | ----: |
| 2537 | AT1G01140 | AT1G01180 |   0.3983536 |     1 |
| 5073 | AT1G01140 | AT1G01210 | \-0.1452497 |     1 |
| 5074 | AT1G01180 | AT1G01210 |   0.4082665 |     1 |
| 7609 | AT1G01140 | AT1G04023 |   0.6753469 |     0 |
| 7610 | AT1G01180 | AT1G04023 |   0.2361129 |     1 |
| 7611 | AT1G01210 | AT1G04023 | \-0.0954803 |     1 |

The generated “Cor.table” object, can be filtered based on absolute
correlation (0.9) and adjusted p-value (0.01)
thresholds:

``` r
Cor.table.filt <- Cor.table [(abs(Cor.table[,3])>0.7 & Cor.table[,4] <0.01 ),]
write.table(Cor.table.filt, "Cor.table.filter.txt", sep="\t", row.names=F, quote=F)

dim(Cor.table.filt)
```

    ## [1] 633795      4

At this point, we generated the “Cor.table.filter.txt” file containing
the statistically significant correlations across the whole data set for
every pair of differentially expressed genes.

For this example, we will create a transcriptional network using a
TF-Target database. We are going to cross the pairs that have a major
correlation to 0.9 with the information obtained from the DAP-seq
database (<http://neomorph.salk.edu/dap_web/pages/browse_table_aj.php>),
which contains information from TF-TARGET pairs with evidence of
possible binding.

Para este ejemplo utilizaremos un subset de esta base de datos, que la
podemos descargar del github del curso (dapseq.subset.txt).

``` r
dapseq<-read.table("dapseq.subset.txt.gz",header=T)
kable(head(dapseq))
```

| TF.id     | target    | source |
| :-------- | :-------- | :----- |
| AT1G01060 | AT1G01010 | DAPSEQ |
| AT1G01060 | AT1G01020 | DAPSEQ |
| AT1G01060 | AT1G01030 | DAPSEQ |
| AT1G01060 | AT1G01040 | DAPSEQ |
| AT1G01060 | AT1G01073 | DAPSEQ |
| AT1G01060 | AT1G01140 | DAPSEQ |

Once loaded the TF-TARGET table we intersect this table with the
correlation data and write to a
file.

``` r
output.dapseq.filt<- unique(rbind(merge(dapseq,Cor.table.filt,by.x=c(1,2),by.y=c(1,2)),merge(dapseq,Cor.table.filt,by.x=c(1,2),by.y=c(2,1))))



write.table(output.dapseq.filt,"output.dapseq.filt.txt", sep="\t",quote=F,row.names = F)

kable(head(output.dapseq.filt))
```

| TF.id     | target    | source |       cor | p.adj |
| :-------- | :-------- | :----- | --------: | ----: |
| AT1G09540 | AT1G21890 | DAPSEQ | 0.7022179 |     0 |
| AT1G09540 | AT1G51680 | DAPSEQ | 0.7076866 |     0 |
| AT1G09540 | AT1G58170 | DAPSEQ | 0.7848574 |     0 |
| AT1G09540 | AT1G61820 | DAPSEQ | 0.8098619 |     0 |
| AT1G09540 | AT2G15230 | DAPSEQ | 0.8595091 |     0 |
| AT1G09540 | AT3G05390 | DAPSEQ | 0.8311129 |     0 |

With this instrucction we will identify the transcription factor used in
the network.

``` r
TFs<-unique(output.dapseq.filt[,1])
TFs<-cbind(TFs,"TF")
colnames(TFs)<-c("id","TF")
write.table(TFs,"TFs.txt",sep="\t",col.names=NA,quote=F)
```

##### Network Visualization

Co-expression networks help associate genes that are involved in similar
biological functions. The analysis and visualization of gene networks is
a key and powerful step to identify relationships and discover important
elements in the network. Analysis of gene networks also offers us the
opportunity to formulate hypotheses about key genes and implicated
biological functions. Below, we describe some simple steps to visualize
a gene co-expression network.

Before we generate a network view, we will calculate a few useful
network statistics. The user can later graphically represent network
statistics on the same network. For instance, in Cytoscape (a popular
software platform to view and analyze networks), statistics can be
calculated and added as attributes to the nodes or edges. The number of
connections of a node in a network is known as “degree”, and it is a
useful statistic to identify relevant nodes in biological networks,
which are typically highly connected nodes or hubs. Another important
node statistic is the number of times that a path passes through the
node, which represents the amount of control that this node exerts over
the interactions of other nodes in the networks known as “Betweenness
centrality”. To calculate network statistics in this guide we used
“igraph” (<http://igraph.org> ). This R package can process simple
graphs and network analysis handling large graphs if necessary. In
addition, it offers functions for generating graph visualization and
complete network statistics.

The basic statistics of the network, degree and betweenness, can be
calculated using “igraph” R package. Download and install “igraph” with
following commands:

``` r
#If it is not installed, install with:
#install.packages("igraph")
library(igraph)
```

Select the columns 1 and 2 in the “Cor.table.filt” object. These columns
contain all the DEGs pairs with high correlation. Store them in “g”
object, which is the network:

``` r
g <- graph.data.frame( output.dapseq.filt, directed=TRUE)
g
```

    ## IGRAPH aa28f2c DN-- 1237 2198 -- 
    ## + attr: name (v/c), source (e/c), cor (e/n), p.adj (e/n)
    ## + edges from aa28f2c (vertex names):
    ##  [1] AT1G09540->AT1G21890 AT1G09540->AT1G51680 AT1G09540->AT1G58170
    ##  [4] AT1G09540->AT1G61820 AT1G09540->AT2G15230 AT1G09540->AT3G05390
    ##  [7] AT1G09540->AT3G21240 AT1G09540->AT3G49190 AT1G09540->AT3G49930
    ## [10] AT1G09540->AT3G63010 AT1G09540->AT4G31100 AT1G09540->AT5G07010
    ## [13] AT1G09540->AT5G19080 AT1G09540->AT5G27000 AT1G09540->AT5G47560
    ## [16] AT1G09540->AT5G56610 AT1G12630->AT1G13300 AT1G12630->AT1G16300
    ## [19] AT1G12630->AT1G19700 AT1G12630->AT1G24290 AT1G12630->AT1G27050
    ## [22] AT1G12630->AT1G35560 AT1G12630->AT1G37130 AT1G12630->AT1G49230
    ## + ... omitted several edges

Calculate the degree and betweenness for each node in “g” object and
store the results in the “degree” and “betweenness” objects
respectively:

``` r
degree <- degree(g)
betweenness <- betweenness(g)
Node_nw_st <- data.frame( degree, betweenness)
```

The “Node\_nw\_st” object contains all the calculated statistics for
each node in “g” object. To integrate the two parameters, we generated
two different ranking for degree and betweenness for all nodes in
“Node\_nw\_st”, and then calculated the mean of the rankings at each
node. The nodes with the higher degree and betweenness will be scored
with the higher value. The combined ranking can be added to the former
table:

``` r
Rank_stat <- rowMeans(cbind(rank(Node_nw_st[,1]), rank(Node_nw_st[,2])))
Node_nw_st <- cbind(Node_nw_st, Rank_stat)
write.table(Node_nw_st,file="Node_nw_st.txt", sep="\t", col.names = NA, quote=F)
kable(head(Node_nw_st))
```

|           | degree | betweenness | Rank\_stat |
| --------- | -----: | ----------: | ---------: |
| AT1G09540 |     18 |        41.5 |    1214.00 |
| AT1G12630 |     52 |         0.0 |     912.25 |
| AT1G25550 |     77 |       158.0 |    1221.50 |
| AT1G29860 |     73 |       192.0 |    1222.00 |
| AT1G32870 |     84 |         0.0 |     916.25 |
| AT1G44830 |      8 |        12.0 |    1212.25 |

The “Node\_nw\_st.txt” file contains all the calculated statistics for
each node. These statistics will be used to complement the visualization
of the “output.dapseq.filt.txt” network created previously.

We will use Cytoscape to visualize the network and corresponding
statistics. Cytoscape is an open source software platform. Cytoscape can
be downloaded from <http://www.cytoscape.org/> and requires JAVA™ JRE or
JDK.

First, we launch Cytoscape, and then we import the network table
“output.dapseq.filt.txt”. This can be done by selecting File \> Import
\> Network \> File. After selecting the file, the user should indicate
where relevant information is stored in the file. In our case, “Source
Interaction” is in the first column which is labeled “TF” and “Target
Interaction” is in the second column which is labeled “Target”. These
columns contain the gene ID information that Cytoscape will use to
identify the nodes in each interacting pair in the network.

It is useful to keep the information contained in the table shown in the
“Preview” window as an edge, (e.g. correlation value, adjusted
p-values), the user should click the corresponding column header to
activate it. After the network has been displayed, and to get a better
visualization, nodes can be arranged following different layouts.

One of the most common displays is the “organic layout” which can be
found under the Layout \> yFiles Layout \> Organic menu within
Cytoscape. This layout will display nodes based in repulsive forces
between the nodes and attractive forces induced by edges (55),
facilitating identification of highly connected nodes.

Since both degree and betweenness centrality are measures of the
function of a node in network connectivity, we will use the combined
ranking in the node attribute file “Node\_nw\_st.txt” generated early in
this section to map the size of the nodes to represent node importance,
the file is located in the “WorkDir” folder.

To load the node information calculated before, the user should go
through File \> Import \> Table \> File and select “Node\_nw\_st.txt”.
In the pop-up window, the dropdown list “Network Collection”, select the
imported network, in this case “output.dapseq.filt.txt”. Import the
table as “Node table columns” in the “Import Data as:” dropdown list and
be sure all columns are checked. Node size can be set accordingly by
selecting the “Style” tab under the “Control Panel” box. The “Style” tab
also allows the user to set graphic properties on edges and the whole
network in their corresponding sections in the lower part of the box.

Under the “Node” section in the “Style” tab of the “Control Panel” box,
set the size of the node by checking “Lock node width and height”, and
then select “Size”. On the new menu displayed, click on “Column” to show
a new menu in which select “Rank\_stat” as the attribute to determine
the size. “Mapping type” should be in “Continuous Mapping” option to
distribute the sizes continuously along the previously generated
combined ranking for each node . The size variation can be adjusted in
the graph that appears in the “Current Mapping” section.

The edges appearance can also be customized in a similar way to nodes.
In the same “Style” tab in the “Control Panel”, we select the “Edge” tab
in the lower part of the section. Click the option “Edge color to
arrows”. In the option “Column” select “cor” as the mapping attribute.
Then in “Mapping Type” select “Continuous Mapping” option. Clicking the
graph is possible to adjust the colors and intensities. The user must
set the minimal and maximum values every time, in this case -1 and 1
respectively.

\#Performing cluster analysis of the network

Analysis of network topology or the structure that determine the way in
which genes are connected is useful to derive biological insights. For
instance, subnetworks of genes that are highly connected (network
cluster) in a co-expression network are usually involved in similar
biological functions. To find groups of genes that may be acting in a
coordinated manner, the user can perform a cluster analysis of the
network.

Numerous network clustering apps can be found at Cytoscape Application
Store web (<http://> apps.cytoscape.org/). For the purpose of this
example, download, install and use “clusterMaker” plugin. clusterMaker
offers many options to perform cluster analysis that can be addressed in
its manual. We use “Community Clustering” (GLay) with the default
options, because it provides an optimized layout for large networks and
a structured visualization for more efficient exploration and analysis
of biological networks . To display the clustering results as a new
network we select “Create new clustered network”, and we check “Assume
edges are undirected”. In this case, directionality cannot be assumed,
because relationships were established based on correlation of gene
expression. To keep the gene-gene connection we select “Restore
inter-cluster edges after layout”. This analysis generates a new network
which contains the same nodes than the previous but arranged according
to the connectivity of the nodes in the different clusters.
