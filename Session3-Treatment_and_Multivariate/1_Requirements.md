# Session 3 | Requirements

### Pipeline 1

For this session, we will use data published in the paper ["Temporal transcriptional logic of dynamic regulatory
networks underlying nitrogen signaling and use
in plants"](https://www.pnas.org/content/pnas/115/25/6494.full.pdf). This time, we will use RNA sequencing data from samples of Arabidopsis thaliana's shoot taken at time 120 min, growing with 5mm of nitrate (no treatment) and no nitrate addition (KCl treatment).

We will not use the original files for this session, but if you'd like to start the analysis from the beginning, the original files and can be found here:

    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR544/FOLDER/*
     
  |FOLDER |Tissue|Treatment|Replicate|Time|
  |----------------|-------|------|---|-----|
  | ```000/SRR5440840``` | Shoot | None | 1 | 120 |
  | ```001/SRR5440831``` | Shoot | None | 2 | 120 |
  | ```002/SRR5440822```| Shoot | None | 3 | 120 |
  | ```003/SRR5440813``` | Shoot | KCl  | 1 | 120 |
  | ```004/SRR5440804``` | Shoot | KCl  | 2 | 120 |
  | ```005/SRR5440795``` | Shoot | KCl  | 3 | 120 |

Following the **Data Preprocessing** tutorial from Session 1, we generated a counts table from these files. You can download it [here](https://github.com/ibioChile/Transcriptomics-R-Workshop-public/blob/master/Session3-Treatment_and_Multivariate/Data/fc0.original.counts.session2-1.txt). 

The associated metadata can be found [here](https://github.com/ibioChile/Transcriptomics-R-Workshop-public/blob/master/Session3-Treatment_and_Multivariate/Data/metadata_session2-1.txt).

Download Arabidopsis Gene Ontology data:

```wget http://current.geneontology.org/annotations/tair.gaf.gz```

### Pipeline 2

For this pipeline, we will use RNA sequencing data from samples of Arabidopsis thaliana's shoot and root, taken at different times and growing with 5mm of nitrate and no nitrate addition. 

The original files and can be found here:

    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR544/FOLDER/*
    
|FOLDER |Tissue|Treatment|Replicate|Time|
|----------------------|-------|------|---|-----|
| ```004/SRR5440784``` | Shoot | None | 3 | 0   |
| ```005/SRR5440785``` | Shoot | None | 2 | 0   |
| ```006/SRR5440786``` | Shoot | None | 1 | 0   |
| ```007/SRR5440787``` | Shoot | KCl  | 3 | 5   |
| ```008/SRR5440788``` | Shoot | KCl  | 3 | 10  |
| ```009/SRR5440789``` | Shoot | KCl  | 3 | 15  |
| ```000/SRR5440790``` | Shoot | KCl  | 3 | 20  |
| ```001/SRR5440791``` | Shoot | KCl  | 3 | 30  |
| ```002/SRR5440792``` | Shoot | KCl  | 3 | 45  |
| ```003/SRR5440793``` | Shoot | KCl  | 3 | 60  |
| ```004/SRR5440794``` | Shoot | KCl  | 3 | 90  |
| ```005/SRR5440795``` | Shoot | KCl  | 3 | 120 |
| ```006/SRR5440796``` | Shoot | KCl  | 2 | 5   |
| ```007/SRR5440797``` | Shoot | KCl  | 2 | 10  |
| ```008/SRR5440798``` | Shoot | KCl  | 2 | 15  |
| ```009/SRR5440799``` | Shoot | KCl  | 2 | 20  |
| ```000/SRR5440800``` | Shoot | KCl  | 2 | 30  |
| ```001/SRR5440801``` | Shoot | KCl  | 2 | 45  |
| ```002/SRR5440802``` | Shoot | KCl  | 2 | 60  |
| ```003/SRR5440803``` | Shoot | KCl  | 2 | 90  |
| ```004/SRR5440804``` | Shoot | KCl  | 2 | 120 |
| ```005/SRR5440805``` | Shoot | KCl  | 1 | 5   |
| ```006/SRR5440806``` | Shoot | KCl  | 1 | 10  |
| ```007/SRR5440807``` | Shoot | KCl  | 1 | 15  |
| ```008/SRR5440808``` | Shoot | KCl  | 1 | 20  |
| ```009/SRR5440809``` | Shoot | KCl  | 1 | 30  |
| ```000/SRR5440810``` | Shoot | KCl  | 1 | 45  |
| ```001/SRR5440811``` | Shoot | KCl  | 1 | 60  |
| ```002/SRR5440812``` | Shoot | KCl  | 1 | 90  |
| ```003/SRR5440813``` | Shoot | KCl  | 1 | 120 |
| ```004/SRR5440814``` | Shoot | None | 1 | 5   |
| ```005/SRR5440815``` | Shoot | None | 1 | 10  |
| ```006/SRR5440816``` | Shoot | None | 3 | 15  |
| ```007/SRR5440817``` | Shoot | None | 3 | 20  |
| ```008/SRR5440818``` | Shoot | None | 3 | 30  |
| ```009/SRR5440819``` | Shoot | None | 3 | 45  |
| ```000/SRR5440820``` | Shoot | None | 3 | 60  |
| ```001/SRR5440821``` | Shoot | None | 3 | 90  |
| ```002/SRR5440822``` | Shoot | None | 3 | 120 |
| ```003/SRR5440823``` | Shoot | None | 3 | 5   |
| ```004/SRR5440824``` | Shoot | None | 3 | 10  |
| ```005/SRR5440825``` | Shoot | None | 2 | 15  |
| ```006/SRR5440826``` | Shoot | None | 2 | 20  |
| ```007/SRR5440827``` | Shoot | None | 2 | 30  |
| ```008/SRR5440828``` | Shoot | None | 2 | 45  |
| ```009/SRR5440829``` | Shoot | None | 2 | 60  |
| ```000/SRR5440830``` | Shoot | None | 2 | 90  |
| ```001/SRR5440831``` | Shoot | None | 2 | 120 |
| ```002/SRR5440832``` | Shoot | None | 2 | 5   |
| ```003/SRR5440833``` | Shoot | None | 2 | 10  |
| ```004/SRR5440834``` | Shoot | None | 1 | 15  |
| ```005/SRR5440835``` | Shoot | None | 1 | 20  |
| ```006/SRR5440836``` | Shoot | None | 1 | 30  |
| ```007/SRR5440837``` | Shoot | None | 1 | 45  |
| ```008/SRR5440838``` | Shoot | None | 1 | 60  |
| ```009/SRR5440839``` | Shoot | None | 1 | 90  |
| ```000/SRR5440840``` | Shoot | None | 1 | 120 |
| ```001/SRR5440841``` | Root  | None | 1 | 0   |
| ```002/SRR5440842``` | Root  | None | 3 | 0   |
| ```003/SRR5440843``` | Root  | None | 2 | 0   |
| ```004/SRR5440844``` | Root  | KCl  | 3 | 5   |
| ```005/SRR5440845``` | Root  | KCl  | 3 | 10  |
| ```006/SRR5440846``` | Root  | KCl  | 3 | 15  |
| ```007/SRR5440847``` | Root  | KCl  | 3 | 20  |
| ```008/SRR5440848``` | Root  | KCl  | 3 | 30  |
| ```009/SRR5440849``` | Root  | KCl  | 3 | 45  |
| ```000/SRR5440850``` | Root  | KCl  | 3 | 60  |
| ```001/SRR5440851``` | Root  | KCl  | 3 | 90  |
| ```002/SRR5440852``` | Root  | KCl  | 3 | 120 |
| ```003/SRR5440853``` | Root  | KCl  | 2 | 5   |
| ```004/SRR5440854``` | Root  | KCl  | 2 | 10  |
| ```005/SRR5440855``` | Root  | KCl  | 2 | 15  |
| ```006/SRR5440856``` | Root  | KCl  | 2 | 20  |
| ```007/SRR5440857``` | Root  | KCl  | 2 | 30  |
| ```008/SRR5440858``` | Root  | KCl  | 2 | 45  |
| ```009/SRR5440859``` | Root  | KCl  | 2 | 60  |
| ```000/SRR5440860``` | Root  | KCl  | 2 | 90  |
| ```001/SRR5440861``` | Root  | KCl  | 2 | 120 |
| ```002/SRR5440862``` | Root  | KCl  | 1 | 5   |
| ```003/SRR5440863``` | Root  | KCl  | 1 | 10  |
| ```004/SRR5440864``` | Root  | KCl  | 1 | 15  |
| ```005/SRR5440865``` | Root  | KCl  | 1 | 20  |
| ```006/SRR5440866``` | Root  | KCl  | 1 | 30  |
| ```007/SRR5440867``` | Root  | KCl  | 1 | 45  |
| ```008/SRR5440868``` | Root  | KCl  | 1 | 60  |
| ```009/SRR5440869``` | Root  | KCl  | 1 | 90  |
| ```000/SRR5440870``` | Root  | KCl  | 1 | 120 |
| ```001/SRR5440871``` | Root  | None | 3 | 5   |
| ```002/SRR5440872``` | Root  | None | 3 | 10  |
| ```003/SRR5440873``` | Root  | None | 3 | 15  |
| ```004/SRR5440874``` | Root  | None | 3 | 20  |
| ```005/SRR5440875``` | Root  | None | 3 | 30  |
| ```006/SRR5440876``` | Root  | None | 3 | 45  |
| ```007/SRR5440877``` | Root  | None | 3 | 60  |
| ```008/SRR5440878``` | Root  | None | 3 | 90  |
| ```009/SRR5440879``` | Root  | None | 3 | 120 |
| ```000/SRR5440880``` | Root  | None | 2 | 5   |
| ```001/SRR5440881``` | Root  | None | 2 | 10  |
| ```002/SRR5440882``` | Root  | None | 2 | 15  |
| ```003/SRR5440883``` | Root  | None | 2 | 20  |
| ```004/SRR5440884``` | Root  | None | 2 | 30  |
| ```005/SRR5440885``` | Root  | None | 2 | 45  |
| ```006/SRR5440886``` | Root  | None | 2 | 60  |
| ```007/SRR5440887``` | Root  | None | 2 | 90  |
| ```008/SRR5440888``` | Root  | None | 2 | 120 |
| ```009/SRR5440889``` | Root  | None | 1 | 5   |
| ```000/SRR5440890``` | Root  | None | 1 | 10  |
| ```001/SRR5440891``` | Root  | None | 1 | 15  |
| ```002/SRR5440892``` | Root  | None | 1 | 20  |
| ```003/SRR5440893``` | Root  | None | 1 | 30  |
| ```004/SRR5440894``` | Root  | None | 1 | 45  |
| ```005/SRR5440895``` | Root  | None | 1 | 60  |
| ```006/SRR5440896``` | Root  | None | 1 | 90  |
| ```007/SRR5440897``` | Root  | None | 1 | 120 |


Following the **Data Preprocessing** tutorial from Session 1, we generated a counts table located [here](https://github.com/ibioChile/Transcriptomics-R-Workshop-public/blob/master/Session3-Treatment_and_Multivariate/Data/fc0.original.counts.session2-2.txt). 

The associated metadata can be found [here](https://github.com/ibioChile/Transcriptomics-R-Workshop-public/blob/master/Session3-Treatment_and_Multivariate/Data/metadata_session2-2.txt).

### Install ViSEAGO in R versions = 3.6.1 (or any incompatible version)

In your terminal, clone ViSEAGO repository:

    git clone https://forgemia.inra.fr/umr-boa/viseago.git
    
From R console:

    >BiocManager::install(c("BiocStyle","heatmaply","plotly","webshot","GOSemSim","DiagrammeR"))
    # build package 
    >devtools::build("/Users/pamelacamejo/viseago/") # Use path where you clone the ViSEAGO repo
    # install package
    >install.packages("/Users/pamelacamejo/ViSEAGO_1.3.5.tar.gz", repos = NULL, type = "source") # Use path where you build package
    >library(ViSEAGO)

    


