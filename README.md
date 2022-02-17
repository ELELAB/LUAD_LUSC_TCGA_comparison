Computational Biology Laboratory, Danish Cancer Society Research Center, Strandboulevarden 49, 2100, Copenhagen, Denmark

Repository associated to the publication:

Distinct signatures of lung cancer types: aberrant mucin O-glycosylation and compromised immune response.
Lucchetta M, da Piedade I, Mounir M, Vabistsevits M, Terkelsen T, Papaleo E*.
BMC Cancer. 2019 Aug 20;19(1):824. doi: 10.1186/s12885-019-5965-x

corresponding author: Elena Papaleo, elenap@cancer.dk

This repository contains curated RNASEQ data from tumor and normal solid tissue samples obtained from TCGA, and the recount2 initiative along  with their analyses. The repository was made with the intent of openly sharing both the raw input data used at the time of the analyses and the R-scripts employed to carry out the study.

The repository contains the following folders which need to be used sequentially. Each folders contain a different README file with further information:

(1) 1-download - here the download and preparation of the datasets is carried out 
(2) 2-DEA_TCGA - here the differential expression analyses and the related plots are reported
etc...

PLEASE, CITE THE PUBLICATION ABOVE IF YOU USE THE SCRIPTS FOR YOUR OWN RESEARCH

Requirements:

R version 3.3.1 or higher
Rstudio version 1.1.383 or higher        
Bioconductor version 3.6 or higher	

Other packages required:

CRAN:

VennDiagram
UpSetR
ggplot2

BIOCONDUCTOR:

TCGAbiolinks
SummarizedExperiment
edgeR
limma
sva
CEMiTools
MoonlightR
pamr
ReactomePA
ClusterProfiler

NOTES:

a) We suggest to use Rstudio to run the scripts of interest so that you can follow the analyses one line at the time and digest the results.


b) N.B. It is vital that the script TCGAbiolinks_functions.R is always run as the initial script as this scripts contains packages and costum functions needed for running the rest of the code.

c) The folder with the original data downloaded by TCGA in October 2016 for this study is too big to be shared on Github and it is available here for download:
