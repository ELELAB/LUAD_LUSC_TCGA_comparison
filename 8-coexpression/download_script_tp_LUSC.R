library(TCGAbiolinks)
library(SummarizedExperiment)

###########################
# download all samples 
###########################

query_legacy <- GDCquery(project = "TCGA-LUSC",
                         data.category = "Gene expression",
                         data.type = "Gene expression quantification",
                         platform = "Illumina HiSeq",
                         file.type = "results",
                         legacy = TRUE,
                         experimental.strategy = "RNA-Seq",
                         sample.type = c("Primary solid Tumor"))
#GDCdownload(query_legacy)
data_legacy <- GDCprepare(query_legacy,save = TRUE, save.filename = "elena/LUSC_Illumina_HiSeq_TP.rda")
dim(data_legacy)
length(which(colData(data_legacy)$shortLetterCode =="TP"))

