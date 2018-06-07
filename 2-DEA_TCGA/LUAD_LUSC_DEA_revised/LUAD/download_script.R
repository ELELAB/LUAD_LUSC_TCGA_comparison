library(TCGAbiolinks)
library(SummarizedExperiment)

# download "all" samples 
query_legacy <- GDCquery(project = "TCGA-LUAD",
                         data.category = "Gene expression",
                         data.type = "Gene expression quantification",
                         platform = "Illumina HiSeq",
                         file.type = "results",
                         legacy = TRUE,
                         experimental.strategy = "RNA-Seq",
                         sample.type = c("Primary solid Tumor","Solid Tissue Normal"))

GDCdownload(query_legacy)
data_legacy <- GDCprepare(query_legacy,save = TRUE, save.filename = "LUAD_Illumina_HiSeq_all.rda")
dim(data_legacy)
length(which(colData(data_legacy)$shortLetterCode =="TP"))
length(which(colData(data_legacy)$shortLetterCode =="NT"))


#---------------------------------------------------
# download paired samples
#---------------------------------------------------

query.exp <- GDCquery(project = "TCGA-LUAD", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type = "results",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Primary solid Tumor","Solid Tissue Normal"))
GDCdownload(query.exp)

# Which samples are primary solid tumor
dataSmTP <- TCGAquery_SampleTypes(query.exp$results[[1]]$cases,"TP") 
# which samples are solid tissue normal
dataSmNT <- TCGAquery_SampleTypes(query.exp$results[[1]]$cases,"NT")

matched <- TCGAquery_MatchedCoupledSampleTypes(c(dataSmTP, dataSmNT), c("TP","NT"))

query.exp.p <- GDCquery(project = "TCGA-LUAD", 
                        legacy = TRUE,
                        data.category = "Gene expression",
                        data.type = "Gene expression quantification",
                        platform = "Illumina HiSeq", 
                        file.type = "results",
                        experimental.strategy = "RNA-Seq",
                        sample.type = c("Primary solid Tumor","Solid Tissue Normal"),
                        barcode = matched)
GDCdownload(query.exp.p)

brca.exp.p <- GDCprepare(query = query.exp.p, save = TRUE, save.filename = "LUAD_Illumina_HiSeq_paired.rda")







