####Recount2 data preparation LUAD/LUNG DATASETS###############
######################################################
devtools::install_github('BioinformaticsFMRP/TCGAbiolinks')
#setwd("~/6-recount")

###install if needed#####
devtools::install_github("waldronlab/TCGAutils")


library(devtools)
library(TCGAbiolinks)
library(TCGAutils)
library(SummarizedExperiment)
library(recount)


###Query from Recount2 platform
luad.recount.tcga<-TCGAquery_recount2(project = "TCGA", tissue = "lung")
luad.recount.gtex<-TCGAquery_recount2(project="GTEX", tissue="lung")
#the following lines are only to query the TCGA data 
query.luad<- GDCquery(project = "TCGA-LUAD",
                             data.category = "Transcriptome Profiling",
                             data.type = "Gene Expression Quantification", 
                             workflow.type = "HTSeq - Counts")


samplesDown.luad <- getResults(query.luad,cols=c("cases"))


dataSmTP.luad <- TCGAquery_SampleTypes(barcode = samplesDown.luad,
                                             typesample = "TP")

dataSmNT.luad <- TCGAquery_SampleTypes(barcode = samplesDown.luad,
                                             typesample = "NT")


#####getting samples with more than 60% tumor purity and removing discordant samples

tcga.barcodes<-c(dataSmTP.luad, dataSmNT.luad)
purityinfo.R.luad<-TCGAtumor_purity(tcga.barcodes, 0, 0, 0, 0, 0.6)
dataSmTP.luad.pure.R<-purityinfo.R.luad$pure
tcga.pure_barcodes<-c(dataSmTP.luad.pure.R, dataSmNT.luad)


uuid.tcga<-barcodeToUUID(tcga.pure_barcodes)$case_id
uuid.tcga.normal<-barcodeToUUID(dataSmNT.luad)$case_id
uuid.tcga.cancer<-barcodeToUUID(dataSmTP.luad.pure.R)$case_id


#####Count tables as queried but transformed according to Rail-RNA workflow#####
eset.gtex<-assays(scale_counts(luad.recount.gtex$GTEX_lung, round = TRUE))$counts
eset.tcga<-assays(scale_counts(luad.recount.tcga$TCGA_lung, round = TRUE))$counts

####replacing SRR ids for colnames with GTEX barcodes 
colnames(eset.gtex)<- colData(luad.recount.gtex$GTEX_lung)$sampid

tcga.metadata<-data.frame(uid=rownames(colData(luad.recount.tcga$TCGA_lung)), 
                     gdc=colData(luad.recount.tcga$TCGA_lung)$gdc_cases.case_id,
                     cancer=colData(luad.recount.tcga$TCGA_lung)$gdc_cases.project.name,
                     tss=colData(luad.recount.tcga$TCGA_lung)$gdc_cases.tissue_source_site.code,
                     tss_name= colData(luad.recount.tcga$TCGA_lung)$gdc_cases.tissue_source_site.name)

gtex.metadata<-data.frame(id=rownames(colData(luad.recount.gtex$GTEX_lung)),
                          sampid = colData(luad.recount.gtex$GTEX_lung)$sampid,
                          tss=colData(luad.recount.gtex$GTEX_lung)$smcenter)

rowsEset.tcga.normal<-tcga.metadata[which(tcga.metadata$gdc %in% uuid.tcga.normal & tcga.metadata$cancer =="Lung Adenocarcinoma"),]$uid
rowsEset.tcga.cancer<-tcga.metadata[which(tcga.metadata$gdc %in% uuid.tcga.cancer & tcga.metadata$cancer =="Lung Adenocarcinoma"),]$uid

eset.tcga.normal<-eset.tcga[,toupper(rowsEset.tcga.normal)]
eset.tcga.cancer<-eset.tcga[,toupper(rowsEset.tcga.cancer)]


###ALL DATA ready for normalization, filtering, and differential expression analysis
eset.all<-cbind(eset.tcga.normal, eset.tcga.cancer, eset.gtex)

####Removing whats after "." in ENSG (ENSMBL) ids
rownames(eset.all)<-gsub("\\..*", "",rownames(eset.all))
#colnames(eset.all)

####Saving data
save(eset.tcga.normal, file="tcga_normal_luad.RData")
save(eset.tcga.cancer, file="tcga_cancer_luad.RData")
save(eset.gtex, file="gtex.RData")
save(tcga.metadata, file="tcga_meta_luad.RData")
save(gtex.metadata, file="gtex_meta_luad.RData")
head(gtex.metadata)

