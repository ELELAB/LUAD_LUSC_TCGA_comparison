####Recount2 data preparation LUSC/LUNG DATASETS###############
######################################################
#devtools::install_github('BioinformaticsFMRP/TCGAbiolinks')
#setwd("~/6-recount")

###install if needed#####
devtools::install_github("waldronlab/TCGAutils")

library(devtools)
library(TCGAbiolinks)
library(TCGAutils)
library(SummarizedExperiment)


###Query from Recount2 platform
lusc.recount.tcga<-TCGAquery_recount2(project = "TCGA", tissue = "lung")
lusc.recount.gtex<-TCGAquery_recount2(project="GTEX", tissue="lung")
#the following lines are only to query the TCGA data 
query.lusc<- GDCquery(project = "TCGA-LUSC",
                             data.category = "Transcriptome Profiling",
                             data.type = "Gene Expression Quantification", 
                             workflow.type = "HTSeq - Counts")


samplesDown.lusc <- getResults(query.lusc,cols=c("cases"))


dataSmTP.lusc <- TCGAquery_SampleTypes(barcode = samplesDown.lusc,
                                             typesample = "TP")

dataSmNT.lusc <- TCGAquery_SampleTypes(barcode = samplesDown.lusc,
                                             typesample = "NT")



#####getting samples with more than 60% tumor purity and removing discordant samples
discordant.lusc<-read.table("./LUSC/discordant_samples.txt",stringsAsFactors=FALSE)$V1
tcga.barcodes<-c(dataSmTP.lusc, dataSmNT.lusc)
tcga.barcodes<-setdiff(tcga.barcodes, discordant.lusc)
purityinfo.R.lusc<-TCGAtumor_purity(tcga.barcodes, 0, 0, 0, 0, 0.6)
dataSmTP.lusc.pure.R<-purityinfo.R.lusc$pure
tcga.pure_barcodes<-c(dataSmTP.lusc.pure.R, dataSsmNT.lusc)

#TCGAtranslateID is deprecated
#uuid.tcga<-TCGAtranslateID(TCGAbarcode(tcga.pure_barcodes), type="file_id")
uuid.tcga<-barcodeToUUID(tcga.pure_barcodes)$case_id
uuid.tcga.normal<-barcodeToUUID(TCGAbarcode(dataSmNT.luad))$case_id
uuid.tcga.cancer<-barcodeToUUID(TCGAbarcode(dataSmTP.luad.pure.R))$case_id


#####Count tables as queried but transformed according to Rail-RNA workflow#####
eset.gtex<-assays(scale_counts(lusc.recount.gtex$GTEX_lung, round = TRUE))$counts
eset.tcga<-assays(scale_counts(lusc.recount.tcga$TCGA_lung, round = TRUE))$counts

####replacing SRR ids for colnames with GTEX barcodes 
colnames(eset.gtex)<- colData(lusc.recount.gtex$GTEX_lung)$sampid

tcga.metadata<-data.frame(uid=rownames(colData(lusc.recount$TCGA_lung)), 
                     gdc=colData(lusc.recount.tcga$TCGA_lung)$gdc_cases.case_id,
                     cancer=colData(lusc.recount.tcga$TCGA_lung)$gdc_cases.project.name,
                     tss=colData(lusc.recount.tcga$TCGA_lung)$gdc_cases.tissue_source_site.code,
                     tss_name= colData(lusc.recount.tcga$TCGA_lung)$gdc_cases.tissue_source_site.name)

gtex.metadata<-data.frame(id=rownames(colData(lusc.recount.gtex$GTEX_lung)),
                          sampid = colData(lusc.recount.gtex$GTEX_lung)$sampid,
                          tss=colData(lusc.recount.gtex$GTEX_lung)$smcenter)

rowsEset.tcga.normal<-tcga.metadata[which(tcga.metadata$gdc %in% uuid.tcga.normal & tcga.metadata$cancer =="Lung Squamous Cell Carcinoma"),]$uid
rowsEset.tcga.cancer<-tcga.metadata[which(tcga.metadata$gdc %in% uuid.tcga.cancer & tcga.metadata$cancer =="Lung Squamous Cell Carcinoma"),]$uid

eset.tcga.normal<-eset.tcga[,toupper(rowsEset.tcga.normal)]
eset.tcga.cancer<-eset.tcga[,toupper(rowsEset.tcga.cancer)]


###ALL DATA ready for normalization, filtering, and differential expression analysis
eset.all<-cbind(eset.tcga.normal, eset.tcga.cancer, eset.gtex)

####Removing whats after "." in ENSG (ENSMBL) ids
rownames(eset.all)<-gsub("\\..*", "",rownames(eset.all))
#colnames(eset.all)

####Saving data
save(eset.tcga.normal, file="tcga_normal_lusc.RData")
save(eset.tcga.cancer, file="tcga_cancer_lusc.RData")
save(eset.gtex, file="gtex.RData")
save(tcga.metadata, file="tcga_meta_lusc.RData")
save(gtex.metadata, file="gtex_meta_lusc.RData")

