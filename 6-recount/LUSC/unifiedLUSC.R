####Recount2 data preparation LUSC/LUNG DATASETS###############
######################################################
#devtools::install_github('BioinformaticsFMRP/TCGAbiolinks')
#setwd("~/Desktop/Thesis/Data_Thesis/")

###install if needed#####
#install_github("waldronlab/TCGAutils")

library(devtools)
library(TCGAbiolinks)
library(TCGAutils)
library(SummarizedExperiment)


###Query from Recount2 platform
lusc.recount<-TCGAquery_recount2(project = "TCGA", tissue = "lung")
lusc.recount.gtex<-TCGAquery_recount2(project="GTEX", tissue="lung")

query.lusc<- GDCquery(project = "TCGA-LUSC",
                             data.category = "Transcriptome Profiling",
                             data.type = "Gene Expression Quantification", 
                             workflow.type = "HTSeq - Counts")


samplesDown.lusc <- getResults(query.lusc,cols=c("cases"))


dataSmTP.lusc <- TCGAquery_SampleTypes(barcode = samplesDown.lusc,
                                             typesample = "TP")

dataSmNT.lusc <- TCGAquery_SampleTypes(barcode = samplesDown.lusc,
                                             typesample = "NT")

query.lusc2 <- GDCquery(project = "TCGA-LUSC",
                              data.category = "Transcriptome Profiling",
                              data.type = "Gene Expression Quantification", 
                              workflow.type = "HTSeq - Counts",
                              barcode = c(dataSmTP.lusc, dataSmNT.lusc))


GDCdownload(query=query.lusc2)

dataPrep1.lusc <- GDCprepare(query = query.lusc2, 
                                   save = TRUE )



#####getting samples with more than 60% tumor purity and removing discordant samples
discordant.lusc<-read.table("discordant_samples.txt",stringsAsFactors=FALSE)$V1
tcga.barcodes<-c(dataSmTP.lusc, dataSmNT.lusc)
tcga.barcodes<-setdiff(tcga.barcodes, discordant.lusc)
purityinfo.R.lusc<-TCGAtumor_purity(tcga.barcodes, 0, 0, 0, 0, 0.6)
dataSmTP.lusc.pure.R<-purityinfo.R.lusc$pure
tcga.pure_barcodes<-c(dataSmTP.lusc.pure.R, dataSmNT.lusc)

uuid.tcga<-TCGAtranslateID(TCGAbarcode(tcga.pure_barcodes), type="file_id")
uuid.tcga.normal<-TCGAtranslateID(TCGAbarcode(dataSmNT.lusc), type="file_id")
uuid.tcga.cancer<-TCGAtranslateID(TCGAbarcode(dataSmTP.lusc.pure.R), type="file_id")



#####Count tables as queried#####
eset.gtex<-assays(lusc.recount.gtex$GTEX_lung)$counts
eset.tcga<-assays(lusc.recount$TCGA_lung)$counts

####replacing SRR ids for colnames with GTEX barcodes 
colnames(eset.gtex)<- colData(lusc.recount.gtex$GTEX_lung)$sampid

tcga.metadata<-data.frame(uid=rownames(colData(lusc.recount$TCGA_lung)), 
                     gdc=colData(lusc.recount$TCGA_lung)$gdc_cases.case_id,
                     cancer=colData(lusc.recount$TCGA_lung)$gdc_cases.project.name,
                     tss=colData(lusc.recount$TCGA_lung)$gdc_cases.tissue_source_site.code,
                     tss_name= colData(lusc.recount$TCGA_lung)$gdc_cases.tissue_source_site.name)

gtex.metadata<-data.frame(id=rownames(colData(lusc.recount.gtex$GTEX_lung)),
                          sampid = colData(lusc.recount.gtex$GTEX_lung)$sampid,
                          tss=colData(lusc.recount.gtex$GTEX_lung)$smcenter)

rowsEset.tcga.normal<-gdccases[which(gdccases$gdc %in% uuid.tcga.normal & gdccases$cancer =="Lung Squamous Cell Carcinoma"),]$uid
rowsEset.tcga.cancer<-gdccases[which(gdccases$gdc %in% uuid.tcga.cancer & gdccases$cancer =="Lung Squamous Cell Carcinoma"),]$uid

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
head(gtex.metadata)

####Normalization- GC content######
###RECOUNT2####
###WORKS#######


eset.all.N<-TCGAanalyze_Normalization(tabDF = eset.all,
                                               geneInfo = geneInfoHT,
                                               method = "gcContent")

eset.all.NF<- TCGAanalyze_Filtering(tabDF = eset.all.N,
                                    method = "quantile", 
                                    qnt.cut =  0.25)
dim(eset.all)




