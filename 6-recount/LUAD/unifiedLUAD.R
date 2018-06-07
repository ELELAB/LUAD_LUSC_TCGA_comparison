####Recount2 data preparation LUAD/LUNG DATASETS###############
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
luad.recount<-TCGAquery_recount2(project = "TCGA", tissue = "lung")
luad.recount.gtex<-TCGAquery_recount2(project="GTEX", tissue="lung")

query.luad<- GDCquery(project = "TCGA-LUAD",
                             data.category = "Transcriptome Profiling",
                             data.type = "Gene Expression Quantification", 
                             workflow.type = "HTSeq - Counts")


samplesDown.luad <- getResults(query.luad,cols=c("cases"))


dataSmTP.luad <- TCGAquery_SampleTypes(barcode = samplesDown.luad,
                                             typesample = "TP")

dataSmNT.luad <- TCGAquery_SampleTypes(barcode = samplesDown.luad,
                                             typesample = "NT")

query.luad2 <- GDCquery(project = "TCGA-LUAD",
                              data.category = "Transcriptome Profiling",
                              data.type = "Gene Expression Quantification", 
                              workflow.type = "HTSeq - Counts",
                              barcode = c(dataSmTP.luad, dataSmNT.luad))


GDCdownload(query=query.luad2)

dataPrep1.luad <- GDCprepare(query = query.luad2, 
                                   save = TRUE )



#####getting samples with more than 60% tumor purity and removing discordant samples
discordant.luad<-read.table("discordant_samples.txt",stringsAsFactors=FALSE)$V1
tcga.barcodes<-c(dataSmTP.luad, dataSmNT.luad)
tcga.barcodes<-setdiff(tcga.barcodes, discordant.luad)
purityinfo.R.luad<-TCGAtumor_purity(tcga.barcodes, 0, 0, 0, 0, 0.6)
dataSmTP.luad.pure.R<-purityinfo.R.luad$pure
tcga.pure_barcodes<-c(dataSmTP.luad.pure.R, dataSmNT.luad)

uuid.tcga<-TCGAtranslateID(TCGAbarcode(tcga.pure_barcodes), type="file_id")
uuid.tcga.normal<-TCGAtranslateID(TCGAbarcode(dataSmNT.luad), type="file_id")
uuid.tcga.cancer<-TCGAtranslateID(TCGAbarcode(dataSmTP.luad.pure.R), type="file_id")



#####Count tables as queried#####
eset.gtex<-assays(luad.recount.gtex$GTEX_lung)$counts
eset.tcga<-assays(luad.recount$TCGA_lung)$counts

####replacing SRR ids for colnames with GTEX barcodes 
colnames(eset.gtex)<- colData(luad.recount.gtex$GTEX_lung)$sampid

tcga.metadata<-data.frame(uid=rownames(colData(luad.recount$TCGA_lung)), 
                     gdc=colData(luad.recount$TCGA_lung)$gdc_cases.case_id,
                     cancer=colData(luad.recount$TCGA_lung)$gdc_cases.project.name,
                     tss=colData(luad.recount$TCGA_lung)$gdc_cases.tissue_source_site.code,
                     tss_name= colData(luad.recount$TCGA_lung)$gdc_cases.tissue_source_site.name)

gtex.metadata<-data.frame(id=rownames(colData(luad.recount.gtex$GTEX_lung)),
                          sampid = colData(luad.recount.gtex$GTEX_lung)$sampid,
                          tss=colData(luad.recount.gtex$GTEX_lung)$smcenter)

rowsEset.tcga.normal<-gdccases[which(gdccases$gdc %in% uuid.tcga.normal & gdccases$cancer =="Lung Adenoarcinoma"),]$uid
rowsEset.tcga.cancer<-gdccases[which(gdccases$gdc %in% uuid.tcga.cancer & gdccases$cancer =="Lung Adenocarcinoma"),]$uid

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

####Normalization- GC content######
###RECOUNT2####

eset.all.N<-TCGAanalyze_Normalization(tabDF = eset.all,
                                               geneInfo = geneInfoHT,
                                               method = "gcContent")

eset.all.NF<- TCGAanalyze_Filtering(tabDF = eset.all.N,
                                    method = "quantile", 
                                    qnt.cut =  0.25)




