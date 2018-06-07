source("/data/user/marta/pipeline/DE/limma/TCGAbiolinks_functions.R")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(limma)

#---------------------------------------------------------------------------------
#                           LUAD-all
#------------------------------------------------------------------------------

SE_LUAD<-get(load("/data/user/marta/pipeline/DE/limma/LUAD/LUAD_Illumina_HiSeq_all.rda"))

#check the samples number
length(which(colData(SE_LUAD)$shortLetterCode =="TP"))
length(which(colData(SE_LUAD)$shortLetterCode =="NT"))

#remove samples with low tumor purity
list <- TCGAtumor_purity(colnames(SE_LUAD),0,0,0,0,0.6)
length(list$pure_barcodes)
length(list$filtered)
SE_LUAD <- SE_LUAD[,colnames(SE_LUAD) %in% union(list$pure_barcodes,list$filtered)]
save(SE_LUAD, file="./LUAD_tss_tumorPurity/all/LUAD_Illumina_HiSeq_all_tumorPurity.rda")
genes_SE <- gsub("\\|.*","",rownames(SE_LUAD))

#Preprocessing
dataPrep <- TCGAanalyze_Preprocessing(object = SE_LUAD, cor.cut = 0.6)


dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")                

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)


save(dataFilt,file="./all/LUAD_PreprocessedData_all.rda")
save(dataFilt, file = "./all/LUAD_PreprocessedData_all_noTumorPurity.rda")

#----------------------------------------------------------------------------------
#                             LUAD unpaired
#------------------------------------------------------------------------------

SE_LUAD<-get(load("/data/user/marta/pipeline/DE/limma/LUAD/LUAD_Illumina_HiSeq_all.rda"))

#check the samples number
length(which(colData(SE_LUAD)$shortLetterCode =="TP"))
length(which(colData(SE_LUAD)$shortLetterCode =="NT"))

#remove samples with low tumor purity
list <- TCGAtumor_purity(colnames(SE_LUAD),0,0,0,0,0.6)
length(list$pure_barcodes)
length(list$filtered)
SE_LUAD <- SE_LUAD[,colnames(SE_LUAD) %in% union(list$pure_barcodes,list$filtered)]

# This step only removes the tumor paired samples).
sampleMetaData <- colData(SE_LUAD)
paired <- TCGAquery_MatchedCoupledSampleTypes(sampleMetaData$barcode,c("NT","TP")) # find paired samples
paired_tumor_index <- lapply(paired, FUN = function(x) substr(unlist(strsplit(x, split = "-"))[[4]],1,2) == "01") # find index for TP part of paired samples
paired_tumor <- paired[unlist(paired_tumor_index)] # find TP samples of paired samples
sampleMetaData <- sampleMetaData[!sampleMetaData$barcode %in% paired_tumor, ] # update sampleMetaData with only unpaired samples
# save sampleMetaData with only unpaired samples
write.table(subset(sampleMetaData, select = c("barcode", "shortLetterCode")),"sampleMetaData_pairedSamplesRemoved.txt", quote=F,sep = ",", row.names = F)

# update SE_LUAD with only unpaired samples
SE_LUAD <- SE_LUAD[,!colnames(SE_LUAD) %in% paired_tumor]

save(SE_LUAD, file="/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/unpaired/LUAD_Illumina_HiSeq_unpaired_tumorPurity.rda")

#Preprocessing
dataPrep <- TCGAanalyze_Preprocessing(object = SE_LUAD, cor.cut = 0.6)


dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")                

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

save(dataFilt,file="./LUAD_tss_tumorPurity/unpaired/LUAD_PreprocessedData_unpaired_TumorPurity.rda")
save(dataFilt,file="./unpaired/LUAD_PreprocessedData_unpaired.rda")


#-------------------------------------------------------------------------------------
#                                   LUAD paired
#------------------------------------------------------------------------------------

SE_LUAD<-get(load("/data/user/marta/pipeline/DE/limma/LUAD/LUAD_Illumina_HiSeq_paired.rda"))

#check the samples number
length(which(colData(SE_LUAD)$shortLetterCode =="TP"))
length(which(colData(SE_LUAD)$shortLetterCode =="NT"))

#remove samples with low tumor purity
list <- TCGAtumor_purity(colnames(SE_LUAD),0,0,0,0,0.6)
length(list$pure_barcodes)
length(list$filtered)
SE_LUAD <- SE_LUAD[,colnames(SE_LUAD) %in% union(list$pure_barcodes,list$filtered)]
#get paired samples only
paired <- TCGAquery_MatchedCoupledSampleTypes(colnames(SE_LUAD),c("NT","TP"))
SE_LUAD <- SE_LUAD[,colnames(SE_LUAD) %in% paired]
save(SE_LUAD, file = "./LUAD_tss_tumorPurity/paired/LUAD_Iluumina_HiSeq_paired_tumorPurity.rda")

#Preprocessing
dataPrep <- TCGAanalyze_Preprocessing(object = SE_LUAD, cor.cut = 0.6)


dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")                

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

save(dataFilt, file = "./LUAD_tss_tumorPurity/paired/LUAD_PreprocessedData_paired_TumorPurity.rda")
