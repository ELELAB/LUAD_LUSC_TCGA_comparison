source("TCGAbiolinks_functions.R")
library(TCGAbiolinks)
library(SummarizedExperiment)


#---------------------------------------------------------------------------------
#                           LUAD-all
#------------------------------------------------------------------------------

SE_LUAD <- get(load("./LUAD/all/LUAD_Illumina_HiSeq_all.rda"))

#check the samples number
length(which(colData(SE_LUAD)$shortLetterCode =="TP"))
length(which(colData(SE_LUAD)$shortLetterCode =="NT"))

#remove samples with low tumor purity
list <- TCGAtumor_purity(colnames(SE_LUAD),0,0,0,0,0.6)
length(list$pure_barcodes)
length(list$filtered)
SE_LUAD <- SE_LUAD[,colnames(SE_LUAD) %in% union(list$pure_barcodes,list$filtered)]
save(SE_LUAD, file="./LUAD/all/LUAD_Illumina_HiSeq_all_tumorPurity.rda")

#Preprocessing
dataPrep <- TCGAanalyze_Preprocessing(object = SE_LUAD, cor.cut = 0.6)


dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")                

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

save(dataFilt, file = "./LUAD/all/LUAD_PreprocessedData_all_tumorPurity.rda")

#----------------------------------------------------------------------------------
#                             LUAD unpaired
#------------------------------------------------------------------------------

SE_LUAD <- get(load("./LUAD/all/LUAD_Illumina_HiSeq_all.rda"))

#check the samples number
length(which(colData(SE_LUAD)$shortLetterCode =="TP"))
length(which(colData(SE_LUAD)$shortLetterCode =="NT"))

#remove samples with low tumor purity
list <- TCGAtumor_purity(colnames(SE_LUAD),0,0,0,0,0.6)
length(list$pure_barcodes)
length(list$filtered)
SE_LUAD <- SE_LUAD[,colnames(SE_LUAD) %in% union(list$pure_barcodes,list$filtered)]

# This step only removes the tumor paired samples.
sampleMetaData <- colData(SE_LUAD)
paired <- TCGAquery_MatchedCoupledSampleTypes(sampleMetaData$barcode,c("NT","TP")) # find paired samples
paired_tumor_index <- lapply(paired, FUN = function(x) substr(unlist(strsplit(x, split = "-"))[[4]],1,2) == "01") # find index for TP part of paired samples
paired_tumor <- paired[unlist(paired_tumor_index)] # find TP samples of paired samples

# update SE_LUAD with only unpaired samples
SE_LUAD <- SE_LUAD[,!colnames(SE_LUAD) %in% paired_tumor]

save(SE_LUAD, file="./LUAD/unpaired/LUAD_Illumina_HiSeq_unpaired_tumorPurity.rda")

#Preprocessing
dataPrep <- TCGAanalyze_Preprocessing(object = SE_LUAD, cor.cut = 0.6)


dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")                

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

save(dataFilt,file="./LUAD/unpaired/LUAD_PreprocessedData_unpaired_tumorPurity.rda")


#-------------------------------------------------------------------------------------
#                                   LUAD paired
#------------------------------------------------------------------------------------

SE_LUAD<-get(load("./LUAD/paired/LUAD_Illumina_HiSeq_paired.rda"))

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
save(SE_LUAD, file = "./LUAD/paired/LUAD_Illumina_HiSeq_paired_tumorPurity.rda")

#Preprocessing
dataPrep <- TCGAanalyze_Preprocessing(object = SE_LUAD, cor.cut = 0.6)

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")                

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

save(dataFilt, file = "./LUAD/paired/LUAD_PreprocessedData_paired_tumorPurity.rda")
