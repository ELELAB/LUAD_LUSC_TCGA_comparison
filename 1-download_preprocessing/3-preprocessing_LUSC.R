source("TCGAbiolinks_functions.R")
library(TCGAbiolinks)
library(SummarizedExperiment)


#---------------------------------------------------------------------------------
#                           LUSC-all
#------------------------------------------------------------------------------

SE_LUSC<-get(load("./LUSC/all/LUSC_Illumina_HiSeq_all.rda"))

#remove discordant-LUSC samples
discordant.samples <- read.table("discordant_LUSC.txt",col.names = "barcode")
SE_LUSC <- discordant_removing(discordant.samples,SE_LUSC)

#remove low tumor purity samples
list <- TCGAtumor_purity(colnames(SE_LUSC),0,0,0,0,0.6)
length(list$pure_barcodes)
length(list$filtered)
SE_LUSC <- SE_LUSC[,colnames(SE_LUSC) %in% union(list$pure_barcodes,list$filtered)]
#check the samples number
length(which(colData(SE_LUSC)$shortLetterCode =="TP"))
length(which(colData(SE_LUSC)$shortLetterCode =="NT"))
save(SE_LUSC, file = "./LUSC/all/LUSC_Illumina_HiSeq_all_tumorPurity.rda")


#Preprocessing

dataPrep <- TCGAanalyze_Preprocessing(object = SE_LUSC, cor.cut = 0.6)

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,geneInfo = geneInfo,method = "gcContent")                

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,method = "quantile",qnt.cut =  0.25)


save(dataFilt, file = "./LUSC/all/LUSC_PreprocessedData_all_tumorPurity.rda")

#----------------------------------------------------------------------------------
#                             LUSC unpaired
#------------------------------------------------------------------------------

SE_LUSC<-get(load("./LUSC/all/LUSC_Illumina_HiSeq_all.rda"))

#check the samples number
length(which(colData(SE_LUSC)$shortLetterCode =="TP"))
length(which(colData(SE_LUSC)$shortLetterCode =="NT"))

#remove discordant-LUSC samples
discordant.samples <- read.table("discordant_LUSC.txt",col.names = "barcode")
SE_LUSC <- discordant_removing(discordant.samples,SE_LUSC)

list <- TCGAtumor_purity(colnames(SE_LUSC),0,0,0,0,0.6)
length(list$pure_barcodes)
length(list$filtered)
SE_LUSC <- SE_LUSC[,colnames(SE_LUSC) %in% union(list$pure_barcodes,list$filtered)]

# This step only removes the tumor paired samples).
sampleMetaData <- colData(SE_LUSC)
paired <- TCGAquery_MatchedCoupledSampleTypes(sampleMetaData$barcode,c("NT","TP")) # find paired samples
paired_tumor_index <- lapply(paired, FUN = function(x) substr(unlist(strsplit(x, split = "-"))[[4]],1,2) == "01") # find index for TP part of paired samples
paired_tumor <- paired[unlist(paired_tumor_index)] # find TP samples of paired samples

# update dataPrep with only unpaired samples
SE_LUSC <- SE_LUSC[,!colnames(SE_LUSC) %in% paired_tumor]
save(SE_LUSC, file="./LUSC/unpaired/LUSC_Illumina_HiSeq_unpaired_tumorPurity.rda")


#Preprocessing
dataPrep <- TCGAanalyze_Preprocessing(object = SE_LUSC, cor.cut = 0.6)

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")                

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)
dim(dataFilt)

save(dataFilt, file = "./LUSC/unpaired/LUSC_PreprocessedData_unpaired_tumorPurity.rda")

#-----------------------------------------------------------------------------------
#                                   LUSC paired
#-----------------------------------------------------------------------------------

SE_LUSC<-get(load("./LUSC/paired/LUSC_Illumina_HiSeq_paired.rda"))
length(which(colData(SE_LUSC)$shortLetterCode =="TP"))
length(which(colData(SE_LUSC)$shortLetterCode =="NT"))

#remove discordant-LUSC samples
discordant.samples <- read.table("discordant_LUSC.txt",col.names = "barcode")
SE_LUSC <- discordant_removing(discordant.samples,SE_LUSC)

#remove samples with low tumor purity
list <- TCGAtumor_purity(colnames(SE_LUSC),0,0,0,0,0.6)
length(list$pure_barcodes)
length(list$filtered)
SE_LUSC <- SE_LUSC[,colnames(SE_LUSC) %in% union(list$pure_barcodes,list$filtered)]
#retrieve paired samples only
paired <- TCGAquery_MatchedCoupledSampleTypes(colnames(SE_LUSC),c("NT","TP"))
SE_LUSC <- SE_LUSC[,colnames(SE_LUSC) %in% paired]
save(SE_LUSC,file = "./LUSC/paired/LUSC_Illumina_HiSeq_paired_tumorPurity.rda")

#Preprocessing
dataPrep <- TCGAanalyze_Preprocessing(object = SE_LUSC, cor.cut = 0.6)


dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")                

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)
save(dataFilt, file = "./LUSC/paired/LUSC_PreprocessedData_paired_tumorPurity.rda")

