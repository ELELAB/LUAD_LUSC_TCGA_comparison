#the user needs to run the script in a folder where already the LUAD/all LUAD/paired LUAD/unpaired LUSC/all LUSC/paired
#and LUSC/unpaired sub-folders have been created

source("TCGAbiolinks_functions.R")
library(TCGAbiolinks)
library(SummarizedExperiment)

#---------------------------------------------------
# download all samples - LUAD
#---------------------------------------------------

query_LUAD <- GDCquery(project = "TCGA-LUAD",
                         data.category = "Gene expression",
                         data.type = "Gene expression quantification",
                         platform = "Illumina HiSeq",
                         file.type = "results",
                         legacy = TRUE,
                         experimental.strategy = "RNA-Seq",
                         sample.type = c("Primary solid Tumor","Solid Tissue Normal"))
##if the results from our study needs to be fully reproduced, the user would need to download the original GDC folder
##corresponding from the data that we obtained in October 2016. Uncommenting the following line will download a more recent
##version of the data where the number of genes/samples could be different
#GDCdownload(query_LUAD)
data_LUAD <- GDCprepare(query_LUAD,save = TRUE, save.filename = "./LUAD/all/LUAD_Illumina_HiSeq_all.rda")
dim(data_LUAD)
length(which(colData(data_LUAD)$shortLetterCode =="TP"))
length(which(colData(data_LUAD)$shortLetterCode =="NT"))


#---------------------------------------------------
# download paired samples LUAD
#---------------------------------------------------

query_LUADpaired <- GDCquery(project = "TCGA-LUAD", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type = "results",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Primary solid Tumor","Solid Tissue Normal"))
#GDCdownload(query_LUADpaired)

# Which samples are primary solid tumor
dataSmTP <- TCGAquery_SampleTypes(query_LUADpaired$results[[1]]$cases,"TP") 
# which samples are solid tissue normal
dataSmNT <- TCGAquery_SampleTypes(query_LUADpaired$results[[1]]$cases,"NT")

matched <- TCGAquery_MatchedCoupledSampleTypes(c(dataSmTP, dataSmNT), c("TP","NT"))

query_LUADpaired2 <- GDCquery(project = "TCGA-LUAD", 
                        legacy = TRUE,
                        data.category = "Gene expression",
                        data.type = "Gene expression quantification",
                        platform = "Illumina HiSeq", 
                        file.type = "results",
                        experimental.strategy = "RNA-Seq",
                        sample.type = c("Primary solid Tumor","Solid Tissue Normal"),
                        barcode = matched)
#GDCdownload(query_LUADpaired2)

data_LUADpaired <- GDCprepare(query = query_LUADpaired2, save = TRUE, save.filename = "./LUAD/paired/LUAD_Illumina_HiSeq_paired.rda")

### the unpaired samples will be aggregated in the script 2-preprocessingLUAD.R since they requires a curation of the dataset all ###

#the user needs to run the script in a folder where already the LUAD/all LUAD/paired LUAD/unpaired LUSC/all LUSC/paired
#and LUSC/unpaired sub-folders have been created


#---------------------------------------------------
# download all samples - LUSC
#---------------------------------------------------

query_LUSC <- GDCquery(project = "TCGA-LUSC",
                       data.category = "Gene expression",
                       data.type = "Gene expression quantification",
                       platform = "Illumina HiSeq",
                       file.type = "results",
                       legacy = TRUE,
                       experimental.strategy = "RNA-Seq",
                       sample.type = c("Primary solid Tumor","Solid Tissue Normal"))
##if the results from our study needs to be fully reproduced, the user would need to download the original GDC folder
##corresponding from the data that we obtained in October 2016. Uncommenting the following line will download a more recent
##version of the data where the number of genes/samples could be different
#GDCdownload(query_LUSC)
data_LUSC <- GDCprepare(query_LUSC,save = TRUE, save.filename = "./LUSC/all/LUSC_Illumina_HiSeq_all.rda")
dim(data_LUSC)
length(which(colData(data_LUSC)$shortLetterCode =="TP"))
length(which(colData(data_LUSC)$shortLetterCode =="NT"))


#---------------------------------------------------
# download paired samples LUSC
#---------------------------------------------------

query_LUSCpaired <- GDCquery(project = "TCGA-LUSC", 
                             legacy = TRUE,
                             data.category = "Gene expression",
                             data.type = "Gene expression quantification",
                             platform = "Illumina HiSeq", 
                             file.type = "results",
                             experimental.strategy = "RNA-Seq",
                             sample.type = c("Primary solid Tumor","Solid Tissue Normal"))
#GDCdownload(query_LUSCpaired)

# Which samples are primary solid tumor
dataSmTP <- TCGAquery_SampleTypes(query_LUSCpaired$results[[1]]$cases,"TP") 
# which samples are solid tissue normal
dataSmNT <- TCGAquery_SampleTypes(query_LUSCpaired$results[[1]]$cases,"NT")

matched <- TCGAquery_MatchedCoupledSampleTypes(c(dataSmTP, dataSmNT), c("TP","NT"))

query_LUSCpaired2 <- GDCquery(project = "TCGA-LUSC", 
                              legacy = TRUE,
                              data.category = "Gene expression",
                              data.type = "Gene expression quantification",
                              platform = "Illumina HiSeq", 
                              file.type = "results",
                              experimental.strategy = "RNA-Seq",
                              sample.type = c("Primary solid Tumor","Solid Tissue Normal"),
                              barcode = matched)
#GDCdownload(query_LUSCpaired2)

data_LUSCpaired <- GDCprepare(query = query_LUSCpaired2, save = TRUE, save.filename = "./LUSC/paired/LUSC_Illumina_HiSeq_paired.rda")

### the unpaired samples will be aggregated in the script 2-preprocessingLUSC.R since they requires a curation of the dataset all ###
#### the 19 discordant LUSC samples will be removed in the 2-preprocessingLUSC.R steps ####