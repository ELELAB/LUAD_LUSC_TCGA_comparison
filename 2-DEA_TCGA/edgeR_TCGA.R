##the user needs to run the script in a folder where already the LUAD/all LUAD/paired LUAD/unpaired LUSC/all LUSC/paired
#and LUSC/unpaired sub-folders have been created - see also README file
#the user also will need to work with the same structure of directories used in this repository to be able to run the scripts as they are
#setwd("~/2-DEA_TCGA")

source("TCGAbiolinks_functions.R")

library(TCGAbiolinks)

#-------------------------------------------------------------------------------
#                     edgeR-TCGAbiolinks
#--------------------------------------------------------------------------------


edgeR_TCGA <- function(samplesNT,samplesTP,dataFilt,edgeR_TCGA_name,up_name,down_name){

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")  

# add information about the genes modulation
my.tags <- dataDEGs
index.up <- which(my.tags$logFC >= 1 & my.tags$FDR < 0.01)
index.down <- which(my.tags$logFC <= -1 & my.tags$FDR < 0.01)
direction <- c()
direction[index.up] <- "up"
direction[index.down] <- "down"
direction[!(1:nrow(my.tags) %in% union(index.up,index.down))] <- "no DE"
my.tags <- cbind(my.tags,direction)
write.csv(my.tags,edgeR_TCGA_name, quote = FALSE)

#number of up-regulated genes cancer vs normal
print(paste0("the analysis has detected ",length(which(my.tags$direction == "up"))," up-regulated genes"))
#number of down-regulated genes cancer vs normal
print(paste0("the analysis has detected ",length(which(my.tags$direction== "down")), " down-regulated genes"))

# up and down regulated genes
up <- data.frame(rownames(my.tags[my.tags$direction == "up", ]))
down <- data.frame(rownames(my.tags[my.tags$direction == "down", ]))
write.table(up, up_name, sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
write.table(down,down_name, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}


#----------------------------------------------------------------------------------
# LUSC paired
#--------------------------------------------------------------
dataframe_LUSC <- get(load("../1-download_preprocessing/LUSC/paired/LUSC_PreprocessedData_paired_tumorPurity.rda"))
# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("TP"))
edgeR_TCGA_name <- "./LUSC/paired/edgeR_TCGA_LUSC_paired_tumorPurity.csv"
up_name <- "./LUSC/paired/up_edgeR_TCGA_LUSC_paired_tumorPurity.txt"
down_name <- "./LUSC/paired/down_edgeR_TCGA_LUSC_paired_tumorPurity.txt"

edgeR_TCGA(samplesNT,samplesTP,dataframe_LUSC,edgeR_TCGA_name,up_name,down_name)


#----------------------------------------------------------------------------
# LUAD paired
#--------------------------------------------------------------------------

dataframe_LUAD <- get(load("../1-download_preprocessing/LUAD/paired/LUAD_PreprocessedData_paired_tumorPurity.rda"))

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("TP"))
edgeR_TCGA_name <- "./LUAD/paired/edgeR_TCGA_LUAD_paired_tumorPurity.csv"
up_name <- "./LUAD/up_edgeR_TCGA_LUAD_paired_tumorPurity.txt"
down_name <- "./LUAD/paired/down_edgeR_TCGA_LUAD_paired_tumorPurity.txt"
edgeR_TCGA(samplesNT,samplesTP,dataframe_LUAD,edgeR_TCGA_name,up_name,down_name)


#--------------------------------------------------------------------------------
# LUAD all
#--------------------------------------------------------------------------------
dataframe_LUAD <- get(load("../1-download_preprocessing/LUAD/all/LUAD_PreprocessedData_all_tumorPurity.rda"))
# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("TP"))
edgeR_TCGA_name <- "./LUAD/all/edgeR_TCGA_LUAD_all_tumorPurity.csv"
up_name <- "./LUAD/all/up_edgeR_TCGA_LUAD_all_tumorPurity.txt"
down_name <- "./LUAD/all/down_edgeR_TCGA_LUAD_all_tumorPurity.txt"
edgeR_TCGA(samplesNT,samplesTP,dataframe_LUAD,edgeR_TCGA_name,up_name,down_name)


#---------------------------------------------------------------------------------
# LUAD unpaired
#---------------------------------------------------------------------------------

dataframe_LUAD <- get(load("../1-download_preprocessing/LUAD/unpaired/LUAD_PreprocessedData_unpaired_tumorPurity.rda"))
# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("TP"))
edgeR_TCGA_name <- "./LUAD/unpaired/edgeR_TCGA_LUAD_unpaired_tumorPurity.csv"
up_name <- "./LUAD/unpaired/up_edgeR_TCGA_LUAD_unpaired_tumorPurity.txt"
down_name <- "./LUAD/unpaired/down_edgeR_TCGA_LUAD_unpaired_tumorPurity.txt"
edgeR_TCGA(samplesNT,samplesTP,dataframe_LUAD,edgeR_TCGA_name,up_name,down_name)

#-------------------------------------------------------------------------------
# LUSC all
#------------------------------------------------------------------------

dataframe_LUSC <- get(load("../1-download_preprocessing/LUSC/all/LUSC_PreprocessedData_all_tumorPurity.rda"))
# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("TP"))
edgeR_TCGA_name <- "./LUSC/all/edgeR_TCGA_LUSC_all_tumorPurity.csv"
up_name <- "./LUSC/all/up_edgeR_TCGA_LUSC_all_tumorPurity.txt"
down_name <- "./LUSC/all/down_edgeR_TCGA_LUSC_all_tumorPurity.txt"

edgeR_TCGA(samplesNT,samplesTP,dataframe_LUSC,edgeR_TCGA_name,up_name,down_name)

#----------------------------------------------------------------
# LUSC unpaired
#----------------------------------------------------------------

dataframe_LUSC <- get(load("../1-download_preprocessing/LUSC/unpaired/LUSC_PreprocessedData_unpaired_tumorPurity.rda"))

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("TP"))
edgeR_TCGA_name <- "./LUSC/unpaired/edgeR_TCGA_LUSC_unpaired_tumorPurity.csv"
up_name <- "./LUSC/unpaired/up_edgeR_TCGA_LUSC_unpaired_tumorPurity.txt"
down_name <- "./LUSC/unpaired/down_edgeR_TCGA_LUSC_unpaired_tumorPurity.txt"

edgeR_TCGA(samplesNT,samplesTP,dataframe_LUSC,edgeR_TCGA_name,up_name,down_name)



