source("TCGAbiolinks_functions.R")

#library(plyr)
library(SummarizedExperiment)

limma_paired <- function(my_IDs,dataframe,limma_name,up_name,down_name){
  
  condition <- as.factor(my_IDs$condition)
  patientID <- as.factor(my_IDs$participant)
  
  #design matrix
  design.matrix <- model.matrix(~0+condition+patientID)
  colnames(design.matrix)[c(1,2)] <- c("cancer","normal")
  
  #voom transformation
  dataframe <- voom(dataframe,design.matrix,plot=TRUE)
  
  # Making group contrasts 
  N_C_cotr <- makeContrasts("cancer-normal", levels= design.matrix)
  
  # Filter for significance - set log fold change and fdr cutoff
  N_C_L <- DE_limma(N_C_cotr, dataframe, design.matrix, 1, 0.01)
  
  # differentially expressed genes file
  write.csv(N_C_L, limma_name, quote = FALSE)
  
  #number of up-regulated genes cancer vs normal
  print(paste0("the analysis detected ",length(which(N_C_L$direction == "up"))," up-regulated genes"))
  #number of down-regulated genes cancer vs normal
  print(paste0("the analysis detected ",length(which(N_C_L$direction== "down")), " down-regulated genes"))
  
  # up and down regulated genes
  up <- data.frame(rownames(N_C_L[N_C_L$direction == "up", ]))
  down <- data.frame(rownames(N_C_L[N_C_L$direction == "down", ]))
  
  write.table(up, up_name, sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  write.table(down, down_name, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
}

#--------------------------------------------------------------------------
# LUAD paired after removing low tumor purity samples
#-------------------------------------------------------------------------

dataframe_LUAD <- get(load("LUAD/paired/LUAD_PreprocessedData_paired_TumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUAD)
#check the tumor samples number
length(which(my_IDs$condition=="cancer"))
#check the normal samples number
length(which(my_IDs$condition=="normal"))
limma_name <- "limma_LUAD_paired_tumorPurity.csv"
up_name <- "limma/up_limma_LUAD_paired_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/paired/limma/down_limma_LUAD_paired_tumorPurity.txt"
limma_paired(my_IDs,dataframe_LUAD,limma_name,up_name,down_name)


#-------------------------------------------------------------------------------
# LUSC paired after removing low tumor purity
#------------------------------------------------------------------------------

dataframe_LUSC <- get(load("/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/paired/LUSC_PreprocessedData_paired_TumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUSC)
#check the tumor samples number
length(which(my_IDs$condition=="cancer"))
#check the normal samples number
length(which(my_IDs$condition=="normal"))
limma_name <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/paired/limma/limma_LUSC_paired_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/paired/limma/up_limma_LUSC_paired_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/paired/limma/down_limma_LUSC_paired_tumorPurity.txt"
limma_paired(my_IDs,dataframe_LUSC,limma_name,up_name,down_name)

#----------------------------------------------------------------------------------------------------------------------
# limma for "all" and "unpaired" samples where the tss information is included 
# into design matrix (for batch-correction)
#----------------------------------------------------------------------------------------------------------------------

limma_tss <- function(my_IDs,dataframe,limma_name,up_name,down_name){
  
  condition <- as.factor(my_IDs$condition)
  tss <- as.factor(my_IDs$tss)
  
  #design matrix
  design.matrix <- model.matrix(~0+condition+tss)
  colnames(design.matrix)[c(1,2)] <- c("cancer","normal")
  
  #voom transformation
  dataframe <- voom(dataframe,design.matrix,plot=TRUE)
  
  # Making group contrasts 
  N_C_cotr <- makeContrasts("cancer-normal", levels= design.matrix)
  
  # Filter for significance - set log fold change and fdr cutoff
  N_C_L <- DE_limma(N_C_cotr, dataframe, design.matrix, 1, 0.01)
  
  # differentially expressed genes file
  write.csv(N_C_L, limma_name, quote = FALSE)
  
  #number of up-regulated genes cancer vs normal
  print(paste0("the analysis detected ",length(which(N_C_L$direction == "up"))," up-regulated genes"))
  #number of down-regulated genes cancer vs normal
  print(paste0("the analysis detected ",length(which(N_C_L$direction== "down")), " down-regulated genes"))
  
  # up and down regulated genes
  up <- data.frame(rownames(N_C_L[N_C_L$direction == "up", ]))
  down <- data.frame(rownames(N_C_L[N_C_L$direction == "down", ]))
  
  write.table(up, up_name, sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  write.table(down, down_name, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
}


#------------------------------------------------------------------------------
# LUAD all after removing low tumor purity samples
#--------------------------------------------------------------------------
dataframe_LUAD <- get(load("/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/all/LUAD_PreprocessedData_all_TumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUAD)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
limma_name <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/all/limma/limma_LUAD_all_tss_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/all/limma/up_limma_LUAD_all_tss_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/all/limma/down_limma_LUAD_all_tss_tumorPurity.txt"

limma_tss(my_IDs,dataframe_LUAD,limma_name,up_name,down_name)


#------------------------------------------------------------------------------
# LUAD unpaired after removing low tumor purity samples
#----------------------------------------------------------------------------
dataframe_LUAD <- get(load("/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/unpaired/LUAD_PreprocessedData_unpaired_TumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUAD)
#check the tumor samples number
length(which(my_IDs$condition=="cancer"))
#check the normal samples number
length(which(my_IDs$condition=="normal"))
limma_name <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/unpaired/limma/limma_LUAD_unpaired_tss_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/unpaired/limma/up_limma_LUAD_unpaired_tss_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/unpaired/limma/down_limma_LUAD_unpaired_tss_tumorPurity.txt"
limma_tss(my_IDs,dataframe_LUAD,limma_name,up_name,down_name)


#-------------------------------------------------------------------------------
# LUSC all after removing low tumor purity samples
#-------------------------------------------------------------------------------

dataframe_LUSC <- get(load("/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/all/LUSC_PreprocessedData_all_TumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUSC)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
limma_name <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/all/limma/limma_LUSC_all_tss_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/all/limma/up_limma_LUSC_all_tss_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/all/limma/down_limma_LUSC_all_tss_tumorPurity.txt"

limma_tss(my_IDs,dataframe_LUSC,limma_name,up_name,down_name)

#---------------------------------------------------------------------------------
# LUSC unpaired after removing low tumor purity samples
#-------------------------------------------------------------------------------

dataframe_LUSC <- get(load("/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/unpaired/LUSC_PreprocessedData_unpaired_TumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUSC)
#check the tumor samples number
length(which(my_IDs$condition=="cancer"))
#check the normal samples number
length(which(my_IDs$condition=="normal"))
limma_name <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/unpaired/limma/limma_LUSC_unpaired_tss_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/unpaired/limma/up_limma_LUSC_unpaired_tss_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/unpaired/limma/down_limma_LUSC_unpaired_tss_tumorPurity.txt"
limma_tss(my_IDs,dataframe_LUSC,limma_name,up_name,down_name)
