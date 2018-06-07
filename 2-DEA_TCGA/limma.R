source("TCGAbiolinks_functions.R")

#library(plyr)
library(SummarizedExperiment)

limma_paired <- function(my_IDs,dataframe,limma_name,up_name,down_name){
  
  condition <- as.factor(my_IDs$condition)
  patientID <- as.factor(my_IDs$participant)
  
  #design matrix
  design.matrix1 <- model.matrix(~condition)
  design.matrix <- model.matrix(~0+condition+patientID)
  colnames(design.matrix)[c(1,2)] <- c("cancer","normal")
  
  #voom transformation
  dataframe <- voom(dataframe,design.matrix1,plot=TRUE)
  dataframe <- dataframe$E
  
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
# LUAD paired
#-------------------------------------------------------------------------

dataframe_LUAD <- get(load("./1-download_preprocessing/LUAD/paired/LUAD_PreprocessedData_paired_tumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUAD)
#check the tumor samples number
length(which(my_IDs$condition=="cancer"))
#check the normal samples number
length(which(my_IDs$condition=="normal"))
limma_name <- "./LUAD/paired/limma_LUAD_paired_tumorPurity.csv"
up_name <- "./LUAD/paired/up_limma_LUAD_paired_tumorPurity.txt"
down_name <- "./LUAD/paired/down_limma_LUAD_paired_tumorPurity.txt"
limma_paired(my_IDs,dataframe_LUAD,limma_name,up_name,down_name)

#-------------------------------------------------------------------------------
# LUSC paired
#------------------------------------------------------------------------------

dataframe_LUSC <- get(load("./1-download_preprocessing/LUSC/paired/LUSC_PreprocessedData_paired_tumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUSC)
#check the tumor samples number
length(which(my_IDs$condition=="cancer"))
#check the normal samples number
length(which(my_IDs$condition=="normal"))
limma_name <- "./LUSC/paired/limma_LUSC_paired_tumorPurity.csv"
up_name <- "./LUSC/paired/up_limma_LUSC_paired_tumorPurity.txt"
down_name <- "./LUSC/paired/down_limma_LUSC_paired_tumorPurity.txt"
limma_paired(my_IDs,dataframe_LUSC,limma_name,up_name,down_name)

#----------------------------------------------------------------------------------------------------------------------
# limma for "all" and "unpaired" samples where the tss information is included 
# into design matrix (for batch-correction)
#----------------------------------------------------------------------------------------------------------------------

limma_tss <- function(my_IDs,dataframe,limma_name,up_name,down_name){
  
  condition <- as.factor(my_IDs$condition)
  tss <- as.factor(my_IDs$tss)
  
  #design matrix
  design.matrix1 <- model.matrix(~condition)
  design.matrix <- model.matrix(~0+condition+tss)
  colnames(design.matrix)[c(1,2)] <- c("cancer","normal")
  
  #voom transformation
  dataframe <- voom(dataframe,design.matrix1,plot=TRUE)
  dataframe <- dataframe$E
  
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
# LUAD all
#--------------------------------------------------------------------------
dataframe_LUAD <- get(load("./1-download_preprocessing/LUAD/all/LUAD_PreprocessedData_all_tumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUAD)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
limma_name <- "./LUAD/all/limma_LUAD_all_tss_tumorPurity.csv"
up_name <- "./LUAD/all/up_limma_LUAD_all_tss_tumorPurity.txt"
down_name <- "./LUAD/all/down_limma_LUAD_all_tss_tumorPurity.txt"

limma_tss(my_IDs,dataframe_LUAD,limma_name,up_name,down_name)


#------------------------------------------------------------------------------
# LUAD unpaired
#----------------------------------------------------------------------------
dataframe_LUAD <- get(load("./1-download_preprocessing/LUAD/unpaired/LUAD_PreprocessedData_unpaired_tumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUAD)
#check the tumor samples number
length(which(my_IDs$condition=="cancer"))
#check the normal samples number
length(which(my_IDs$condition=="normal"))
limma_name <- "./LUAD/unpaired/limma_LUAD_unpaired_tss_tumorPurity.csv"
up_name <- "./LUAD/unpaired/up_limma_LUAD_unpaired_tss_tumorPurity.txt"
down_name <- "./LUAD/unpaired/down_limma_LUAD_unpaired_tss_tumorPurity.txt"
limma_tss(my_IDs,dataframe_LUAD,limma_name,up_name,down_name)

#-------------------------------------------------------------------------------
# LUSC all
#-------------------------------------------------------------------------------

dataframe_LUSC <- get(load("./1-download_preprocessing/LUSC/all/LUSC_PreprocessedData_all_tumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUSC)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
limma_name <- "./LUSC/all/limma_LUSC_all_tss_tumorPurity.csv"
up_name <- "./LUSC/all/up_limma_LUSC_all_tss_tumorPurity.txt"
down_name <- "./LUSC/all/down_limma_LUSC_all_tss_tumorPurity.txt"

limma_tss(my_IDs,dataframe_LUSC,limma_name,up_name,down_name)

#---------------------------------------------------------------------------------
# LUSC unpaired 
#-------------------------------------------------------------------------------

dataframe_LUSC <- get(load("./1-download_preprocessing/LUSC/unpaired/LUSC_PreprocessedData_unpaired_TumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUSC)
#check the tumor samples number
length(which(my_IDs$condition=="cancer"))
#check the normal samples number
length(which(my_IDs$condition=="normal"))
limma_name <- "./LUSC/unpaired/limma_LUSC_unpaired_tss_tumorPurity.csv"
up_name <- "./LUSC/unpaired/up_limma_LUSC_unpaired_tss_tumorPurity.txt"
down_name <- "./LUSC/unpaired/down_limma_LUSC_unpaired_tss_tumorPurity.txt"
limma_tss(my_IDs,dataframe_LUSC,limma_name,up_name,down_name)
