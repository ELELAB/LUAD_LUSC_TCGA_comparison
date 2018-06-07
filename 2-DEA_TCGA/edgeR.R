source("TCGAbiolinks_functions.R")

#------------------------------------------------------------------------------------

#DIFFERENTIAL EXPRESSION ANALYSIS - edgeR

# edgeR_paired function takes into account the patient information (a paired dateset is used)
# edgeR_tss function takes into account the tss information (for batch-correction)
#-------------------------------------------------------------------------------------

# dataframe: .rda file after normalization and filtering
# edgeR_name: path and file name of output
# up_name and down_name: path and file name where the up and down regulated genes
#                        are saved, respectively

edgeR_paired <- function(dataframe,edgeR_name,up_name,down_name){
my_IDs <- get_IDs(dataframe)
# edgeR object
y <- DGEList(counts=dataframe)

# model matrix info
y$samples$condition <- as.factor(my_IDs$condition)
y$samples$patientID <- as.factor(my_IDs$participant)

# relevel data
y$samples$condition = relevel(y$samples$condition, ref="normal")

# design matrix
design.mat <- model.matrix(~condition+patientID, data=y$samples)

# estimate dispersion
y <- estimateDisp(y,design.mat)

# fit overdispersed poisson model
my.fit <- glmFit(y, design.mat)

# Performing likelihood ratio test
N_C_coef <- glmLRT(my.fit, coef = 2)

# Filter for significance - set log fold change and fdr cutoff
N_C_E <- DE_edgeR(N_C_coef, y, 1, 0.01)

# differentially expressed genes file

write.csv(N_C_E, edgeR_name, quote = FALSE)

#number of up-regulated genes cancer vs normal
print(paste0("the analysis detected ",length(which(N_C_E$direction == "up"))," up-regulated genes"))
#number of down-regulated genes cancer vs normal
print(paste0("the analysis detected ",length(which(N_C_E$direction == "down"))," down-regulated genes"))

# up and down regulated genes
up <- data.frame(rownames(N_C_E[N_C_E$direction == "up", ]))
down <- data.frame(rownames(N_C_E[N_C_E$direction == "down", ]))

write.table(up, up_name, sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
write.table(down,down_name, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

}


#--------------------------------------------------------------------------------
# LUSC paired
#-------------------------------------------------------------------------------

dataframe_LUSC <- get(load("./1-download_preprocessing/LUSC/paired/LUSC_PreprocessedData_paired_tumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUSC)
#check the number of normal and tumor samples after tumor replicates removing
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "./LUSC/paired/edgeR_LUSC_paired_tumorPurity.csv"
up_name <- "./LUSC/paired/up_edgeR_LUSC_paired_tumorPurity.txt"
down_name <- "./LUSC/paired/down_edgeR_LUSC_paired_tumorPurity.txt"

edgeR_paired(dataframe_LUSC,edgeR_name,up_name,down_name)

#-------------------------------------------------------------------------------
# LUAD paired
#-------------------------------------------------------------------------------

dataframe_LUAD <- get(load("./1-download_preprocessing/LUAD/paired/LUAD_PreprocessedData_paired_tumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUAD)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "./LUAD/paired/edgeR_LUAD_paired_tumorPurity.csv"
up_name <- "./LUAD/paired/up_edgeR_LUAD_paired_tumorPurity.txt"
down_name <- "./LUAD/paired/down_edgeR_LUAD_paired_tumorPurity.txt"

edgeR_paired(dataframe_LUAD,edgeR_name,up_name,down_name)

#----------------------------------------------------------------------------------------------------------------------
# edgeR for "all" and "unpaired" samples where the tss information is included into design matrix (for batch-correction)
#----------------------------------------------------------------------------------------------------------------------

edgeR_tss <- function(dataframe,edgeR_name,up_name,down_name){
  
  my_IDs <- get_IDs(dataframe)
  # edgeR object
  y <- DGEList(counts=dataframe)
  
  # model matrix info
  y$samples$condition <- as.factor(my_IDs$condition)
  y$samples$tss <- as.factor(my_IDs$tss)
  
  # relevel data
  y$samples$condition = relevel(y$samples$condition, ref="normal")
  
  # design matrix
  design.mat <- model.matrix(~condition+tss, data=y$samples)
  
  # estimate dispersion
  y <- estimateDisp(y,design.mat)
  
  # fit overdispersed poisson model
  my.fit <- glmFit(y, design.mat)
  
  # Performing likelihood ratio test
  N_C_coef <- glmLRT(my.fit, coef = 2)
  
  # Filter for significance - set log fold change and fdr cutoff
  N_C_E <- DE_edgeR(N_C_coef, y, 1, 0.01)
  
  # differentially expressed genes file
  
  write.csv(N_C_E, edgeR_name, quote = FALSE)
  
  #number of up-regulated genes cancer vs normal
  print(paste0("the analysis detected ",length(which(N_C_E$direction == "up"))," up-regulated genes"))
  #number of down-regulated genes cancer vs normal
  print(paste0("the analysis detected ",length(which(N_C_E$direction == "down"))," down-regulated genes"))
  
  # up and down regulated genes
  up <- data.frame(rownames(N_C_E[N_C_E$direction == "up", ]))
  down <- data.frame(rownames(N_C_E[N_C_E$direction == "down", ]))
  
  write.table(up, up_name, sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  write.table(down,down_name, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
}

#------------------------------------------------------------------------------
# LUAD all 
#--------------------------------------------------------------------------

dataframe_LUAD <- get(load("./1-download_preprocessing/LUAD/all/LUAD_PreprocessedData_all_tumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUAD)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "./LUAD/all/edgeR_LUAD_all_tss_tumorPurity.csv"
up_name <- "./LUAD/all/up_edgeR_LUAD_all_tss_tumorPurity.txt"
down_name <- "./LUAD/all/down_edgeR_LUAD_all_tss_tumorPurity.txt"

edgeR_tss(dataframe_LUAD,edgeR_name,up_name,down_name)


#------------------------------------------------------------------------------
# LUAD unpaired
#--------------------------------------------------------------------------
dataframe_LUAD <- get(load("./1-download_preprocessing/LUAD/unpaired/LUAD_PreprocessedData_unpaired_tumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUAD)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "./LUAD/unpaired/edgeR_LUAD_unpaired_tss_tumorPurity.csv"
up_name <- "./LUAD/unpaired/up_edgeR_LUAD_unpaired_tss_tumorPurity.txt"
down_name <- "./LUAD/unpaired/down_edgeR_LUAD_unpaired_tss_tumorPurity.txt"

edgeR_tss(dataframe_LUAD,edgeR_name,up_name,down_name)


#-------------------------------------------------------------------------------
# LUSC all
#-------------------------------------------------------------------------------
dataframe_LUSC <- get(load("./1-download_preprocessing/LUSC/all/LUSC_PreprocessedData_all_tumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUSC)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "./LUSC/all/edgeR_LUSC_all_tss_tumorPurity.csv"
up_name <- "./LUSC/all/up_edgeR_LUSC_all_tss_tumorPurity.txt"
down_name <- "./LUSC/all/down_edgeR_LUSC_all_tss_tumorPurity.txt"

edgeR_tss(dataframe_LUSC,edgeR_name,up_name,down_name)


#------------------------------------------------------------------------------
# LUSC unpaired
#---------------------------------------------------------------------------

dataframe_LUSC <- get(load("./1-download_preprocessing/LUSC/unpaired/LUSC_PreprocessedData_unpaired_tumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUSC)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "./LUSC/unpaired/edgeR_LUSC_unpaired_tss_tumorPurity.csv"
up_name <- "./LUSC/unpaired/up_edgeR_LUSC_unpaired_tss_tumorPurity.txt"
down_name <- "./LUSC/unpaired/down_edgeR_LUSC_unpaired_tss_tumorPurity.txt"

edgeR_tss(dataframe_LUSC,edgeR_name,up_name,down_name)
