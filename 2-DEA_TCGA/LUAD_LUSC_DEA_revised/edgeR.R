
source("/data/user/marta/pipeline/DE/limma/TCGAbiolinks_functions.R")

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

dataframe_LUSC <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/paired/LUSC_PreprocessedData_woBatch-voom_paired.rda"))
my_IDs <- get_IDs(dataframe_LUSC)
#check the number of normal and tumor samples after tumor replicates removing
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/paired/edgeR/edgeR_LUSC_paired.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/paired/edgeR/up_edgeR_LUSC_paired.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/paired/edgeR/down_edgeR_LUSC_paired.txt"

edgeR_paired(dataframe_LUSC,edgeR_name,up_name,down_name)

#--------------------------------------------------------------------------------
# LUSC paired after removing low tumor purity samples
#-------------------------------------------------------------------------------

dataframe_LUSC <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/paired/LUSC_PreprocessedData_paired_TumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUSC)
#check the number of normal and tumor samples after tumor replicates removing
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/paired/edgeR/edgeR_LUSC_paired_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/paired/edgeR/up_edgeR_LUSC_paired_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/paired/edgeR/down_edgeR_LUSC_paired_tumorPurity.txt"

edgeR_paired(dataframe_LUSC,edgeR_name,up_name,down_name)

#-------------------------------------------------------------------------------
# LUAD paired after removing low tumor purity samples
#-------------------------------------------------------------------------------

dataframe_LUAD <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/paired/LUAD_PreprocessedData_paired_TumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUAD)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/paired/edgeR/edgeR_LUAD_paired_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/paired/edgeR/up_edgeR_LUAD_paired_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/paired/edgeR/down_edgeR_LUAD_paired_tumorPurity.txt"

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
dataframe_LUAD <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/all/LUAD_PreprocessedData_all.rda"))
my_IDs <- get_IDs(dataframe_LUAD)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/all/edgeR/edgeR_LUAD_all.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/all/edgeR/up_edgeR_LUAD_all.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/all/edgeR/down_edgeR_LUAD_all.txt"

edgeR_tss(dataframe_LUAD,edgeR_name,up_name,down_name)

#------------------------------------------------------------------------------
# LUAD all after removing low tumor purity samples
#--------------------------------------------------------------------------
dataframe_LUAD <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/all/LUAD_PreprocessedData_all_TumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUAD)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/all/edgeR/edgeR_LUAD_all_tss_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/all/edgeR/up_edgeR_LUAD_all_tss_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/all/edgeR/down_edgeR_LUAD_all_tss_tumorPurity.txt"

edgeR_tss(dataframe_LUAD,edgeR_name,up_name,down_name)

#------------------------------------------------------------------------------
# LUAD unpaired
#--------------------------------------------------------------------------
dataframe_LUAD <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/unpaired/LUAD_PreprocessedData_unpaired.rda"))
my_IDs <- get_IDs(dataframe_LUAD)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/unpaired/edgeR/edgeR_LUAD_unpaired.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/unpaired/edgeR/up_edgeR_LUAD_unpaired.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/unpaired/edgeR/down_edgeR_LUAD_unpaired.txt"

edgeR_tss(dataframe_LUAD,edgeR_name,up_name,down_name)

#------------------------------------------------------------------------------
# LUAD unpaired after removing low tumor purity samples
#--------------------------------------------------------------------------
dataframe_LUAD <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/unpaired/LUAD_PreprocessedData_unpaired_TumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUAD)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/unpaired/edgeR/edgeR_LUAD_unpaired_tss_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/unpaired/edgeR/up_edgeR_LUAD_unpaired_tss_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/unpaired/edgeR/down_edgeR_LUAD_unpaired_tss_tumorPurity.txt"

edgeR_tss(dataframe_LUAD,edgeR_name,up_name,down_name)


#-------------------------------------------------------------------------------
#LUSC all
#----------------------------------------------------------------------------
dataframe_LUSC <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/all/LUSC_PreprocessedData_all.rda"))
my_IDs <- get_IDs(dataframe_LUSC)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/all/edgeR/edgeR_LUSC_all.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/all/edgeR/up_edgeR_LUSC_all.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/all/edgeR/down_edgeR_LUSC_all.txt"

edgeR_tss(dataframe_LUSC,edgeR_name,up_name,down_name)

#-------------------------------------------------------------------------------
# LUSC all after removing low tumor purity samples
#-------------------------------------------------------------------------------
dataframe_LUSC <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/all/LUSC_PreprocessedData_all_TumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUSC)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/all/edgeR/edgeR_LUSC_all_tss_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/all/edgeR/up_edgeR_LUSC_all_tss_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/all/edgeR/down_edgeR_LUSC_all_tss_tumorPurity.txt"

edgeR_tss(dataframe_LUSC,edgeR_name,up_name,down_name)



#------------------------------------------------------------------------------
# LUSC unpaired
#---------------------------------------------------------------------------

dataframe_LUSC <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/unpaired/LUSC_PreprocessedData_unpaired.rda"))
my_IDs <- get_IDs(dataframe_LUSC)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/unpaired/edgeR/edgeR_LUSC_unpaired.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/unpaired/edgeR/up_edgeR_LUSC_unpaired.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/unpaired/edgeR/down_edgeR_LUSC_unpaired.txt"

edgeR_tss(dataframe_LUSC,edgeR_name,up_name,down_name)

#------------------------------------------------------------------------------
# LUSC unpaired after removing low tumor purity
#---------------------------------------------------------------------------

dataframe_LUSC <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/unpaired/LUSC_PreprocessedData_unpaired_TumorPurity.rda"))
my_IDs <- get_IDs(dataframe_LUSC)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
edgeR_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/unpaired/edgeR/edgeR_LUSC_unpaired_tss_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/unpaired/edgeR/up_edgeR_LUSC_unpaired_tss_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/unpaired/edgeR/down_edgeR_LUSC_unpaired_tss_tumorPurity.txt"

edgeR_tss(dataframe_LUSC,edgeR_name,up_name,down_name)


#----------------------------------------------------------------------------------
# edgeR correcting for plates
#---------------------------------------------------------------------------------

edgeR_plate <- function(dataframe,edgeR_name,up_name,down_name){
  
  my_IDs <- get_IDs(dataframe)
  # edgeR object
  y <- DGEList(counts=dataframe)
  
  # model matrix info
  y$samples$condition <- as.factor(my_IDs$condition)
  y$samples$plate <- as.factor(my_IDs$plate)
  
  # relevel data
  y$samples$condition = relevel(y$samples$condition, ref="normal")
  
  # design matrix
  design.mat <- model.matrix(~condition+plate, data=y$samples)
  
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

#------------------------------------------------------ 
#LUSC all without low tumor purity samples
#-------------------------------------------------------

dataframe_LUSC <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_plate_tumorPurity/LUSC_PreprocessedData_all_noTumorPurity.rda"))
dim(dataframe_LUSC)
edgeR_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_plate_tumorPurity/all/edgeR/edgeR_LUSC_all_plate_TumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_plate_tumorPurity/all/edgeR/up_edgeR_LUSC_all_plate_TumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_plate_tumorPurity/all/edgeR/down_edgeR_LUSC_all_plate_TumorPurity.txt"

edgeR_plate(dataframe_LUSC,edgeR_name,up_name,down_name)
