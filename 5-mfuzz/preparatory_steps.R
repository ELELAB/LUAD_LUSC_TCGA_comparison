library(sva)
library(limma)
source("TCGAbiolinks_functions.R")
#setwd("~/5-mfuzz/")
#-----------------------------------------------------------------------------------
#  matrices with DE genes (after removing low tumor purity samples) to be used for Mfuzz
#-----------------------------------------------------------------------------------

#LUAD paired
LUAD_paired <- get(load("../1-download_preprocessing/LUAD/paired/LUAD_PreprocessedData_paired_tumorPurity.rda"))
meta<-get_IDs(LUAD_paired)
mod_design <- model.matrix(~as.factor(meta$condition))
colnames(mod_design) <- c("cancer", "normal")
#voom-transformation
v <- voom(LUAD_paired, mod_design, plot=TRUE)
LUAD_paired <- v$E
#batch-correction
batch_corr <- ComBat(LUAD_paired, as.factor(meta$tss), mod_design, par.prior=TRUE,prior.plots=FALSE)
LUAD_paired <- batch_corr
save(LUAD_paired,file = "./preparatory_files/LUAD_PreprocessedData_voom_tss_paired.rda")

down_paired_LUAD <- read.table("../4-final_gene_list_DEA_TCGA/final_gene_list/new_down_paired_LUAD.txt")
up_paired_LUAD <- read.table("../4-final_gene_list_DEA_TCGA/final_gene_list/new_up_paired_LUAD.txt")
LUAD_paired <- LUAD_paired[rownames(LUAD_paired) %in% union(down_paired_LUAD$V1,up_paired_LUAD$V1),]
write.csv(LUAD_paired,"./preparatory_files/data_LUAD_paired_DE.csv", quote = FALSE)
LUAD_paired <-get(load("./preparatory_files/LUAD_PreprocessedData_voom_tss_paired.rda"))
LUAD_paired_noDE <- LUAD_paired[!rownames(LUAD_paired) %in% union(down_paired_LUAD$V1,up_paired_LUAD$V1),]
write.csv(LUAD_paired_noDE,"./preparatory_files/data_LUAD_paired_noDE.csv", quote = FALSE)


#---------------------------------------------------------------
#LUSC paired
#-----------------------------------------------------------------
LUSC_paired <- get(load("../1-download_preprocessing/LUSC/paired/LUSC_PreprocessedData_paired_TumorPurity.rda"))
meta<-get_IDs(LUSC_paired)
mod_design <- model.matrix(~as.factor(meta$condition))
colnames(mod_design) <- c("cancer", "normal")
#voom-transformation
v <- voom(LUSC_paired, mod_design, plot=TRUE)
LUSC_paired <- v$E
#batch-correction
batch_corr <- ComBat(LUSC_paired, as.factor(meta$tss), mod_design, par.prior=TRUE,prior.plots=FALSE)
LUSC_paired <- batch_corr
save(LUSC_paired,file = "./preparatory_files/LUSC_PreprocessedData_voom_tss_paired.rda")
down_paired_LUSC <- read.table("../4-final_gene_list_DEA_TCGA/final_gene_list/new_down_paired_LUSC.txt")
up_paired_LUSC <- read.table("../4-final_gene_list_DEA_TCGA/final_gene_list/new_up_paired_LUSC.txt")
LUSC_paired <- LUSC_paired[rownames(LUSC_paired) %in% union(down_paired_LUSC$V1,up_paired_LUSC$V1),]
write.csv(LUSC_paired,"./preparatory_files/data_LUSC_paired_DE.csv", quote = FALSE)

LUSC_paired <- get(load("./preparatory_files/LUSC_PreprocessedData_voom_tss_paired.rda"))
LUSC_paired_noDE <- LUSC_paired[!rownames(LUSC_paired) %in% union(down_paired_LUSC$V1,up_paired_LUSC$V1),]
write.csv(LUSC_paired_noDE,"./preparatory_files/data_LUSC_paired_noDE.csv", quote = FALSE)

#---------------------------------------------------------------------------------
#LUAD all
#------------------------------------------------------------------------------
LUAD_all <- get(load("../1-download_preprocessing/LUAD/all/LUAD_PreprocessedData_all_TumorPurity.rda"))
meta<-get_IDs(LUAD_all)
mod_design <- model.matrix(~as.factor(meta$condition))
colnames(mod_design) <- c("cancer", "normal")
#voom-transformation
v <- voom(LUAD_all, mod_design, plot=TRUE)
LUAD_all <- v$E
#batch-correction
batch_corr <- ComBat(LUAD_all, as.factor(meta$tss), mod_design, par.prior=TRUE,prior.plots=FALSE)
LUAD_all <- batch_corr
save(LUAD_all,file = "./prepatory_files/LUAD_PreprocessedData_voom_tss_all.rda")

down_all_LUAD <- read.table("../4-final_gene_list_DEA_TCGA/new_down_all_LUAD.txt")
up_all_LUAD <- read.table("../4-final_gene_list_DEA_TCGA/new_up_all_LUAD.txt")
LUAD_all <- LUAD_all[rownames(LUAD_all) %in% union(down_all_LUAD$V1,up_all_LUAD$V1),]
write.csv(LUAD_all,"./prepartory_files/data_LUAD_all_DE.csv", quote = FALSE)

LUAD_all <- get(load("./preparatory_files/LUAD_PreprocessedData_voom_tss_all.rda"))
LUAD_all_noDE <- LUAD_all[!rownames(LUAD_all) %in% union(down_all_LUAD$V1,up_all_LUAD$V1),]
write.csv(LUAD_all_noDE,"./preparatory_files/data_LUAD_all_noDE.csv", quote = FALSE)

#------------------------------------------------------------------------------
# LUSC all
#-----------------------------------------------------------------------------
LUSC_all <- get(load("../1-download_preprocessing/LUSC/all/LUSC_PreprocessedData_all_TumorPurity.rda"))
meta<-get_IDs(LUSC_all)
mod_design <- model.matrix(~as.factor(meta$condition))
colnames(mod_design) <- c("cancer", "normal")
#voom-transformation
v <- voom(LUSC_all, mod_design, plot=TRUE)
LUSC_all <- v$E
#batch-correction
batch_corr <- ComBat(LUSC_all, as.factor(meta$tss), mod_design, par.prior=TRUE,prior.plots=FALSE)
LUSC_all <- batch_corr
save(LUSC_all,file = "./prepatory_files/LUSC_PreprocessedData_voom_tss_all.rda")

down_all_LUSC <- read.table("../4-final_gene_list_DEA_TCGA/final_gene_list/new_down_all_LUSC.txt")
up_all_LUSC <- read.table("../4-final_gene_list_DEA_TCGA/final_gene_list/new_up_all_LUSC.txt")
LUSC_all <- LUSC_all[rownames(LUSC_all) %in% union(down_all_LUSC$V1,up_all_LUSC$V1),]
write.csv(LUSC_all,"./prepatory_files/data_LUSC_all_DE.csv", quote = FALSE)

LUSC_all <- get(load("./preparatory_files/LUSC_PreprocessedData_voom_tss_all.rda"))
LUSC_all_noDE <- LUSC_all[!rownames(LUSC_all) %in% union(down_all_LUSC$V1,up_all_LUSC$V1),]
write.csv(LUSC_all_noDE,"./preparatory_files/data_LUSC_all_noDE.csv", quote = FALSE)
