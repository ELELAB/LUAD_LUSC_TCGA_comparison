#### survival analysis
library(TCGAbiolinks)
library(survival)
library(plyr)
library(ggfortify)
library(ggplot2)
library(survminer)
source("/data/user/marta/pipeline/DE/new_LUAD-LUSC/survival_analysis/functions_survival.R")
source("/data/user/marta/pipeline/DE/limma/TCGAbiolinks_functions.R")

#------------- LUSC----------------------------------------------------------------

setwd("./LUSC")
gene_list <- c("MUC5B","HABP2","MUC21","KCNK5","ICA1","ITGA6","CSTA","P2RY1",
               "ANXA8","FZD7","CHST7","RND3","ACOX2","ALDOC","AQP5","ARSE",
               "FABP5","SIPA1L2","SLC1A3","SLC2A9","NRCAM","AGR2", "SPDEF")
#get dataFilt
dataFilt <- get(load("../../LUSC/all/LUSC_PreprocessedData_all_tumorPurity.rda"))

for(i in 1:length(gene_list)){
  
  clinical <- read.csv("TCGA-LUSC_clinical.csv")
  clinical <- get_survival_table(dataFilt,clinical,gene_list[i], 0.25,0.75,"KM")
  survival_plot(clinical,gene_list[i],"LUSC",25)

}

#------------ LUAD-------------------------------------------------------------------
setwd("../LUAD")

gene_list <- c("MUC5B","HABP2","MUC21","KCNK5","ICA1","ITGA6","CSTA","P2RY1",
               "ANXA8","FZD7","CHST7","RND3","ACOX2","ALDOC","AQP5","ARSE",
               "FABP5","SIPA1L2","SLC1A3","SLC2A9","NRCAM","AGR2", "SPDEF")

dataFilt <- get(load("../../LUAD/all/LUAD_PreprocessedData_all_tumorPurity.rda"))
my_IDs <- get_IDs(dataFilt)
#remove tumor duplicates calculating the mean of them
dataFilt <- mean.duplicated.tumor(dataFilt,my_IDs)

for(i in 1:length(gene_list)){
  
  clinical <- read.csv("TCGA-LUAD_clinical.csv")
  clinical <- get_survival_table(dataFilt,clinical,gene_list[i], 0.25,0.75,"KM")
  survival_plot(clinical,gene_list[i],"LUAD",25)
}

