library(TCGAbiolinks)
library(survival)
library(plyr)
library(ggfortify)
library(ggplot2)
library(survminer)
source("functions_survival.R")
source("TCGAbiolinks_functions.R")



############################################################################
# LUAD

dataFilt <- get(load("../1-download_preprocessing/LUAD/all/LUAD_PreprocessedData_all_tumorPurity.rda"))
my_IDs <- get_IDs(dataFilt)
#remove tumor duplicates calculating the mean of them
dataFilt <- mean.duplicated.tumor(dataFilt,my_IDs)
gene_list <- c("MUC5B","HABP2","MUC21","KCNK5","ICA1","ITGA6","CSTA","P2RY1",
               "ANXA8","FZD7","CHST7","ACOX2","ALDOC","AQP5","ARSE",
               "FABP5","SLC2A9","NRCAM","AGR2", "SPDEF")

list <- c()
for(i in 1:length(gene_list)){
  
  clinical <- read.csv("./LUAD/TCGA-LUAD_clinical.csv")
  list <- get_signif_list(dataFilt,clinical,gene_list[i],0.25, 0.75,"cox",list)
}

num <- length(list) #11
final_table <- c()
for(i in 1:length(list)){
  
  clinical <- read.csv("./LUAD/TCGA-LUAD_clinical.csv")
  final_table <- cox_regression(clinical,dataFilt,list[i],0.25, 0.75,"cox",final_table,num)
}
write.csv(final_table,"cox_table_LUAD_final_table.csv")



#############################################################################
#LUSC

dataFilt <- get(load("../1-download_preprocessing/LUSC/all/LUSC_PreprocessedData_all_tumorPurity.rda"))
gene_list <- c("MUC5B","HABP2","MUC21","KCNK5","ICA1","ITGA6","CSTA","P2RY1",
               "ANXA8","FZD7","CHST7","ACOX2","ALDOC","AQP5","ARSE",
               "FABP5","SLC2A9","NRCAM","AGR2", "SPDEF")
list <- c()
for(i in 1:length(gene_list)){
  
  clinical <- read.csv("./LUSC/TCGA-LUSC_clinical.csv")
  list <- get_signif_list(dataFilt,clinical,gene_list[i],0.25, 0.75,"cox",list)
}

num <- length(list)#13
final_table <- c()

for(i in 1:length(list)){
  
  clinical <- read.csv("./LUSC/TCGA-LUSC_clinical.csv")
  final_table <- cox_regression(clinical,dataFilt,list[i],0.25,0.75,"cox",final_table,num)
  
}

write.csv(final_table,"cox_table_LUSC_final_table.csv")
