
#--------------------------------------------------------------------------------------
# function to create a dataframe with up-regulated genes detected by only edgeR_TCGAbiolinks
# in the rows and logFC and FDR for the three methods in the columns

# dataframe.edgeR: output dataframe obtained after DEA using edgeR
# dataframe.limma: the same but after limma
# dataframe.edgeRTCGA: after edgeR_TCGAbiolinks
# genes: list of genes that you want to get logFC and FDR
# nameFile: name of the output
#--------------------------------------------------------------------------------------

comparison_FC_FDR_3methods <- function(dataframe.edgeR,dataframe.limma,dataframe.edgeRTCGA,genes,nameFile){
  #get logFC and FDR of edgeR 
  columns.edgeR <- subset(dataframe.edgeR,rownames(dataframe.edgeR) %in% genes$V1)
  columns.edgeR <- columns.edgeR[order(rownames(columns.edgeR)),]
  columns.edgeR <- columns.edgeR[,colnames(columns.edgeR) %in% c("logFC","FDR")]
  colnames(columns.edgeR)<- c("logFC_edgeR","FDR_edgeR")
  #get logFC and FDR of limma 
  columns.limma <- subset(dataframe.limma,rownames(dataframe.limma) %in% genes$V1)
  columns.limma <- columns.limma[order(rownames(columns.limma)),]
  columns.limma <- columns.limma[,colnames(columns.limma) %in% c("logFC","adj.P.Val")]
  colnames(columns.limma)<- c("logFC_limma","FDR_limma")
  #get logFC and FDR of edgeR_TCGA 
  columns.edgeRTCGA <- subset(dataframe.edgeRTCGA,rownames(dataframe.edgeRTCGA) %in% genes$V1)
  columns.edgeRTCGA <- columns.edgeRTCGA[order(rownames(columns.edgeRTCGA)),]
  columns.edgeRTCGA <- columns.edgeRTCGA[,colnames(columns.edgeRTCGA) %in% c("logFC","FDR")]
  colnames(columns.edgeRTCGA)<- c("logFC_edgeRTCGA","FDR_edgeRTCGA")
  
  #final dataframe
  final.dataframe <- cbind(columns.edgeRTCGA,columns.edgeR,columns.limma)
  write.csv(final.dataframe,paste0(nameFile,".csv"),quote=FALSE,row.names = TRUE)
  return(final.dataframe)
}

#--------------------------------------------------------------------------------
# up-regulated genes detected only by edgeR_TCGA and not by the other methods
# LUAD paired
#-------------------------------------------------------------------------------
dataframe.edgeR <- read.csv("./LUAD/paired/edgeR_LUAD_paired_tumorPurity.csv",row.names = 1)
dataframe.limma <- read.csv("./LUAD/paired/limma_LUAD_paired_tumorPurity.csv",row.names = 1)
dataframe.edgeRTCGA<- read.csv("./LUAD/paired/edgeR_TCGA_LUAD_paired_tumorPurity.csv",row.names = 1)
# get the up-regulated genes detected by edgeR_TCGAbiolinks in paired dataset
genes <- read.table("paired/no_overlap_up_paired_TCGAedgeR_LUAD.txt")
nameFile <- "comparison_FC_FDR_up_paired_LUAD"
output<-comparison_FC_FDR_3methods(dataframe.edgeR,dataframe.limma,dataframe.edgeRTCGA,genes,nameFile)

#---------------------------------------------------------------------------------
# up-regulated genes detected by only edgeR_TCGA and not by the other methods
# LUSC paired
#--------------------------------------------------------------------------------
dataframe.edgeR<- read.csv("./LUSC/paired/edgeR_LUSC_paired_tumorPurity.csv",row.names = 1)
dataframe.limma <- read.csv("./LUSC/paired/limma_LUSC_paired_tumorPurity.csv",row.names = 1)
dataframe.edgeRTCGA <- read.csv("./LUSC/paired/edgeR_TCGA_LUSC_paired_tumorPurity.csv",row.names = 1)
genes <- read.table("paired/no_overlap_up_paired_TCGAedgeR_LUSC.txt")
nameFile <- "comparison_FC_FDR_up_paired_LUSC"
output<-comparison_FC_FDR_3methods(dataframe.edgeR,dataframe.limma,dataframe.edgeRTCGA,genes,nameFile)


