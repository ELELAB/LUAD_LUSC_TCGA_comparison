#this script needs to be run after performing the DEA with the three pipelines
##the user also will need to work with the same structure of directories used in this repository to be able to run the scripts as they are
#setwd("~/2-DEA_TCGA")
#to run after the script comparison_FC_FDR.R
library(ggplot2)

#-----------------------------------------------------------------------
# function to create a dataframe that it will be used to generate scatterplot
# (compare edgeR_TCGA with edgeR or limma)

# method_toCompare: limma or edgeR
# dataframe: dataframe with genes detected by only edgeR_TCGAbiolinks in the rows
#            and logFC and FDR for each method in the columns
#-----------------------------------------------------------------------

dataframe_compare_edgeRTCGA <- function(method_toCompare,dataframe){
  if(method_toCompare=="limma"){
    FDR<- "FDR_limma"
    logFC <- "logFC_limma"
  }
  else{
    FDR <- "FDR_edgeR"
    logFC <- "logFC_edgeR"
  }
  dataframe$legend[dataframe[,FDR] >0.01 & dataframe[,logFC]<1 & dataframe[,logFC]>-1] <- "neither"
  dataframe$legend[dataframe[,FDR] >0.01 & dataframe[,logFC]>1]<- "logFC only"
  dataframe$legend[dataframe[,FDR] >0.01 & dataframe[,logFC]< -1]<- "logFC only"
  dataframe$legend[dataframe[,FDR] <0.01 & dataframe[,logFC]<1 & dataframe[,logFC]>-1]<- "FDR only"
  dataframe$legend[dataframe[,FDR] <0.01 & dataframe[,logFC]< -1]<- "down-regulated genes"
  
  return(dataframe)
}

#-----------------------------------------------------------------
# function to create dataframe used in the scatterplot to compare
# logFC and FDR between edgeR and limma
#-----------------------------------------------------------------

dataframe_compare_edgeR_limma <-function(dataframe){
  
  dataframe$legend[dataframe$FDR_edgeR >0.01 & dataframe$FDR_limma >0.01 & dataframe$logFC_edgeR<1 & dataframe$logFC_edgeR>-1 & dataframe$logFC_limma<1 & dataframe$logFC_limma>-1] <- "no DE"
  dataframe$legend[dataframe$FDR_edgeR >0.01 & dataframe$FDR_limma >0.01 & dataframe$logFC_edgeR>1 & dataframe$logFC_limma>1]<- "logFC only"
  dataframe$legend[dataframe$FDR_edgeR >0.01 & dataframe$FDR_limma >0.01 & dataframe$logFC_edgeR< -1 & dataframe$logFC_limma< -1]<- "logFC only"
  dataframe$legend[dataframe$FDR_edgeR <0.01 & dataframe$FDR_limma <0.01 & dataframe$logFC_edgeR<1 & dataframe$logFC_edgeR>-1 & dataframe$logFC_limma<1 & dataframe$logFC_limma>-1]<- "FDR only"
  dataframe$legend[dataframe$FDR_edgeR <0.01 & dataframe$FDR_limma<0.01 & dataframe$logFC_edgeR< -1 & dataframe$logFC_limma< -1]<- "down-regulated genes"
  dataframe$legend[is.na(dataframe$legend=="TRUE")] <- "other"
  return(dataframe)
}

#----------------------------------------------------------------------------------
# limma LUAD paired
#----------------------------------------------------------------------------------
up_LUAD <- read.csv("comparison_FC_FDR_up_paired_LUAD.csv",row.names = 1)
up_LUAD <- dataframe_compare_edgeRTCGA("limma",up_LUAD)
#scatterplot to compare logFC
#verify that you have created the folder scatterplots
pdf("./scatterplots/scatterplot_up_paired_LUAD_limma.pdf", width=10, height=7)
ggplot(up_LUAD, aes(x=logFC_edgeRTCGA,y=logFC_limma,colour=legend))+geom_point()+scale_colour_manual(values=c("red","dark green","dark turquoise","purple"))+
  xlab("logFC_edgeR-TCGAbiolinks")+ylab("logFC_limma")+
  theme(legend.title=element_blank(),axis.title=element_text(size=25),axis.text=element_text(size=25),
        legend.text=element_text(size=20))+ 
        guides(colour = guide_legend(override.aes = list(size=3)))+ylim(-2.5,3)
dev.off()
#scatterplot to compare FDR
png("./scatterplots/scatterplot_up_paired_LUAD_limma_FDR.png", width=700, height=500)
ggplot(up_LUAD,aes(x=FDR_edgeRTCGA,y=FDR_limma))+geom_point()+geom_hline(yintercept = 0.01)
dev.off()


#----------------------------------------------------------------------------------
# edgeR LUAD paired
#----------------------------------------------------------------------------------
up_LUAD <- read.csv("comparison_FC_FDR_up_paired_LUAD.csv",row.names = 1)
up_LUAD <- dataframe_compare_edgeRTCGA("edgeR",up_LUAD)
#scatterplot
pdf("./scatterplots/scatterplot_up_paired_LUAD_edgeR.pdf", width=10, height=7)
ggplot(up_LUAD, aes(x=logFC_edgeRTCGA,y=logFC_edgeR,colour=legend))+geom_point()+scale_colour_manual(values=c("red","dark green","dark turquoise","purple"))+
  xlab("logFC_edgeR-TCGAbiolinks")+ylab("logFC_edgeR")+
  theme(legend.title=element_blank(),axis.title=element_text(size=25),
        axis.text=element_text(size=25),legend.text=element_text(size=20))+ 
        guides(colour = guide_legend(override.aes = list(size=3)))+ylim(-3.5,3)
dev.off()
#scatterplot to compare FDR
png("./scatterplots/scatterplot_up_paired_LUAD_edgeR_FDR.png", width=700, height=500)
ggplot(up_LUAD,aes(x=FDR_edgeRTCGA,y=FDR_edgeR))+geom_point()
dev.off()

#-----------------------------------------------------------------------------------
# compare edgeR and limma LUAD paired
#-----------------------------------------------------------------------------------
up_LUAD <- read.csv("comparison_FC_FDR_up_paired_LUAD.csv",row.names = 1)
up_LUAD <- dataframe_compare_edgeR_limma(up_LUAD)
pdf("./scatterplots/scatterplot_up_paired_LUAD_edgeR-limma.pdf", width=10, height=7)
ggplot(up_LUAD, aes(x=logFC_edgeR,y=logFC_limma, colour=legend))+geom_point()+scale_colour_manual(values=c("red","dark green","dark turquoise","purple","dark grey"))+
  xlab("logFC_edgeR")+ylab("logFC_limma")+theme(legend.title=element_blank(),axis.title=element_text(size=25),axis.text=element_text(size=20),legend.text=element_text(size=20))+ guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()

#------------------------------------------------------------------------------------
# limma LUSC paired
#-----------------------------------------------------------------------------------
up_LUSC <- read.csv("comparison_FC_FDR_up_paired_LUSC.csv",row.names = 1)
up_LUSC <- dataframe_compare_edgeRTCGA("limma",up_LUSC)
#scatterplot
png("./scatterplots/scatterplot_up_paired_LUSC_limma.png", width=700, height=500)
ggplot(up_LUSC, aes(x=logFC_edgeRTCGA,y=logFC_limma,colour=legend))+geom_point()+scale_colour_manual(values=c("red","dark green","dark turquoise","purple"))+
  xlab("logFC_edgeR-TCGAbiolinks")+ylab("logFC_limma")+theme(legend.title=element_blank(),axis.title=element_text(size=25),axis.text=element_text(size=25),legend.text=element_text(size=20))+ guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()

#scatterplot to compare FDR
png("./scatterplots/scatterplot_up_paired_LUSC_limma_FDR.png", width=700, height=500)
ggplot(up_LUSC,aes(x=FDR_edgeRTCGA,y=FDR_limma))+geom_point()
dev.off()

#------------------------------------------------------------------------------------
# edgeR LUSC paired
#-----------------------------------------------------------------------------------
up_LUSC <- read.csv("comparison_FC_FDR_up_paired_LUSC.csv",row.names = 1)
up_LUSC <- dataframe_compare_edgeRTCGA("edgeR",up_LUSC)
#scatterplot
png("./scatterplots/scatterplot_up_paired_LUSC_edgeR.png", width=700, height=500)
ggplot(up_LUSC, aes(x=logFC_edgeRTCGA,y=logFC_edgeR,colour=legend))+geom_point()+scale_colour_manual(values=c("red","dark green","dark turquoise","purple"))+
  xlab("logFC_edgeR-TCGAbiolinks")+ylab("logFC_edgeR")+theme(legend.title=element_blank(),axis.title=element_text(size=25),axis.text=element_text(size=25),legend.text=element_text(size=20))+ guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()

#scatterplot to compare FDR
png("./scatterplots/scatterplot_up_paired_LUSC_edgeR_FDR.png", width=700, height=500)
ggplot(up_LUSC,aes(x=FDR_edgeRTCGA,y=FDR_edgeR))+geom_point()
dev.off()

#-----------------------------------------------------------------------------------
# compare edgeR and limma LUSC paired
#-----------------------------------------------------------------------------------
up_LUSC <- read.csv("comparison_FC_FDR_up_paired_LUSC.csv",row.names = 1)
up_LUSC <- dataframe_compare_edgeR_limma(up_LUSC)
png("./scatterplots/scatterplot_up_paired_LUSC_edgeR-limma.png", width=700, height=500)
ggplot(up_LUSC, aes(x=logFC_edgeR,y=logFC_limma, colour=legend))+geom_point()+scale_colour_manual(values=c("red","dark green","dark turquoise","purple","dark grey"))+
  xlab("logFC_edgeR")+ylab("logFC_limma")+theme(legend.title=element_blank(),axis.title=element_text(size=25),axis.text=element_text(size=20),legend.text=element_text(size=20))+ guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()

