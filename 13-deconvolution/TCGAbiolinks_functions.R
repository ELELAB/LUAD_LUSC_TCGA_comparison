# -----------------------------------------------------------------------------------------------------------------------------
# LOADING PACKAGES
# -----------------------------------------------------------------------------------------------------------------------------

library(TCGAbiolinks)
library(plyr)
library(ggplot2)
library(pamr)
library(limma)
library(sva)
library(edgeR)
library(EDASeq)
library(heatmap.plus)
library(VennDiagram)



# -----------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR EXTRACTING INFORMATION ON BATCH, PATIENT AND ID.
# -----------------------------------------------------------------------------------------------------------------------------

get_IDs <- function(data) {
  IDs <- strsplit(c(colnames(data)), "-")
  IDs <- ldply(IDs, rbind)
  colnames(IDs) <- c('project', 'tss','participant', 'sample', "portion", "plate", "center")
  cols <- c("project", "tss", "participant")
  IDs$patient <- apply(IDs[,cols],1,paste,collapse = "-" )
  barcode <- colnames(data)
  IDs <- cbind(IDs, barcode)
  condition <- gsub("11+[[:alpha:]]", "normal", as.character(IDs$sample))
  condition  <- gsub("01+[[:alpha:]]", "cancer", condition)
  IDs$condition <- condition
  IDs$myorder  <- 1:nrow(IDs)
  return(IDs)
}




# -----------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR EMDS-PLOTTING
# -----------------------------------------------------------------------------------------------------------------------------


myMDSplot <- function(my.data, my.group, my.labels) {
  d<-dist(t(my.data))
  fit <- cmdscale(d,eig=TRUE, k=2)
  res<-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
  p <- ggplot(data=res)
  p + geom_point(aes(x=M1,y=M2,color=my.group)) + geom_text(aes(x=M1,y=M2, label= my.labels, color=my.group)) +
    coord_cartesian(xlim=c(min(res$M1)*1.4,max(res$M1)*1.4)) + theme_gray() + theme(legend.text = element_text(size = 16, face="bold"), axis.title=element_text(size=16,face="bold"), axis.text.x=element_blank()) + guides(colour = guide_legend(override.aes = list(size=6))) + theme(legend.position = "top") + theme(axis.text=element_text(size=16, face="bold")) + theme(axis.text = element_text(colour = "black"))
}



# -----------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR IDENTIFICATION OF SIGNIFICANTLY DIFFERENTIALLY EXPRESSED GENES - EDGER 
# -----------------------------------------------------------------------------------------------------------------------------


DE_edgeR <- function(my.lrt, my.data, coLFC, coFDR) {
  my.tags <- topTags(my.lrt, n=nrow(my.data$counts))
  my.tags <- my.tags$table
  
  #my.tags <- subset(my.tags,abs(my.tags$logFC) >= coLFC & my.tags$FDR < coFDR)
  index.up <- which(my.tags$logFC >= coLFC & my.tags$FDR < coFDR)
  index.down <- which(my.tags$logFC <= -coLFC & my.tags$FDR < coFDR)
  direction <- c()
  direction[index.up] <- "up"
  direction[index.down] <- "down"
  direction[!(1:nrow(my.tags) %in% union(index.up,index.down))] <- "no DE"
  my.tags <- cbind(my.tags,direction)
  
  #up <- data.frame(rownames(my.tags[my.tags$logFC >= coLFC & my.tags$FDR < coFDR, ]))
  #down <- data.frame(rownames(my.tags[my.tags$logFC <= -coLFC & my.tags$FDR < coFDR, ]))
  
  #colnames(up) <- as.character("up")
  #colnames(down) <- as.character("down")
  
  #final <- list(up, down)
  return(my.tags)
}


# -----------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR IDENTIFICATION OF SIGNIFICANTLY DIFFERENTIALLY EXPRESSED GENES - LIMMA
# -----------------------------------------------------------------------------------------------------------------------------


DE_limma <- function(my.contrast, my.data, my.design, coLFC, coFDR) {
  fit3 <- eBayes(contrasts.fit(lmFit(my.data, my.design), my.contrast))
  tt <- toptable(fit3, coef=1, adjust='fdr', number=nrow(my.data))
  
  #up <- data.frame(rownames(tt[tt$logFC >= coLFC & tt$adj.P.Val < coFDR, ]))
  #down <- data.frame(rownames(tt[tt$logFC <= -coLFC & tt$adj.P.Val < coFDR, ]))
  
  #colnames(up) <- as.character("up")
  #colnames(down) <- as.character("down")
  
  #tt <- subset(tt,abs(tt$logFC) >= coLFC & tt$adj.P.Val < coFDR)
  index.up <- which(tt$logFC >= coLFC & tt$adj.P.Val < coFDR)
  index.down <- which(tt$logFC <= -coLFC & tt$adj.P.Val < coFDR)
  direction <- c()
  direction[index.up] <- "up"
  direction[index.down] <- "down"
  direction[!(1:nrow(tt) %in% union(index.up,index.down))] <- "no DE"
  tt <- cbind(tt,direction)
 
   #final <- list(up, down)
  
  return(tt)
}



# -----------------------------------------------------------------------------------------------------------------------------
# FUNCTION THAT GENERATES A COLOR MATRIX FOR HEATMAP.PLUS
# -----------------------------------------------------------------------------------------------------------------------------

get_colors <- function(my.truestatus) {
  hm_col <- data.frame(c("normal", "Basal-like", "HER2-enriched", "Luminal A", "Luminal B", "Normal-like", "cancer"), c("#FD61D1", "#B983FF", "#00BA38", "#00B0F6", "#F8766D", "#E58700", "#00BA38"))
  colnames(hm_col) <- c("status", "mycolor")
  true_status <- data.frame(my.truestatus)
  myorder <- 1:nrow(true_status)
  true_status$order <- myorder
  colnames(true_status) <- c("status", "order")
  col <- merge(true_status, hm_col, by="status", all.x =TRUE)
  col <-col[order(col$order),]
  col <- replicate(2,col$mycolor)
  return(col)
}

#------------------------------------------------------------------------

#FUNCTION THAT REMOVE THE DUPLICATED TUMOR SAMPLES AND REPLACE THEM WITH
# THEIR MEAN (Marta)

# dataframe = PreprocessedData.rda file
# my_IDs = output of get_IDs function
#-------------------------------------------------------------------------

mean.duplicated.tumor <- function(dataframe,my_IDs){

# tumor.samples includes only the tumor samples
  tumor.samples <- my_IDs[my_IDs$condition == "cancer", ]
# tumor.replicates dataframe includes only the tumor duplicates
  tumor.replicates <- tumor.samples[duplicated(tumor.samples$patient, fromLast=FALSE) | duplicated(tumor.samples$patient, fromLast=TRUE),]
# find the patients' names of tumor duplicates
  tumor.replicates.patient <- as.factor(tumor.replicates$patient)
  
  for(i in 1:length(levels(tumor.replicates.patient))){
  # find the index of tumor replicates for each patient
    index <- which(tumor.replicates.patient==levels(tumor.replicates.patient)[i])
    vector.replicates <- as.vector(tumor.replicates$barcode[index]) # find the barcodes of tumor duplicates for each patient
  # make a new matrix (nrow= genes number,ncol=1) that includes the mean between tumor duplicates
    new.matrix <- as.matrix(rowMeans(dataframe[,colnames(dataframe) %in% vector.replicates]))
    rownames(new.matrix) <- rownames(dataframe)
    colnames(new.matrix) <- as.character(vector.replicates[1])
  # remove the tumor replicates from dataframe
    dataframe <- dataframe[,!colnames(dataframe) %in% vector.replicates]
    dataframe <- cbind(new.matrix, dataframe)# add the new column
  }
  
  return (dataframe)
  
}

#------------------------------------------------------------------------

#FUNCTION THAT REMOVE THE DUPLICATED TUMOR SAMPLES AND REPLACE THEM WITH
# THEIR MEDIAN (Marta)

# dataframe = PreprocessedData.rda file
# my_IDs = output of get_IDs function
#-------------------------------------------------------------------------

median.duplicated.tumor <- function(dataframe,my_IDs){
  
  # tumor.samples includes only the tumor samples
 tumor.samples <- my_IDs[my_IDs$condition == "cancer", ]
  #tumor.replicates includes only the tumor duplicates
  tumor.replicates <- tumor.samples[duplicated(tumor.samples$patient, fromLast=FALSE) | duplicated(tumor.samples$patient, fromLast=TRUE),]
  # find the patients' names of tumor duplicates
  tumor.replicates.patient <- as.factor(tumor.replicates$patient)
  
  for(i in 1:length(levels(tumor.replicates.patient))){
    # find the indeces of tumor replicates for each patient
    index <- which(tumor.replicates.patient==levels(tumor.replicates.patient)[i])
    vector.replicates <- as.vector(tumor.replicates$barcode[index])
    # make a new matrix (nrow= genes number,ncol=1) that includes the median between tumor duplicates
    new.matrix <- as.matrix(rowMedians(dataframe[,colnames(dataframe) %in% vector.replicates]))
    rownames(new.matrix) <- rownames(dataframe)
    colnames(new.matrix) <- as.character(vector.replicates[1])
    # remove the tumor replicates from dataframe
    dataframe <- dataframe[,!colnames(dataframe) %in% vector.replicates]
    dataframe <- cbind(new.matrix, dataframe) # add the new column
  }
  
  return (dataframe)
  
}

#--------------------------------------------------------------------

# RANDOM REMOVING OF TUMOR REPLICATES (Marta)

# dataframe = PreprocessedData.rda file
# my_IDs = output of get_IDs function
#--------------------------------------------------------------------

random.removing <- function(my_IDs,dataframe){

  #tumor.replicates includes only the tumor duplicates
  tumor.samples <- my_IDs[my_IDs$condition == "cancer", ]
  tumor.replicates <- tumor.samples[duplicated(tumor.samples$patient, fromLast=FALSE) | duplicated(tumor.samples$patient, fromLast=TRUE),]
  # find the patients' names of tumor duplicates
  tumor.replicates.patient <- as.factor(tumor.replicates$patient)
  
  for(i in 1:length(levels(tumor.replicates.patient))){
    # find the index of tumor replicates for each patient
    index <- which(tumor.replicates.patient==levels(tumor.replicates.patient)[i])
    vector.replicates <- as.vector(tumor.replicates$barcode[index]) # find the barcodes of tumor duplicates
  # get (n-1) random indices where n= number of tumor duplicates  
    index.random <- sample(1:length(vector.replicates),length(vector.replicates)-1)
    #remove n-1 tumor replicates from dataframe
    dataframe <- dataframe[,!colnames(dataframe) %in% vector.replicates[index.random]]
  }
  
  return (dataframe)
  
}


#--------------------------------------------------------------------

# VENN-DIAGRAM comparison between all,paired,no-paired 
#into a same cancer type

#----------------------------------------------------------------------

Venn_Diagram <- function(cancerType, direction, method){
dir <- "/data/user/marta/pipeline/DE/limma/"

paired <- read.table(paste0(dir,cancerType,"/",direction,"_",method,"_paired_",cancerType,".txt"),col.names = "genes")
all <- read.table(paste0(dir,cancerType,"/",direction,"_",method,"_all_",cancerType,".txt"),col.names = "genes")
nonpaired <- read.table(paste0(dir,cancerType,"/",direction,"_",method,"_nonpaired_",cancerType,".txt"),col.names = "genes")

png(paste0(dir,"VennDiagram/","vennDiagram_",direction,"_",cancerType,".png"), height = 500, width = 600)
venn <- venn.diagram(list(paired$genes, all$genes, nonpaired$genes), 
                     category.names = c("paired","all","non-paired"),
                     filename = NULL, lwd = 0.5,
                     main = (paste0(direction,"-regulated genes in ",cancerType,"-project")),
                     fill=c(4,2,3),main.cex=2,sub.cex = 1.5, cat.cex=1.5, cex=1.5)
grid.draw(venn)

dev.off()


overlap <- intersect(paired$genes,all$genes)
overlap <- intersect(overlap,nonpaired$genes)
write.table(overlap,paste0(dir,"VennDiagram/","overlap_",direction,"_",method,"_",cancerType,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
total_genes <- union(paired$genes,all$genes)
total_genes <- union(total_genes,nonpaired$genes)
no_overlap <- total_genes[!(total_genes %in% overlap)]
write.table(no_overlap,paste0(dir,"VennDiagram/","no_overlap_",direction,"_",method,"_",cancerType,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)

return(print("Venn Diagram done!!"))

}

#--------------------------------------------------------------------

# VENN-DIAGRAM comparison between two cancerTypes

#----------------------------------------------------------------------

Venn_Diagram2 <- function(cancerType1,cancerType2,method,direction,datasetType){

  dir <- paste0("/data/user/marta/pipeline/DE/",method,"/")
  
  list1 <- read.table(paste0(dir,cancerType1,"/",direction,"_",method,"_",datasetType,"_",cancerType1,".txt"),col.names = "genes")
  list2 <- read.table(paste0(dir,cancerType2,"/",direction,"_",method,"_",datasetType,"_",cancerType2,".txt"),col.names = "genes")
  
  
  png(paste0(dir,"VennDiagram/","vennDiagram_",method,"_",direction,"_",datasetType,"_",cancerType1,"-",cancerType2,".png"), height = 500, width = 700)
  venn <- venn.diagram(list(list1$genes, list2$genes), 
                       category.names = c(paste0(datasetType,"-",cancerType1),paste0(datasetType,"-",cancerType2)),
                       filename = NULL, lwd = 0.5,
                       main = (paste0(direction,"-regulated genes in ",cancerType1,"-",cancerType2)),
                       fill=c(4,2),main.cex=2,sub.cex = 1.5, cat.cex=1, cex=1)
  grid.draw(venn)
  
  dev.off()
  
  
  overlap <- intersect(list1$genes,list2$genes)
  write.table(overlap,paste0(dir,"VennDiagram/","overlap_",direction,"_",method,"_",datasetType,"_",cancerType1,"-",cancerType2,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  total_genes <- union(list1$genes,list2$genes)
  no_overlap <- total_genes[!(total_genes %in% overlap)]
  write.table(no_overlap,paste0(dir,"VennDiagram/","no_overlap_",direction,"_",method,"_",datasetType,"_",cancerType1,"-",cancerType2,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  return(print("Venn Diagram done!!"))
  
  }

#--------------------------------------------------------------------

# VENN-DIAGRAM comparison between two methods (limma-edgeR)

#----------------------------------------------------------------------

Venn_Diagram3 <- function(cancerType,method1,method2,direction,datasetType,index=""){
  
  dir1 <- paste0("/data/user/marta/pipeline/DE/",method1,"/")
  dir2 <- paste0("/data/user/marta/pipeline/DE/",method2,"/")
  
  list1 <- read.table(paste0(dir1,cancerType,"/",direction,"_",method1,"_",datasetType,"_",cancerType,".txt"),col.names = "genes")
  list2 <- read.table(paste0(dir2,cancerType,"/",direction,"_",method2,"_",datasetType,index,"_",cancerType,".txt"),col.names = "genes")
  
  dir <- paste0("/data/user/marta/pipeline/DE/")
  png(paste0(dir,"VennDiagram_comparison/","vennDiagram_",method1,"-",method2,"_",direction,"_",datasetType,index,"_",cancerType,".png"), height = 500, width = 700)
  venn <- venn.diagram(list(list1$genes, list2$genes), 
                       category.names = c(paste0(datasetType,"-",method1),paste0(datasetType,index,"-",method2)),
                       filename = NULL, lwd = 0.5,
                       main = (paste0(direction,"-regulated genes in ",cancerType,"-",datasetType,"(",method1,"-",method2,")")),
                       fill=c(4,2),main.cex=2,sub.cex = 1.5, cat.cex=1, cex=1)
  grid.draw(venn)
  
  dev.off()
  
  
  overlap <- intersect(list1$genes,list2$genes)
  write.table(overlap,paste0(dir,"VennDiagram_comparison/","overlap_",direction,"_",method1,"-",method2,"_",datasetType,index,"_",cancerType,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  total_genes <- union(list1$genes,list2$genes)
  no_overlap <- total_genes[!(total_genes %in% overlap)]
  write.table(no_overlap,paste0(dir,"VennDiagram_comparison/","no_overlap_",direction,"_",method1,"-",method2,"_",datasetType,index,"_",cancerType,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  return(print("Venn Diagram done!!"))
  
}




#---------------------------------------------------------------------

# FUNCTION that checks and removes the LUSC discordant samples (Marta)  

#----------------------------------------------------------------------

discordant_removing <- function(discordant.samples,dataframe){
  
  index <- which(colnames(dataframe) %in% discordant.samples$barcode)
  print(paste0("The dataset includes ",length(index)," discordant LUSC samples"))
  print(colnames(dataframe)[index])
  if(length(index)!=0){
    dataframe <- dataframe[ ,-(index)]
    print("Removing done!!")
    print(paste0("The total samples number is ", ncol(dataframe)))
  }
  else{
    print("no discordant samples into dataset")
  }
  return(dataframe)
}
