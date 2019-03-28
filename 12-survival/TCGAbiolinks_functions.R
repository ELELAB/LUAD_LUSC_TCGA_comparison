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
# FUNCTION FOR IDENTIFICATION OF SIGNIFICANTLY DIFFERENTIALLY EXPRESSED GENES - EDGER 
# -----------------------------------------------------------------------------------------------------------------------------


DE_edgeR <- function(my.lrt, my.data, coLFC, coFDR) {
  my.tags <- topTags(my.lrt, n=nrow(my.data$counts))
  my.tags <- my.tags$table
  
  index.up <- which(my.tags$logFC >= coLFC & my.tags$FDR < coFDR)
  index.down <- which(my.tags$logFC <= -coLFC & my.tags$FDR < coFDR)
  direction <- c()
  direction[index.up] <- "up"
  direction[index.down] <- "down"
  direction[!(1:nrow(my.tags) %in% union(index.up,index.down))] <- "no DE"
  my.tags <- cbind(my.tags,direction)
  
  return(my.tags)
}


# -----------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR IDENTIFICATION OF SIGNIFICANTLY DIFFERENTIALLY EXPRESSED GENES - LIMMA
# -----------------------------------------------------------------------------------------------------------------------------


DE_limma <- function(my.contrast, my.data, my.design, coLFC, coFDR) {
  fit3 <- eBayes(contrasts.fit(lmFit(my.data, my.design), my.contrast))
  tt <- toptable(fit3, coef=1, adjust='fdr', number=nrow(my.data))
  index.up <- which(tt$logFC >= coLFC & tt$adj.P.Val < coFDR)
  index.down <- which(tt$logFC <= -coLFC & tt$adj.P.Val < coFDR)
  direction <- c()
  direction[index.up] <- "up"
  direction[index.down] <- "down"
  direction[!(1:nrow(tt) %in% union(index.up,index.down))] <- "no DE"
  tt <- cbind(tt,direction)

  return(tt)
}


#------------------------------------------------------------------------

#FUNCTION THAT REMOVE THE DUPLICATED TUMOR SAMPLES AND REPLACE THEM WITH
# THEIR MEAN

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
# THEIR MEDIAN

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

# RANDOM REMOVING OF TUMOR REPLICATES

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

# limma's function for unified data (tcga+gtex); tss-correction 

limma_unified <- function(dataframe, condition, tss, limma_name,up_name, down_name){
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

