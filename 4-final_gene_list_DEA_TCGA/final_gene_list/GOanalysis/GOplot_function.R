library(GOplot)


#--------------------------------------------------------------------------------------
#Function that  makes genes dataframe to input in the circle_dat function

# fileName: name of file that includes logFC. It is usually a DEA output file
# genes_interest: a character vector that includes all the genes that belong to GO terms
# genes (output): A data frame with columns for 'ID' (genes symbol), 'logFC'

#--------------------------------------------------------------------------------------
make_genes_dataframe <- function(logFC_file,genes_interest){
  ID <- logFC_file$X
  logFC <- logFC_file$logFC
  genes <- data.frame(ID,logFC)
  genes <- subset(genes,genes$ID %in% genes_interest)
  rownames(genes) <- 1:nrow(genes)
  return(genes)
}

make_genes_column <- function(genes){ 
  c <- genes[1]
  if(length(genes)>=2){
    for(i in 2:length(genes)){
      c <- paste(c,genes[i],sep=",")
      c <- c
    }
  }
  return(c)
}

make_vector <- function(my.genes,GO.up,HUGO){
  vector <- c()
  for(i in 1:nrow(GO.up)){
    genes <- intersect(my.genes,HUGO[[GO.up$GO.ID[i]]])
    if(length(genes)==0)
      c <- 0
    else{
      c <- make_genes_column(genes)
    }
    vector <- c(vector,c)
  }
  return(vector)
}

make_genes_interest <- function(my.genes,GO.up,HUGO){
  genes_interest <- c()
  for(i in 1:nrow(GO.up)){
    genes <- intersect(my.genes,HUGO[[GO.up$GO.ID[i]]])
    genes_interest <- c(genes_interest,genes)
  }
  return(genes_interest)
}

#------------------------------------------------------------------------------------
# make terms dataframe to input in the circle_dat function
#-------------------------------------------------------------------------------------
make_terms_dataframe <- function(GO.up,vector, category){
  category <- category
  ID <- GO.up$GO.ID
  term <- GO.up$Term
  adj_pval <- GO.up$FDR_fisher
  genes <- vector
  terms <- data.frame(category,ID,term,adj_pval,genes)
  # remove the terms that not include any genes
  terms <- subset(terms,!terms$genes==0)
  return(terms)
}
