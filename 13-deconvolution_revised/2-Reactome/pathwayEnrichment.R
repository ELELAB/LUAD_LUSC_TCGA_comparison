library(ReactomePA)
library(biomaRt)
library(AnnotationDbi)
library(S4Vectors)
library(IRanges)
library(SummarizedExperiment)

# --------------------------------------------------------------------------------------
# PATHWAY ENRICHMENT ANALYSIS

#background: dataframe that includes the universe gene entrez-IDs and the corrisponding
# gene symbols
# my.genes: genes selected as DE
# my.name: the name of the output file
# my.plot: logical (TRUE if you want a output plot)
# value: p-value cutoff (default: 0.05)
# --------------------------------------------------------------------------------------

enrich_pathway <- function(background, my.genes, my.name, my.plot, value=0.05) {
  entrez.data <- background[background$geneName %in% my.genes, ]
  my.pathway <- enrichPathway(as.character(entrez.data$entrez), organism = "human", pvalueCutoff = value, pAdjustMethod = "fdr", qvalueCutoff = value, as.character(background$entrez), minGSSize = 3, readable = T)
  if (my.plot == TRUE) {
    png(filename = paste0(my.name,".png"), height = 1500, width = 1500)
    cnetplot(my.pathway, categorySize="geneNum")
    dev.off()
  }
  return(my.pathway)
}

#------------------------------------------------------------------------------
# function that performs pathway enrichment analysis
# cancerType= LUAD or LUSC
# populationType= B_lineage, Tcells and so on
#-------------------------------------------------------------------------------

pathway_deconvolution <- function(cancerType,populationType){
  dataframe <- get(load(paste0("../",cancerType,"_PreprocessedData_all_TumorPurity.rda")))
  SE <- get(load(paste0("../",cancerType,"_Illumina_HiSeq_all_tumorPurity.rda")))
  background <- subset(rowData(SE), rowData(SE)$gene_id %in% rownames(dataframe), select=c("gene_id","entrezgene"))
  colnames(background) <- c("geneName","entrez")


  up <- read.table(paste0("../1-DEA/",populationType,"/up_",cancerType,"_only.txt"))
  up_genes <- up$V1
  down<- read.table(paste0("../1-DEA/",populationType,"/down_",cancerType,"_only.txt"))
  down_genes <- down$V1

  pathway_up <- enrich_pathway(background,up_genes,paste0(populationType,"/pathway_up_",cancerType,"_",populationType),TRUE,0.05)
  write.csv(pathway_up, paste0(populationType,"/pathway_table_up_",cancerType,"_",populationType,".csv"),row.names=FALSE)
  
  pathway_down <- enrich_pathway(background,down_genes,paste0(populationType,"/pathway_down_",cancerType,"_",populationType),TRUE,0.05)
  write.csv(pathway_down, paste0(populationType,"/pathway_table_down_",cancerType,"_",populationType,".csv"),row.names=FALSE)
  

}


########### B_lineage ###############

pathway_deconvolution('LUAD','B_lineage')
pathway_deconvolution('LUSC','B_lineage')


######### Endothelial_cells ############

pathway_deconvolution('LUAD','Endothelial_cells')
pathway_deconvolution('LUSC','Endothelial_cells')

######### Tcells ############

pathway_deconvolution('LUAD','Tcells')
pathway_deconvolution('LUSC','Tcells')

######### Fibroblasts ############

pathway_deconvolution('LUAD','Fibroblasts')
pathway_deconvolution('LUSC','Fibroblasts')

######### Lymphocytes ############   luad down !!!!

pathway_deconvolution('LUAD','Lymphocytes')
pathway_deconvolution('LUSC','Lymphocytes')

######### Monocytic ########### luad down !!!!

pathway_deconvolution('LUAD','Monocytic') 
pathway_deconvolution('LUSC','Monocytic')

######### Myeloid_cells ############ luad down !!!!

pathway_deconvolution('LUAD','Myeloid_cells')
pathway_deconvolution('LUSC','Myeloid_cells')

######### Neutrophils ############

pathway_deconvolution('LUAD','Neutrophils')
pathway_deconvolution('LUSC','Neutrophils')


setwd("~/Desktop/13-deconvolution/2-Reactome/Pathway_0.08")

#------------------------------------------------------------------------------
# function that performs pathway enrichment analysis with a higher fdr (0.08)
# cancerType= LUAD or LUSC
# populationType= B_lineage, Tcells and so on
#-------------------------------------------------------------------------------

pathway_deconvolution_0.08 <- function(cancerType,populationType){
  dataframe <- get(load(paste0("../../",cancerType,"_PreprocessedData_all_TumorPurity.rda")))
  SE <- get(load(paste0("../../",cancerType,"_Illumina_HiSeq_all_tumorPurity.rda")))
  background <- subset(rowData(SE), rowData(SE)$gene_id %in% rownames(dataframe), select=c("gene_id","entrezgene"))
  colnames(background) <- c("geneName","entrez")
  
  
  up <- read.table(paste0("../../1-DEA/",populationType,"/up_",cancerType,"_only.txt"))
  up_genes <- up$V1
  down<- read.table(paste0("../../1-DEA/",populationType,"/down_",cancerType,"_only.txt"))
  down_genes <- down$V1
  
  pathway_up <- enrich_pathway(background,up_genes,paste0(populationType,"/pathway_up_",cancerType,"_",populationType),TRUE,0.08)
  write.csv(pathway_up, paste0(populationType,"/pathway_table_up_",cancerType,"_",populationType,".csv"),row.names=FALSE)
  
  pathway_down <- enrich_pathway(background,down_genes,paste0(populationType,"/pathway_down_",cancerType,"_",populationType),TRUE,0.08)
  write.csv(pathway_down, paste0(populationType,"/pathway_table_down_",cancerType,"_",populationType,".csv"),row.names=FALSE)
  
  
}

########### B_lineage ###############

pathway_deconvolution_0.08('LUAD','B_lineage')
pathway_deconvolution_0.08('LUSC','B_lineage')


######### Endothelial_cells ############

pathway_deconvolution_0.08('LUAD','Endothelial_cells')
pathway_deconvolution_0.08('LUSC','Endothelial_cells')

######### Tcells ############

pathway_deconvolution_0.08('LUAD','Tcells')
pathway_deconvolution_0.08('LUSC','Tcells')

######### Fibroblasts ############

pathway_deconvolution_0.08('LUAD','Fibroblasts')
pathway_deconvolution_0.08('LUSC','Fibroblasts')

######### Lymphocytes ############   luad down !!!!

pathway_deconvolution_0.08('LUAD','Lymphocytes')
pathway_deconvolution_0.08('LUSC','Lymphocytes')

######### Monocytic ########### luad down !!!!

pathway_deconvolution_0.08('LUAD','Monocytic') 
pathway_deconvolution_0.08('LUSC','Monocytic')

######### Myeloid_cells ############ luad down !!!!

pathway_deconvolution_0.08('LUAD','Myeloid_cells')
pathway_deconvolution_0.08('LUSC','Myeloid_cells')

######### Neutrophils ############

pathway_deconvolution_0.08('LUAD','Neutrophils')
pathway_deconvolution_0.08('LUSC','Neutrophils')
