source("GOplot_function.R")
source("GeneOntology_functions.R")
load("~/Desktop/PhD_CNR/co-expression_April2019/CorePeel_moduliMarco_22-03-19/3-GO/HUGO.RData")

#-----------------------------------------------------------------------------
#ont= "BP" or "CC", "MF"
#nameGO: topGO output name
#logFC_file: file that includes logFC (for example from limma analysis)

GOanalysis_plot <- function(ont,my.gene.univers,my.DE.genes,HUGO,nTerm,nameGO,logFC_file,plot){
  
  GO <- TOPGO(ont, my.gene.univers, my.DE.genes, HUGO, nTerm, nameGO)
  vector <- make_vector(my.DE.genes,GO,HUGO)
  genes_interest <- make_genes_interest(my.DE.genes,GO,HUGO)
  terms <- make_terms_dataframe(GO,vector,ont)
  genes <- make_genes_dataframe(logFC_file,genes_interest)
  circ <- circle_dat(terms, genes)
  process <- terms$term
  chord <- chord_dat(data = circ, genes, process)
  png(filename=plot,height = 1000, width = 1500)
  GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4,lfc.min = min(genes$logFC), lfc.max = max(genes$logFC))
  #dev.off()
  
  #return(dev.off())
}


GO_analysis_deconvolution <- function(cancerType,populationType, modulation){
  
  dataFilt <- get(load(paste0("../",cancerType,"_PreprocessedData_all_TumorPurity.rda")))
  my.gene.univers <- rownames(dataFilt)
  genes <- read.table(paste0("../1-DEA/",populationType,"/",modulation,"_",cancerType,"_only.txt"))
  genes <- genes$V1
  
  logFC_file <- read.csv(paste0('limma_',cancerType,'_all_tss_tumorPurity.csv'))
  plot <- paste0('circle_BP_',modulation,'_',populationType,'_',cancerType,'.png')
  GOanalysis_plot("BP",my.gene.univers,genes,HUGO,20,paste0("GO.",modulation,'.',populationType,'.',cancerType),logFC_file,plot)
}

##################

GO_analysis_deconvolution('LUAD','B_lineage','up')
dev.off()

GO_analysis_deconvolution('LUAD','B_lineage','down')
dev.off()

GO_analysis_deconvolution('LUSC','B_lineage','up')
dev.off()

GO_analysis_deconvolution('LUSC','B_lineage','down')
dev.off()

#------------------------------

GO_analysis_deconvolution('LUAD','Endothelial_cells','up')
dev.off()

GO_analysis_deconvolution('LUAD','Endothelial_cells','down')
dev.off()

GO_analysis_deconvolution('LUSC','Endothelial_cells','up')
dev.off()

GO_analysis_deconvolution('LUSC','Endothelial_cells','down')
dev.off()

#---------------------------

GO_analysis_deconvolution('LUAD','Fibroblasts','up')
dev.off()

GO_analysis_deconvolution('LUAD','Fibroblasts','down')
dev.off()

GO_analysis_deconvolution('LUSC','Fibroblasts','up')
dev.off()

GO_analysis_deconvolution('LUSC','Fibroblasts','down')
dev.off()


#------------------------------------

GO_analysis_deconvolution('LUAD','Lymphocytes','up')
dev.off()

GO_analysis_deconvolution('LUAD','Lymphocytes','down')
dev.off()

GO_analysis_deconvolution('LUSC','Lymphocytes','up')
dev.off()

GO_analysis_deconvolution('LUSC','Lymphocytes','down')
dev.off()

#------------------------------------------

GO_analysis_deconvolution('LUAD','Monocytic','up')
dev.off()

GO_analysis_deconvolution('LUAD','Monocytic','down')
dev.off()

GO_analysis_deconvolution('LUSC','Monocytic','up')
dev.off()

GO_analysis_deconvolution('LUSC','Monocytic','down')
dev.off()

#-----------------------------------
GO_analysis_deconvolution('LUAD','Myeloid_cells','up')
dev.off()

GO_analysis_deconvolution('LUAD','Myeloid_cells','down')
dev.off()

GO_analysis_deconvolution('LUSC','Myeloid_cells','up')
dev.off()

GO_analysis_deconvolution('LUSC','Myeloid_cells','down')
dev.off()


#---------------------------------------------

GO_analysis_deconvolution('LUAD','Neutrophils','up')
dev.off()

GO_analysis_deconvolution('LUAD','Neutrophils','down')
dev.off()

GO_analysis_deconvolution('LUSC','Neutrophils','up')
dev.off()

GO_analysis_deconvolution('LUSC','Neutrophils','down')
dev.off()

#-------------------------------------------
GO_analysis_deconvolution('LUAD','Tcells','up')
dev.off()

GO_analysis_deconvolution('LUAD','Tcells','down')
dev.off()

GO_analysis_deconvolution('LUSC','Tcells','up')
dev.off()

GO_analysis_deconvolution('LUSC','Tcells','down')
dev.off()
