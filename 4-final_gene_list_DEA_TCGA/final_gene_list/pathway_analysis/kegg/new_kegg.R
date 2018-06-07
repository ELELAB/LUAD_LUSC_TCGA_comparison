library(pathview)
library(gage)
library(biomaRt)

# ---------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR KEGG PATHWAY ENRICHMENT ANALYSIS

# 1)  my.keggset: list that includes Kegg-pathways with assigned gene IDs
#     (kegg.hsa.sigmet.gsets.RData loaded as kegg.gs)

# 2)  my.geneset: it includes the logFC values of the corrisponding genes IDs selected
#     as differentially expressed
# ---------------------------------------------------------------------------------------------------------------------------


my.keggs <- function(my.geneset, my.keggset) {
  keggres <- gage(my.geneset, gsets=my.keggset, same.dir=FALSE)
  keggpathways <- data.frame(id=rownames(keggres$greater), keggres$greater)
  keggpathways <- keggpathways[!is.na(keggpathways$p.val),]
  keggpathways <- as.character(keggpathways[keggpathways$p.val< 0.08,]$id)
  keggids <- substr(keggpathways, start=1, stop=8)
  pv.out.list <- sapply(keggids, function(pid) pathview(gene.data=my.geneset, pathway.id=pid, species="hsa", gene.idtype="KEGG", both.dirs = TRUE))
  return(pv.out.list)
}

#----------------------------------------------------------------------------------
# Function that gets the logFC values of the DE genes

# 1)  DE_genes: list of DE genes of interest
# 2)  fileName: name of the file that includes the logFC values
#----------------------------------------------------------------------------------
get.mygeneset <- function(DE_genes,fileName){
  
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  # to get the corrisponding entrez IDs of my DE genes symbols
  DE_gene_entrez <- getBM(attributes=c('hgnc_symbol', 'entrezgene'), filters = 'hgnc_symbol',
                          values = DE_genes , mart = ensembl)
  
  data <- read.csv(fileName)
  data$genesymbol <- data$X
  logFC <- data[data$genesymbol %in% DE_gene_entrez$hgnc_symbol,]
  new <- DE_gene_entrez[match(logFC$genesymbol, DE_gene_entrez$hgnc_symbol),]
  forkegg <- logFC$logFC
  names(forkegg) <- new$entrezgene
  return(forkegg)
}

#------------------------------------------------------------------------------------
# load keggsets
load("../kegg.hsa.sigmet.gsets.RData")

#------------------------------------------------------------------------------
# LUAD paired
#-----------------------------------------------------------------------------

up_paired_LUAD <- read.table("/data/user/marta/pipeline/DE/new_LUAD-LUSC/VennDiagram_three_methods/VennDiagram_LUAD-LUSC/final_genes_list/tumor_Purity/new_up_paired_LUAD.txt")
down_paired_LUAD <- read.table("/data/user/marta/pipeline/DE/new_LUAD-LUSC/VennDiagram_three_methods/VennDiagram_LUAD-LUSC/final_genes_list/tumor_Purity/new_down_paired_LUAD.txt")
DE_genes <- union(up_paired_LUAD$V1,down_paired_LUAD$V1)
fileName <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/paired/limma/limma_LUAD_paired_tumorPurity.csv"
forkegg <- get.mygeneset(DE_genes,fileName)
kegg_pws <- my.keggs(forkegg, kegg.gs)


#-----------------------------------------------------------------------------
# LUAD all
#-----------------------------------------------------------------------------

up_all_LUAD <- read.table("/data/user/marta/pipeline/DE/new_LUAD-LUSC/VennDiagram_three_methods/VennDiagram_LUAD-LUSC/final_genes_list/tumor_Purity/new_up_all_LUAD.txt")
down_all_LUAD <- read.table("/data/user/marta/pipeline/DE/new_LUAD-LUSC/VennDiagram_three_methods/VennDiagram_LUAD-LUSC/final_genes_list/tumor_Purity/new_down_all_LUAD.txt")
DE_genes <- union(up_all_LUAD$V1,down_all_LUAD$V1)
fileName <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/all/limma/limma_LUAD_all_tss_tumorPurity.csv"
forkegg <- get.mygeneset(DE_genes,fileName)
kegg_pws <- my.keggs(forkegg, kegg.gs)

#----------------------------------------------------------------------------
# LUSC all
#---------------------------------------------------------------------------
up_all_LUSC <- read.table("/data/user/marta/pipeline/DE/new_LUAD-LUSC/VennDiagram_three_methods/VennDiagram_LUAD-LUSC/final_genes_list/tumor_Purity/new_up_all_LUSC.txt")
down_all_LUSC <- read.table("/data/user/marta/pipeline/DE/new_LUAD-LUSC/VennDiagram_three_methods/VennDiagram_LUAD-LUSC/final_genes_list/tumor_Purity/new_down_all_LUSC.txt")
DE_genes <- union(up_all_LUSC$V1,down_all_LUSC$V1)
fileName <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/all/limma/limma_LUSC_all_tss_tumorPurity.csv"
forkegg <- get.mygeneset(DE_genes,fileName)
kegg_pws <- my.keggs(forkegg, kegg.gs)

#----------------------------------------------------------------------------------
# LUSC paired
#---------------------------------------------------------------------------------
up_paired_LUSC <- read.table("/data/user/marta/pipeline/DE/new_LUAD-LUSC/VennDiagram_three_methods/VennDiagram_LUAD-LUSC/final_genes_list/tumor_Purity/new_up_paired_LUSC.txt")
down_paired_LUSC <- read.table("/data/user/marta/pipeline/DE/new_LUAD-LUSC/VennDiagram_three_methods/VennDiagram_LUAD-LUSC/final_genes_list/tumor_Purity/new_down_paired_LUSC.txt")
DE_genes <- union(up_paired_LUSC$V1,down_paired_LUSC$V1)
fileName <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/paired/limma/limma_LUSC_paired_tumorPurity.csv"
forkegg <- get.mygeneset(DE_genes,fileName)
kegg_pws <- my.keggs(forkegg, kegg.gs)



