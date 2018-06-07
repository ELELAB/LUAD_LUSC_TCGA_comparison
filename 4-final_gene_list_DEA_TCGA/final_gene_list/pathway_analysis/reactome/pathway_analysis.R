library(ReactomePA)
library(biomaRt)
library(AnnotationDbi)
library(S4Vectors)
library(IRanges)
library(SummarizedExperiment)
library(plyr)

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
  my.pathway <- enrichPathway(as.character(entrez.data$entrez), organism = "human", pvalueCutoff = value, pAdjustMethod = "fdr", qvalueCutoff = value, as.character(background$entrez), minGSSize = 2, readable = T)
  if (my.plot == TRUE) {
    png(filename = paste0(my.name,".png"), height = 1000, width = 1000)
    cnetplot(my.pathway, categorySize="geneNum", showCategory = 2)
    dev.off()
  }
  return(my.pathway)
} 


#-----------------------------------------------------------------------------------
#LUAD all
#-----------------------------------------------------------------------------------
dataframe_all_LUAD <- get(load("../../../../1-download_preprocessing/LUAD/all/LUAD_PreprocessedData_all_tumorPurity.rda"))
# use SummarizedExperiment object to get the entrez IDs
SE_LUAD <- get(load("../../../../1-download_preprocessing/LUAD/all/LUAD_Illumina_HiSeq_all_tumorPurity.rda"))
background_all_LUAD <- subset(rowData(SE_LUAD), rowData(SE_LUAD)$gene_id %in% rownames(dataframe_all_LUAD), select=c("gene_id","entrezgene"))
colnames(background_all_LUAD) <- c("geneName","entrez")
#up
up_all_LUAD <- read.table("../../new_up_all_LUAD.txt")
up_all_LUAD_genes <- up_all_LUAD$V1
# pathway analysis
pathway_up_all_LUAD <- enrich_pathway(background_all_LUAD,up_all_LUAD_genes,"pathway_up_all_LUAD_tumorPurity",TRUE,0.05)
write.csv(pathway_up_all_LUAD, "pathway_table_up_all_LUAD_tumorPurity.csv",row.names=FALSE)

#down
down_all_LUAD <- read.table("../..//new_down_all_LUAD.txt")
down_all_LUAD_genes <- down_all_LUAD$V1
subset <- background_all_LUAD[background_all_LUAD$geneName %in% down_all_LUAD_genes, ]
pathway_down_all_LUAD <- enrich_pathway(background_all_LUAD,down_all_LUAD_genes,"pathway_down_all_LUAD_tumorPurity",TRUE,0.05)
write.csv(pathway_down_all_LUAD, "pathway_table_down_all_LUAD_tumorPurity.csv",row.names = FALSE)

#--------------------------------------------------------------------------------
# LUSC all
#-------------------------------------------------------------------------------
dataframe_all_LUSC <- get(load("../../../../1-download_preprocessing/LUSC/all/LUSC_PreprocessedData_all_tumorPurity.rda"))
SE_LUSC <- get(load("../../../../1-download_preprocessing/LUSC/all/LUSC_Illumina_HiSeq_all_tumorPurity.rda"))
background_all_LUSC <- subset(rowData(SE_LUSC), rowData(SE_LUSC)$gene_id %in% rownames(dataframe_all_LUSC), select=c("gene_id","entrezgene"))
colnames(background_all_LUSC) <- c("geneName","entrez")
#background_all_LUSC <- conversion.background(dataframe_all_LUSC)
#up
up_all_LUSC <- read.table("../../new_up_all_LUSC.txt")
up_all_LUSC_genes <- up_all_LUSC$V1
subset <- background_all_LUSC[background_all_LUSC$geneName %in% up_all_LUSC_genes, ]
pathway_up_all_LUSC <- enrich_pathway(background_all_LUSC,up_all_LUSC_genes,"pathway_up_all_LUSC_tumorPurity",TRUE,0.05)
write.csv(pathway_up_all_LUSC,"pathway_table_up_all_LUSC_tumorPurity.csv",row.names = FALSE)
#down
down_all_LUSC <- read.table("../../new_down_all_LUSC.txt")
down_all_LUSC_genes <- down_all_LUSC$V1
subset <- background_all_LUSC[background_all_LUSC$geneName %in% down_all_LUSC_genes, ]
pathway_down_all_LUSC <- enrich_pathway(background_all_LUSC,down_all_LUSC_genes,"pathway_down_all_LUSC_tumorPurity",TRUE,0.05)
write.csv(pathway_down_all_LUSC,"pathway_table_down_all_LUSC_tumorPurity.csv",row.names = FALSE)

#----------------------------------------------------------------------------------
# LUAD paired
#----------------------------------------------------------------------------------

dataframe_paired_LUAD <- get(load("../../../../1-download_preprocessing/LUAD/paired/LUAD_PreprocessedData_paired_tumorPurity.rda"))
SE_LUAD <- get(load("../../../../1-download_preprocessing/LUAD/paired/LUAD_Illumina_HiSeq_paired_tumorPurity.rda"))
background_paired_LUAD <- subset(rowData(SE_LUAD), rowData(SE_LUAD)$gene_id %in% rownames(dataframe_paired_LUAD), select=c("gene_id","entrezgene"))
colnames(background_paired_LUAD) <- c("geneName","entrez")
#up
up_paired_LUAD <- read.table("../../new_up_paired_LUAD.txt")
up_paired_LUAD_genes <- up_paired_LUAD$V1
pathway_up_paired_LUAD <- enrich_pathway(background_paired_LUAD,up_paired_LUAD_genes,"pathway_up_paired_LUAD_tumorPurity",TRUE,0.05)
write.csv(pathway_up_paired_LUAD, "pathway_table_up_paired_LUAD_tumorPurity.csv",row.names=FALSE)

#down
down_paired_LUAD <- read.table("../../new_down_paired_LUAD.txt")
down_paired_LUAD_genes <- down_paired_LUAD$V1
subset <- background_paired_LUAD[background_paired_LUAD$geneName %in% down_paired_LUAD_genes, ]
pathway_down_paired_LUAD <- enrich_pathway(background_paired_LUAD,down_paired_LUAD_genes,"pathway_down_paired_LUAD_tumorPurity",TRUE,0.05)
write.csv(pathway_down_paired_LUAD, "pathway_table_down_paired_LUAD_tumorPurity.csv",row.names = FALSE)


#--------------------------------------------------------------------------------
# LUSC paired after removing low tumor purity samples
#-------------------------------------------------------------------------------
dataframe_paired_LUSC <- get(load("../../../../1-download_preprocessing/LUSC/paired/LUSC_PreprocessedData_paired_tumorPurity.rda"))
SE_LUSC <- get(load("../../../../1-download_preprocessing/LUSC/paired/LUSC_Illumina_HiSeq_paired_tumorPurity.rda"))
background_paired_LUSC <- subset(rowData(SE_LUSC), rowData(SE_LUSC)$gene_id %in% rownames(dataframe_paired_LUSC), select=c("gene_id","entrezgene"))
colnames(background_paired_LUSC) <- c("geneName","entrez")
#up
up_paired_LUSC <- read.table("../../new_up_paired_LUSC.txt")
up_paired_LUSC_genes <- up_paired_LUSC$V1
subset <- background_paired_LUSC[background_paired_LUSC$geneName %in% up_paired_LUSC_genes, ]
pathway_up_paired_LUSC <- enrich_pathway(background_paired_LUSC,up_paired_LUSC_genes,"pathway_up_paired_LUSC_tumorPurity_test",TRUE,0.05)
write.csv(pathway_up_paired_LUSC,"pathway_table_up_paired_LUSC_tumorPurity_test.csv",row.names = FALSE)
#down
down_paired_LUSC <- read.table("../../new_down_paired_LUSC.txt")
down_paired_LUSC_genes <- down_paired_LUSC$V1
pathway_down_paired_LUSC <- enrich_pathway(background_paired_LUSC,down_paired_LUSC_genes,"pathway_down_paired_LUSC_tumorPurity",TRUE,0.05)
write.csv(pathway_down_paired_LUSC,"pathway_table_down_paired_LUSC_tumorPurity.csv",row.names = FALSE)

