source("GOplot_function.R")
source("GeneOntology_functions.R")
load("HUGO.RData")

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

#-------------------------------------------------------------------------------------------
# LUAD paired
#----------------------------------------------------------------------------------------
#up
load("../../../../LUAD/paired/LUAD_PreprocessedData_paired_tumorPurity.rda")
my.gene.univers <- rownames(dataFilt)
up_paired_LUAD <- read.table("../new_up_paired_LUAD.txt")
my.up.genes <- up_paired_LUAD$V1

logFC_file <- read.csv("../../../../LUAD/paired/limma/limma_LUAD_paired_tumorPurity.csv")
plot <- "circle_BP_up_paired_LUAD.png"
GOanalysis_plot("BP",my.gene.univers,my.up.genes,HUGO,20,
                "GO.up.paired.LUAD",logFC_file,plot)
dev.off()

#down
down_paired_LUAD <- read.table("../new_down_paired_LUAD.txt")
my.down.genes <- down_paired_LUAD$V1
logFC_file <- read.csv("../../../../LUAD/paired/limma/limma_LUAD_paired_tumorPurity.csv")
plot <- "circle_BP_down_paired_LUAD.png"
GOanalysis_plot("BP",my.gene.univers,my.down.genes,HUGO,20,
                "GO.down.paired.LUAD",logFC_file,plot)
dev.off()


#------------------------------------------------------------------------------------
#LUAD all
#-----------------------------------------------------------------------------------
#up
load("../../../../LUAD/all/LUAD_PreprocessedData_all_tumorPurity.rda")
my.gene.univers <- rownames(dataFilt)
up_all_LUAD <- read.table("../new_up_all_LUAD.txt")
my.up.genes <- up_all_LUAD$V1

logFC_file <- read.csv("../../../../LUAD/all/limma/limma_LUAD_all_tss_tumorPurity.csv")
plot <- "circle_BP_up_all_LUAD.png"
GOanalysis_plot("BP",my.gene.univers,my.up.genes,HUGO,20,
                "GO.up.all.LUAD",logFC_file,plot)
dev.off()


#down
load("../../../../LUAD/all/LUAD_PreprocessedData_all_tumorPurity.rda")
my.gene.univers <- rownames(dataFilt)
down_all_LUAD <- read.table("../new_down_all_LUAD.txt")
my.down.genes <- down_all_LUAD$V1

logFC_file <- read.csv("../../../../LUAD/all/limma/limma_LUAD_all_tss_tumorPurity.csv")
plot <- "circle_BP_down_all_LUAD.png"
GOanalysis_plot("BP",my.gene.univers,my.down.genes,HUGO,20,
                "GO.down.all.LUAD",logFC_file,plot)
dev.off()


#--------------------------------------------------------------------------------
# LUSC all
#--------------------------------------------------------------------------------
#up
load("../../../../LUSC/all/LUSC_PreprocessedData_all_tumorPurity.rda")
my.gene.univers <- rownames(dataFilt)
up_all_LUSC <- read.table("../new_up_all_LUSC.txt")
my.up.genes <- up_all_LUSC$V1
logFC_file <- read.csv("../../../../LUSC/all/limma/limma_LUSC_all_tss_tumorPurity.csv")
plot <- "circle_BP_up_all_LUSC.png"
GOanalysis_plot("BP",my.gene.univers,my.up.genes,HUGO,10,
                "GO.up.all.LUSC",logFC_file,plot)
dev.off()

#down
load("../../../../LUSC/all/LUSC_PreprocessedData_all_tumorPurity.rda")
my.gene.univers <- rownames(dataFilt)
down_all_LUSC <- read.table("../new_down_all_LUSC.txt")
my.down.genes <- down_all_LUSC$V1
logFC_file <- read.csv("../../../../LUSC/all/limma/limma_LUSC_all_tss_tumorPurity.csv")
plot <- "circle_BP_down_all_LUSC.png"
GOanalysis_plot("BP",my.gene.univers,my.down.genes,HUGO,7,
                "GO.down.all.LUSC",logFC_file,plot)
dev.off()

#---------------------------------------------------------------------------------
#LUSC paired
#---------------------------------------------------------------------------------
#up
load("../../../../LUSC/paired/LUSC_PreprocessedData_paired_tumorPurity.rda")
my.gene.univers <- rownames(dataFilt)
up_paired_LUSC <- read.table("../new_up_paired_LUSC.txt")
my.up.genes <- up_paired_LUSC$V1
logFC_file <- read.csv("../../../../LUSC/paired/limma/limma_LUSC_paired_tumorPurity.csv")
plot <- "circle_BP_up_paired_LUSC.png"
GOanalysis_plot("BP",my.gene.univers,my.up.genes,HUGO,10,
                "GO.up.paired.LUSC",logFC_file,plot)
dev.off()

#down
load("../../../../LUSC/paired/LUSC_PreprocessedData_paired_tumorPurity.rda")
my.gene.univers <- rownames(dataFilt)
down_paired_LUSC <- read.table("../new_down_paired_LUSC.txt")
my.down.genes <- down_paired_LUSC$V1
logFC_file <- read.csv("../../../../LUSC/paired/limma/limma_LUSC_paired_tumorPurity.csv")
plot <- "circle_BP_down_paired_LUSC.png"
GOanalysis_plot("BP",my.gene.univers,my.down.genes,HUGO,10,
                "GO.down.paired.LUSC",logFC_file,plot)
dev.off()
