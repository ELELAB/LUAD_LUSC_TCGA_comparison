#this script needs to be run after performing the DEA with the three pipelines 
#the user also will need to work with the same structure of directories used in this repository to be able to run the scripts as they are
#setwd("~/2-DEA_TCGA")
#the purpose is to calculate the overlap between the DE genes from the three methods applied to the same curation of the dataset

library(VennDiagram)

#---------------------------------------------------------------------------------

# VennDiagram between genes detected from three methods using the same dataset

#---------------------------------------------------------------------------------

Venn_Diagram_3methods <- function(edgeR,edgeRTCGA,limma,direction,cancerType,dataset){
  
  edgeR <- read.table(edgeR,col.names = "genes")
  edgeRTCGA <- read.table(edgeRTCGA,col.names = "genes")
  limma <- read.table(limma,col.names = "genes")
  
  png(paste0(dataset,"/vennDiagram_edgeR-edgeRTCGA-limma_",direction,"_",dataset,"_",cancerType,".png"), height = 800, width = 1100)
  venn <- venn.diagram(list(edgeR$genes, limma$genes, edgeRTCGA$genes), 
                       category.names = c("edgeR","limma","edgeRTCGA"),lwd = 0.5,
                       height = 5000, width = 4500, col="grey",
                       filename = NULL,
                       main = (paste0(dataset,"-",cancerType)),
                       fill=c(4,2,3),alpha=0.3,main.cex=3.5,cat.cex=3, cex=3)
  grid.draw(venn)
  dev.off()
  
  overlap_edgeR_edgeRTCGA <- intersect(edgeR$genes,edgeRTCGA$genes)
  overlap3 <- intersect(overlap_edgeR_edgeRTCGA,limma$genes)
  write.table(overlap3,paste0(dataset,"/overlap_edgeR-edgeRTCGA-limma_",direction,"_",dataset,"_",cancerType,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  overlap_edgeR_limma <- intersect(edgeR$genes,limma$genes)
  overlap_limma_edgeRTCGA <- intersect(limma$genes,edgeRTCGA$genes)
  
  no_overlap.edgeR <- edgeR$genes[!(edgeR$genes %in% overlap_edgeR_edgeRTCGA)]
  no_overlap.edgeR <- no_overlap.edgeR[!(no_overlap.edgeR %in% overlap_edgeR_limma)]
  write.table(no_overlap.edgeR,paste0(dataset,"/no_overlap_",direction,"_",dataset,"_edgeR_",cancerType,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  no_overlap.edgeRTCGA <- edgeRTCGA$genes[!(edgeRTCGA$genes %in% overlap_limma_edgeRTCGA)]
  no_overlap.edgeRTCGA <- no_overlap.edgeRTCGA[!(no_overlap.edgeRTCGA %in% overlap_edgeR_edgeRTCGA)]
  write.table(no_overlap.edgeRTCGA,paste0(dataset,"/no_overlap_",direction,"_",dataset,"_TCGAedgeR_",cancerType,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  no_overlap.limma <- limma$genes[!(limma$genes %in% overlap_edgeR_limma)]
  no_overlap.limma <- no_overlap.limma[!(no_overlap.limma %in% overlap_limma_edgeRTCGA)]
  write.table(no_overlap.limma,paste0(dataset,"/no_overlap_",direction,"_",dataset,"_limma_",cancerType,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  return ("Venn Diagram done!")
}


#---------------------------------------------------------------------------------
# LUAD paired
#--------------------------------------------------------------------------------
#up
edgeR <- "./LUAD/paired/up_edgeR_LUAD_paired_tumorPurity.txt"
edgeRTCGA <- "./LUAD/paired/up_edgeR_TCGA_LUAD_paired_tumorPurity.txt"
limma<- "./LUAD/paired/up_limma_LUAD_paired_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"up_paired","LUAD","plots_DEA_comparison_methods")


#down
edgeR <- "./LUAD/paired/down_edgeR_LUAD_paired_tumorPurity.txt"
edgeRTCGA <- "./LUAD/paired/down_edgeR_TCGA_LUAD_paired_tumorPurity.txt"
limma<- "./LUAD/paired/down_limma_LUAD_paired_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"down_paired","LUAD","plots_DEA_comparison_methods")


#--------------------------------------------------------------------------------
# LUSC paired
#---------------------------------------------------------------------------------
#up
edgeR <- "./LUSC/paired/up_edgeR_LUSC_paired_tumorPurity.txt"
edgeRTCGA <- "./LUSC/paired/up_edgeR_TCGA_LUSC_paired_tumorPurity.txt"
limma <- "./LUSC/paired/up_limma_LUSC_paired_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"up_paired","LUSC","plots_DEA_comparison_methods")

#down
edgeR <- "./LUSC/paired/down_edgeR_LUSC_paired_tumorPurity.txt"
edgeRTCGA <- "./LUSC/paired/down_edgeR_TCGA_LUSC_paired_tumorPurity.txt"
limma <- "./LUSC/paired/down_limma_LUSC_paired_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"down_paired","LUSC","plots_DEA_comparison_methods")


#--------------------------------------------------------------------------------
# LUAD all
#-------------------------------------------------------------------------------
#up
edgeR <- "./LUAD/all/up_edgeR_LUAD_all_tss_tumorPurity.txt"
edgeRTCGA <- "./LUAD/all/up_edgeR_TCGA_LUAD_all_tumorPurity.txt"
limma<- "./LUAD/all/up_limma_LUAD_all_tss_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"up_all","LUAD","plots_DEA_comparison_methods")

#down
edgeR <- "./LUAD/all/down_edgeR_LUAD_all_tss_tumorPurity.txt"
edgeRTCGA <- "./LUAD/all/down_edgeR_TCGA_LUAD_all_tumorPurity.txt"
limma<- "./LUAD/all/down_limma_LUAD_all_tss_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"down_all","LUAD","plots_DEA_comparison_methods")


#-----------------------------------------------------------------------------------
# LUAD unpaired
#---------------------------------------------------------------------------------
#up
edgeR <- "./LUAD/unpaired/up_edgeR_LUAD_unpaired_tss_tumorPurity.txt"
edgeRTCGA <- "./LUAD/unpaired/up_edgeR_TCGA_LUAD_unpaired_tumorPurity.txt"
limma<- "./LUAD/unpaired/up_limma_LUAD_unpaired_tss_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"up_unpaired","LUAD","plots_DEA_comparison_methods")

#down
edgeR <- "./LUAD/unpaired/down_edgeR_LUAD_unpaired_tss_tumorPurity.txt"
edgeRTCGA <- "./LUAD/unpaired/down_edgeR_TCGA_LUAD_unpaired_tumorPurity.txt"
limma<- "./LUAD/unpaired/down_limma_LUAD_unpaired_tss_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"down(unpaired","LUAD","plots_DEA_comparison_methods")

#--------------------------------------------------------------------------------
#LUSC all
#-------------------------------------------------------------------------------

#up
edgeR <- "./LUSC/all/up_edgeR_LUSC_all_tss_tumorPurity.txt"
edgeRTCGA <- "./LUSC/all/up_edgeR_TCGA_LUSC_all_tumorPurity.txt"
limma <- "./LUSC/all/up_limma_LUSC_all_tss_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"up_all","LUSC","plots_DEA_comparison_methods")

#down
edgeR <- "./LUSC/all/down_edgeR_LUSC_all_tss_tumorPurity.txt"
edgeRTCGA <- "./LUSC/all/down_edgeR_TCGA_LUSC_all_tumorPurity.txt"
limma <- "./LUSC/all/down_limma_LUSC_all_tss_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"down_all","LUSC","plots_DEA_comparison_methods")


#--------------------------------------------------------------------------------
# LUSC unpaired
#---------------------------------------------------------------------------------
#up
edgeR <- "./LUSC/unpaired/up_edgeR_LUSC_unpaired_tss_tumorPurity.txt"
edgeRTCGA <- "./LUSC/unpaired/up_edgeR_TCGA_LUSC_unpaired_tumorPurity.txt"
limma <- "./LUSC/unpaired/up_limma_LUSC_unpaired_tss_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"up_unpaired","LUSC","plots_DEA_comparison_methods")

#down
edgeR <- "./LUSC/unpaired/down_edgeR_LUSC_unpaired_tss_tumorPurity.txt"
edgeRTCGA <- "./LUSC/unpaired/down_edgeR_TCGA_LUSC_unpaired_tumorPurity.txt"
limma <- "./LUSC/unpaired/down_limma_LUSC_unpaired_tss_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"down_unpaired","LUSC","plots_DEA_comparison_methods")
