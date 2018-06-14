#this script needs to be run after performing the DEA with the three pipelines
#the user also will need to work with the same structure of directories used in this repository to be able to run the scripts as they are
#setwd("~/2-DEA_TCGA")

#-------------------------------------------------------------------

# VennDiagram that compare LUAD-LUSC

#-------------------------------------------------------------------
library(VennDiagram)

VennDiagram_LUAD_LUSC <- function(listLUAD,listLUSC,direction,dataset){
  
  listLUAD<- read.table(listLUAD,col.names = "genes")
  listLUSC <- read.table(listLUSC,col.names = "genes")
  
  png(paste0(dir,"vennDiagram_",direction,"_",dataset,"_LUSC-LUAD.png"), height = 800, width = 1000)
  venn <- venn.diagram(list(listLUAD$genes, listLUSC$genes), 
                       category.names = c("LUAD","LUSC"), lwd = 0.5,
                       filename = NULL, main = (paste0(dataset)),
                       fill=c(4,2),alpha=0.3,main.cex=4, sub.cex = 1, cat.cex=3, cex=4)
  grid.draw(venn)
  
  dev.off()
  
  overlap <- intersect(listLUAD$genes,listLUSC$genes)
  write.table(overlap,paste0("overlap_",direction,"_",dataset,"_LUSC-LUAD.txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  no_overlap.LUAD <- listLUAD$genes[!(listLUAD$genes %in% overlap)]
  write.table(no_overlap.LUAD,paste0("no_overlap_",direction,"_",dataset,"_LUAD_only.txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  no_overlap.LUSC <- listLUSC$genes[!(listLUSC$genes %in% overlap)]
  write.table(no_overlap.LUSC,paste0("no_overlap_",direction,"_",dataset,"_LUSC_only.txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  return("Venn Diagram done!")
}


#------------------------------------------------------------------------------
# up paired 
#--------------------------------------------------------------------------------

up_LUAD <- "./plots_DEA_comparison_methods/overlap_edgeR-edgeRTCGA-limma_up_paired_plots_DEA_comparison_methods_LUAD.txt"
up_LUSC <- "./plots_DEA_comparison_methods/overlap_edgeR-edgeRTCGA-limma_up_paired_plots_DEA_comparison_methods_LUSC.txt"
VennDiagram_LUAD_LUSC(up_LUAD,up_LUSC,"up","paired")

#----------------------------------------------------------------------------------
# down paired
#---------------------------------------------------------------------------------

down_LUAD <- "./plots_DEA_comparison_methods/overlap_edgeR-edgeRTCGA-limma_down_paired_plots_DEA_comparison_methods_LUAD.txt"
down_LUSC <- "./plots_DEA_comparison_methods/overlap_edgeR-edgeRTCGA-limma_down_paired_plots_DEA_comparison_methods_LUSC.txt"
VennDiagram_LUAD_LUSC(down_LUAD,down_LUSC,"down","paired")


#----------------------------------------------------------------------------------
#up all
#----------------------------------------------------------------------------------
up_LUAD <- "./plots_DEA_comparison_methods/overlap_edgeR-edgeRTCGA-limma_up_all_plots_DEA_comparison_methods_LUAD.txt"
up_LUSC <- "./plots_DEA_comparison_methods/overlap_edgeR-edgeRTCGA-limma_up_all_plots_DEA_comparison_methods_LUSC.txt"
VennDiagram_LUAD_LUSC(up_LUAD,up_LUSC,"up","all")

#----------------------------------------------------------------------------------
#down all
#----------------------------------------------------------------------------------
down_LUAD <- "./plots_DEA_comparison_methods/overlap_edgeR-edgeRTCGA-limma_down_all_plots_DEA_comparison_methods_LUAD.txt"
down_LUSC <- "./plots_DEA_comparison_methods/overlap_edgeR-edgeRTCGA-limma_down_all_plots_DEA_comparison_methods_LUSC.txt"
VennDiagram_LUAD_LUSC(down_LUAD,down_LUSC,"down","all")

