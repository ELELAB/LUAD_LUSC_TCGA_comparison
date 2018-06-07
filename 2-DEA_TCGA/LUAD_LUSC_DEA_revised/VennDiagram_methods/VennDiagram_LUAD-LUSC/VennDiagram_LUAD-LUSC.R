#-------------------------------------------------------------------

# VennDiagram that compare LUAD-LUSC

#-------------------------------------------------------------------
library(VennDiagram)

VennDiagram_LUAD_LUSC <- function(listLUAD,listLUSC,direction,dataset,dir){
  
  listLUAD<- read.table(listLUAD,col.names = "genes")
  listLUSC <- read.table(listLUSC,col.names = "genes")
  
  png(paste0(dir,"vennDiagram_",direction,"_",dataset,"_LUSC-LUAD.png"), height = 800, width = 1000)
  venn <- venn.diagram(list(listLUAD$genes, listLUSC$genes), 
                       category.names = c("LUAD","LUSC"), lwd = 0.5,
                       #filename =paste0("/data/user/marta/pipeline/DE/VennDiagram_limma_edgeR_edgeRTCGA/VennDiagram_LUSC-LUAD/vennDiagram_",direction,"_",dataset,"_LUSC-LUAD.pdf"), lwd = 0.5,
                       #height = 6700, width = 10000,
                       filename = NULL, main = (paste0(dataset)),
                       fill=c(4,2),alpha=0.3,main.cex=4, sub.cex = 1, cat.cex=3, cex=4)
  grid.draw(venn)
  
  dev.off()
  
  overlap <- intersect(listLUAD$genes,listLUSC$genes)
  write.table(overlap,paste0(dir,"overlap_",direction,"_",dataset,"_LUSC-LUAD.txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  no_overlap.LUAD <- listLUAD$genes[!(listLUAD$genes %in% overlap)]
  write.table(no_overlap.LUAD,paste0(dir,"no_overlap_",direction,"_",dataset,"_LUAD_only.txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  no_overlap.LUSC <- listLUSC$genes[!(listLUSC$genes %in% overlap)]
  write.table(no_overlap.LUSC,paste0(dir,"no_overlap_",direction,"_",dataset,"_LUSC_only.txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  return("Venn Diagram done!")
}


#------------------------------------------------------------------------------
# up paired after removing low tumor purity samples
#--------------------------------------------------------------------------------

up_LUAD <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/paired/overlap_edgeR-edgeRTCGA-limma_up_paired_LUAD.txt"
up_LUSC <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/paired/overlap_edgeR-edgeRTCGA-limma_up_paired_LUSC.txt"
dir <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/VennDiagram_LUAD-LUSC/"
VennDiagram_LUAD_LUSC(up_LUAD,up_LUSC,"up","paired",dir)

#----------------------------------------------------------------------------------
# down paired after removing low tumor purity samples
#---------------------------------------------------------------------------------

down_LUAD <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/paired/overlap_edgeR-edgeRTCGA-limma_down_paired_LUAD.txt"
down_LUSC <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/paired/overlap_edgeR-edgeRTCGA-limma_down_paired_LUSC.txt"
dir <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/VennDiagram_LUAD-LUSC/"
VennDiagram_LUAD_LUSC(down_LUAD,down_LUSC,"down","paired",dir)


#----------------------------------------------------------------------------------
#up all after removing low tumor purity
#----------------------------------------------------------------------------------
up_LUAD <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/all/overlap_edgeR-edgeRTCGA-limma_up_all_LUAD.txt"
up_LUSC <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/all/overlap_edgeR-edgeRTCGA-limma_up_all_LUSC.txt"
dir <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/VennDiagram_LUAD-LUSC/"
VennDiagram_LUAD_LUSC(up_LUAD,up_LUSC,"up","all",dir)

#----------------------------------------------------------------------------------
#down all after removing low tumor purity
#----------------------------------------------------------------------------------
down_LUAD <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/all/overlap_edgeR-edgeRTCGA-limma_down_all_LUAD.txt"
down_LUSC <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/all/overlap_edgeR-edgeRTCGA-limma_down_all_LUSC.txt"
dir <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/VennDiagram_LUAD-LUSC/"
VennDiagram_LUAD_LUSC(down_LUAD,down_LUSC,"down","all",dir)
