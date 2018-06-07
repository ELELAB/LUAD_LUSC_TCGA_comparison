# created by Marta (16/01/18)

VennDiagram_LUAD_LUSC <- function(listLUAD,listLUSC,direction){
  
  listLUAD<- read.table(listLUAD,col.names = "genes")
  listLUSC <- read.table(listLUSC,col.names = "genes")
  
  png(paste0("vennDiagram_",direction,"_LUSC-LUAD_unified.png"), height = 800, width = 1000)
  venn <- venn.diagram(list(listLUAD$genes, listLUSC$genes), 
                       category.names = c("LUAD","LUSC"), lwd = 0.5,
                       #filename =paste0("/data/user/marta/pipeline/DE/VennDiagram_limma_edgeR_edgeRTCGA/VennDiagram_LUSC-LUAD/vennDiagram_",direction,"_",dataset,"_LUSC-LUAD.pdf"), lwd = 0.5,
                       #height = 6700, width = 10000,
                       filename = NULL,
                       fill=c(4,2),alpha=0.3,main.cex=4, sub.cex = 1, cat.cex=3, cex=4)
  grid.draw(venn)
  
  dev.off()
  
  overlap <- intersect(listLUAD$genes,listLUSC$genes)
  write.table(overlap,paste0("overlap_",direction,"_LUSC-LUAD.txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  no_overlap.LUAD <- listLUAD$genes[!(listLUAD$genes %in% overlap)]
  write.table(no_overlap.LUAD,paste0("no_overlap_",direction,"_LUAD_only.txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  no_overlap.LUSC <- listLUSC$genes[!(listLUSC$genes %in% overlap)]
  write.table(no_overlap.LUSC,paste0("no_overlap_",direction,"_LUSC_only.txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  return("Venn Diagram done!")
}



#################
# up-regulated genes

listLUAD <- "../LUAD/limma/up_limma_unified_LUAD_tss.txt"
listLUSC <- "../LUSC/limma/up_limma_unified_LUSC_tss.txt"
VennDiagram_LUAD_LUSC(listLUAD,listLUSC,"up")


# down-regulated genes

listLUAD <- "../LUAD/limma/down_limma_unified_LUAD_tss.txt"
listLUSC <- "../LUSC/limma/down_limma_unified_LUSC_tss.txt"
VennDiagram_LUAD_LUSC(listLUAD,listLUSC,"down")




