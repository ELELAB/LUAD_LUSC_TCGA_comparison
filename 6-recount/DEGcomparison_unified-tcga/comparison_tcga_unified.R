#setwd("~/6-recount/DEGcomparison_unified-tcga/")
library(VennDiagram)
#type1=tcga
#type2=unified
#outputname=specify cancertype
VennDiagram_comparison_unified_tcga <- function(list1,list2,type1,type2,cancerType,direction){

list1 <- read.table(list1,col.names = "genes")
list2 <- read.table(list2,col.names = "genes")
pdf(paste0("./vennDiagram_tcga_unified_",direction,"_",cancerType,".pdf"))
#png(paste0(dir,dataset,"/vennDiagram_edgeR-edgeRTCGA-limma_",direction,"_",dataset,"_",cancerType,".png"), height = 800, width = 1100)
venn <- venn.diagram(list(list1$genes, list2$genes), 
                     category.names = c(type1,type2),
                     lwd = 0.1,
                     #filename =paste0(dir,dataset,"/vennDiagram_edgeR-edgeRTCGA-limma_",direction,"_",dataset,"_",cancerType,".pdf"), 
                     height = 5000, width = 4500, col="grey",
                     filename = NULL,
                     main = (paste0(cancerType," ",direction,"-regulated genes")),
                     fill=c(4,2),alpha=0.3,main.cex=2,cat.cex=2, cex=2)
grid.draw(venn)
dev.off()

overlap <- intersect(list1$genes,list2$genes)
write.table(overlap,paste0("overlap_tcga_unified_",direction,"_",cancerType,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)

no_overlap.list1 <- list1$genes[!(list1$genes %in% overlap)]
write.table(no_overlap.list1,paste0("no_overlap_",direction,"_",cancerType,"_",type1,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)

no_overlap.list2 <- list2$genes[!(list2$genes %in% overlap)]
write.table(no_overlap.list2,paste0("no_overlap_",direction,"_",cancerType,"_",type2,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)

return ("Venn Diagram done!")
}

####################

# LUAD up
setwd("~/6-recount/DEGcomparison_unified-tcga/LUAD/")

up_luad_tcga <- "../../../2-DEA_TCGA/LUAD/all/up_limma_LUAD_all_tss_tumorPurity.txt"
up_luad_unified <- "../../LUAD/limma/up_limma_unified_LUAD_tss.txt"
VennDiagram_comparison_unified_tcga(up_luad_tcga,up_luad_unified,"tcga","unified","LUAD","up")


# LUAD down
down_luad_tcga <- "../../../2-DEA_TCGA/LUAD/all/down_limma_LUAD_all_tss_tumorPurity.txt"
down_luad_unified <- "../../LUAD/limma/down_limma_unified_LUAD_tss.txt"
VennDiagram_comparison_unified_tcga(down_luad_tcga,down_luad_unified,"tcga","unified","LUAD","down")


##################
setwd("~/6-recount/DEGcomparison_unified-tcga/LUSC/")
# LUSC up

up_lusc_tcga <- "../../../2-DEA_TCGA/LUSC/all/up_limma_LUSC_all_tss_tumorPurity.txt"
up_lusc_unified <- "../../LUSC/limma/up_limma_unified_LUSC_tss.txt"
VennDiagram_comparison_unified_tcga(up_lusc_tcga,up_lusc_unified,"tcga","unified","LUSC","up")

#down

down_lusc_tcga <- "../../../2-DEA_TCGA/LUSC/all/limma/down_limma_LUSC_all_tss_tumorPurity.txt"
down_lusc_unified <- "../../LUSC/limma/down_limma_unified_LUSC_tss.txt"
VennDiagram_comparison_unified_tcga(down_lusc_tcga,down_lusc_unified,"tcga","unified","LUSC","down")


