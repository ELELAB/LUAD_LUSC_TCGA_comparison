library(VennDiagram)
#---------------------------------------------------------------------------------

# VennDiagram between genes detected from three methods using the same datset

#---------------------------------------------------------------------------------

Venn_Diagram_3methods <- function(edgeR,edgeRTCGA,limma,direction,cancerType,dataset,dir){
  
  edgeR <- read.table(edgeR,col.names = "genes")
  edgeRTCGA <- read.table(edgeRTCGA,col.names = "genes")
  limma <- read.table(limma,col.names = "genes")
  
  #pdf(paste0(dir,dataset,"/vennDiagram_edgeR-edgeRTCGA-limma_",direction,"_",dataset,"_",cancerType,".pdf"))
  png(paste0(dir,dataset,"/vennDiagram_edgeR-edgeRTCGA-limma_",direction,"_",dataset,"_",cancerType,".png"))
  venn <- venn.diagram(list(edgeR$genes, limma$genes, edgeRTCGA$genes), 
                       category.names = c("edgeR","limma","edgeRTCGA"),
                        lwd = 0.5,
                       #filename =paste0(dir,dataset,"/vennDiagram_edgeR-edgeRTCGA-limma_",direction,"_",dataset,"_",cancerType,".pdf"), 
                       height = 5000, width = 4500, col="grey",
                       filename = NULL,
                       main = (paste0(dataset,"-",cancerType)),
                       fill=c(4,2,3),alpha=0.3,main.cex=3,cat.cex=2.5, cex=3)
  grid.draw(venn)
  dev.off()
  
  overlap_edgeR_edgeRTCGA <- intersect(edgeR$genes,edgeRTCGA$genes)
  overlap3 <- intersect(overlap_edgeR_edgeRTCGA,limma$genes)
  write.table(overlap3,paste0(dir,dataset,"/overlap_edgeR-edgeRTCGA-limma_",direction,"_",dataset,"_",cancerType,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  overlap_edgeR_limma <- intersect(edgeR$genes,limma$genes)
  overlap_limma_edgeRTCGA <- intersect(limma$genes,edgeRTCGA$genes)
  
  no_overlap.edgeR <- edgeR$genes[!(edgeR$genes %in% overlap_edgeR_edgeRTCGA)]
  no_overlap.edgeR <- no_overlap.edgeR[!(no_overlap.edgeR %in% overlap_edgeR_limma)]
  write.table(no_overlap.edgeR,paste0(dir,dataset,"/no_overlap_",direction,"_",dataset,"_edgeR_",cancerType,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  no_overlap.edgeRTCGA <- edgeRTCGA$genes[!(edgeRTCGA$genes %in% overlap_limma_edgeRTCGA)]
  no_overlap.edgeRTCGA <- no_overlap.edgeRTCGA[!(no_overlap.edgeRTCGA %in% overlap_edgeR_edgeRTCGA)]
  write.table(no_overlap.edgeRTCGA,paste0(dir,dataset,"/no_overlap_",direction,"_",dataset,"_TCGAedgeR_",cancerType,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  no_overlap.limma <- limma$genes[!(limma$genes %in% overlap_edgeR_limma)]
  no_overlap.limma <- no_overlap.limma[!(no_overlap.limma %in% overlap_limma_edgeRTCGA)]
  write.table(no_overlap.limma,paste0(dir,dataset,"/no_overlap_",direction,"_",dataset,"_limma_",cancerType,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  return ("Venn Diagram done!")
}



#---------------------------------------------------------------------------------
# LUAD paired after removing low tumor purity samples
#--------------------------------------------------------------------------------
#up
edgeR <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/paired/edgeR/up_edgeR_LUAD_paired_tumorPurity.txt"
edgeRTCGA <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/paired/edgeR_TCGA/up_edgeR_TCGA_LUAD_paired_tumorPurity.txt"
limma<- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/paired/limma/up_limma_LUAD_paired_tumorPurity.txt"
dir <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/Venn_Diagram_methods/"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"up","LUAD","paired",dir)


#down
edgeR <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/paired/edgeR/down_edgeR_LUAD_paired_tumorPurity.txt"
edgeRTCGA <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/paired/edgeR_TCGA/down_edgeR_TCGA_LUAD_paired_tumorPurity.txt"
limma<- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/paired/limma/down_limma_LUAD_paired_tumorPurity.txt"
dir <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"down","LUAD","paired",dir)


#--------------------------------------------------------------------------------
# LUSC paired after removing low tumor purity samples
#---------------------------------------------------------------------------------
#up
dir <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/"
edgeR <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/paired/edgeR/up_edgeR_LUSC_paired_tumorPurity.txt"
edgeRTCGA <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/paired/edgeR_TCGA/up_edgeR_TCGA_LUSC_paired_tumorPurity.txt"
limma <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/paired/limma/up_limma_LUSC_paired_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"up","LUSC","paired",dir)

#down
dir <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/"
edgeR <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/paired/edgeR/down_edgeR_LUSC_paired_tumorPurity.txt"
edgeRTCGA <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/paired/edgeR_TCGA/down_edgeR_TCGA_LUSC_paired_tumorPurity.txt"
limma <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/paired/limma/down_limma_LUSC_paired_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"down","LUSC","paired",dir)


#--------------------------------------------------------------------------------
# LUAD all after removing low tumor purity samples
#-------------------------------------------------------------------------------
#up
edgeR <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/all/edgeR/up_edgeR_LUAD_all_tss_tumorPurity.txt"
edgeRTCGA <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/all/edgeR_TCGA/up_edgeR_TCGA_LUAD_all_tumorPurity.txt"
limma<- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/all/limma/up_limma_LUAD_all_tss_tumorPurity.txt"
dir <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"up","LUAD","all",dir)

#down
edgeR <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/all/edgeR/down_edgeR_LUAD_all_tss_tumorPurity.txt"
edgeRTCGA <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/all/edgeR_TCGA/down_edgeR_TCGA_LUAD_all_tumorPurity.txt"
limma<- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/all/limma/down_limma_LUAD_all_tss_tumorPurity.txt"
dir <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"down","LUAD","all",dir)



#-----------------------------------------------------------------------------------
# LUAD unpaired after removing low tumor purity samples
#---------------------------------------------------------------------------------
#up
edgeR <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/unpaired/edgeR/up_edgeR_LUAD_unpaired_tss_tumorPurity.txt"
edgeRTCGA <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/unpaired/edgeR_TCGA/up_edgeR_TCGA_LUAD_unpaired_tumorPurity.txt"
limma<- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/unpaired/limma/up_limma_LUAD_unpaired_tss_tumorPurity.txt"
dir <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"up","LUAD","unpaired",dir)

#down
edgeR <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/unpaired/edgeR/down_edgeR_LUAD_unpaired_tss_tumorPurity.txt"
edgeRTCGA <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/unpaired/edgeR_TCGA/down_edgeR_TCGA_LUAD_unpaired_tumorPurity.txt"
limma<- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUAD/unpaired/limma/down_limma_LUAD_unpaired_tss_tumorPurity.txt"
dir <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/VennDiagram_methods/"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"down","LUAD","unpaired",dir)


#--------------------------------------------------------------------------------
#LUSC all after removing low tumor purity samples
#-------------------------------------------------------------------------------

#up
edgeR <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/all/edgeR/up_edgeR_LUSC_all_tss_tumorPurity.txt"
edgeRTCGA <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/all/edgeR_TCGA/up_edgeR_TCGA_LUSC_all_tumorPurity.txt"
limma <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/all/limma/up_limma_LUSC_all_tss_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"up","LUSC","all",dir)

#down
edgeR <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/all/edgeR/down_edgeR_LUSC_all_tss_tumorPurity.txt"
edgeRTCGA <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/all/edgeR_TCGA/down_edgeR_TCGA_LUSC_all_tumorPurity.txt"
limma <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/all/limma/down_limma_LUSC_all_tss_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"down","LUSC","all",dir)


#--------------------------------------------------------------------------------
# LUSC unpaired after removing low tumor purity samples
#---------------------------------------------------------------------------------
#up
edgeR <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/unpaired/edgeR/up_edgeR_LUSC_unpaired_tss_tumorPurity.txt"
edgeRTCGA <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/unpaired/edgeR_TCGA/up_edgeR_TCGA_LUSC_unpaired_tumorPurity.txt"
limma <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/unpaired/limma/up_limma_LUSC_unpaired_tss_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"up","LUSC","unpaired",dir)

#down
edgeR <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/unpaired/edgeR/down_edgeR_LUSC_unpaired_tss_tumorPurity.txt"
edgeRTCGA <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/unpaired/edgeR_TCGA/down_edgeR_TCGA_LUSC_unpaired_tumorPurity.txt"
limma <- "/data/user/marta/pipeline/DE/LUAD_LUSC_May2018/LUSC/unpaired/limma/down_limma_LUSC_unpaired_tss_tumorPurity.txt"
Venn_Diagram_3methods(edgeR,edgeRTCGA,limma,"down","LUSC","unpaired",dir)
