#this script needs to be run after performing the DEA with the three pipelines and the different curations of datasets
#the user also will need to work with the same structure of directories used in this repository to be able to run the scripts as they are
#setwd("~/2-DEA_TCGA")
#the purpose is to calculate the overlap between the DE genes from each curation of the dataset


#---------------------------------------------------------------------------------------------

# Function to compare the DE genes detected by the same method using three different dataset

# pairedFile, allFile, unpairedFile: path where the DEGs are stored for paired, all, unpaired dataset, respectively
# method: "edgeR", "edgeR_TCGA", "limma"
# direction: "up" or "down"
# cancer_type: "LUAD" or "LUSC" 
#--------------------------------------------------------------------------------------------

library(VennDiagram)

Venn_Diagram_dataset <- function(pairedFile,allFile,unpairedFile,method,direction,cancerType){
  
  paired <- read.table(pairedFile,col.names = "genes")
  all <- read.table(allFile,col.names = "genes")
  unpaired <- read.table(unpairedFile,col.names = "genes")
  
  #pdf(paste0(dir,method,"/vennDiagram_",method,"_",direction,"_",cancerType,".pdf"))
  png(paste0(method,"/vennDiagram_",method,"_",direction,"_",cancerType,".png"), height = 700, width = 900)
  venn <- venn.diagram(list(paired$genes, all$genes, unpaired$genes), 
                       #category.names = c(paste0(cancerType,"_paired"),paste0(cancerType,"_all"),paste0(cancerType,"_unpaired")),
                       category.names = c("paired","all","unpaired"),
                       #lwd = 0,
                       filename = NULL, col= "grey",
                       #height = 5000, width = 4500,
                       main = (paste0(method,"-",cancerType)),
                       fill=c(4,2,3),alpha=0.3,main.cex=3.5,sub.cex = 1, cat.cex=3, cex=3)
  grid.draw(venn)
  dev.off()
  
  overlap_paired_all <- intersect(paired$genes,all$genes)
  overlap3 <- intersect(overlap_paired_all,unpaired$genes)
  write.table(overlap3,paste0(method,"/overlap_",direction,"_",method,"_",cancerType,".txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  overlap_paired_unpaired <- intersect(paired$genes,unpaired$genes)
  overlap_all_unpaired <- intersect(all$genes,unpaired$genes)
  
  no_overlap.paired <- paired$genes[!(paired$genes %in% overlap_paired_all)]
  no_overlap.paired <- no_overlap.paired[!(no_overlap.paired %in% overlap_paired_unpaired)]
  write.table(no_overlap.paired,paste0(method,"/no_overlap_",method,"_",direction,"_",cancerType,"_","pairedOnly.txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  no_overlap.all <- all$genes[!(all$genes %in% overlap_paired_all)]
  no_overlap.all <- no_overlap.all[!(no_overlap.all %in% overlap_all_unpaired)]
  write.table(no_overlap.all,paste0(method,"/no_overlap_",method,"_",direction,"_",cancerType,"_","allOnly.txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  no_overlap.unpaired <- unpaired$genes[!(unpaired$genes %in% overlap_all_unpaired)]
  no_overlap.unpaired <- no_overlap.unpaired[!(no_overlap.unpaired %in% overlap_paired_unpaired)]
  write.table(no_overlap.unpaired,paste0(method,"/no_overlap_",method,"_",direction,"_",cancerType,"_","unpairedOnly.txt"),sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  
  return(print("Venn Diagram done!!"))
  
}

#--------------------------------------------------------------------------------
# LUAD edgeR
#-------------------------------------------------------------------------------
#up
paired <- "./LUAD/paired/up_edgeR_LUAD_paired_tumorPurity.txt"
all <- "./LUAD/all/up_edgeR_LUAD_all_tss_tumorPurity.txt"
unpaired <- "./LUAD/unpaired/up_edgeR_LUAD_unpaired_tss_tumorPurity.txt"
Venn_Diagram_dataset(paired,all,unpaired,"plots_DEA_comparison_dataset","up_edgeR","LUAD")

#down
paired <- "./LUAD/paired/down_edgeR_LUAD_paired_tumorPurity.txt"
all <- "./LUAD/all/down_edgeR_LUAD_all_tss_tumorPurity.txt"
unpaired <- "./LUAD/unpaired/down_edgeR_LUAD_unpaired_tss_tumorPurity.txt"
Venn_Diagram_dataset(paired,all,unpaired,"plots_DEA_comparison_dataset","down_edgeR","LUAD")

#----------------------------------------------------------------------------------
# LUAD edgeR_TCGA
#--------------------------------------------------------------------------------
#up
paired <- "./LUAD/paired/up_edgeR_TCGA_LUAD_paired_tumorPurity.txt"
all <- "./LUAD/all/up_edgeR_TCGA_LUAD_all_tumorPurity.txt"
unpaired <- "./LUAD/unpaired/up_edgeR_TCGA_LUAD_unpaired_tumorPurity.txt"
Venn_Diagram_dataset(paired,all,unpaired,"plots_DEA_comparison_dataset","edgeR_TCGA_up","LUAD")

#down
paired <- "./LUAD/paired/down_edgeR_TCGA_LUAD_paired_tumorPurity.txt"
all <- "./LUAD/all/down_edgeR_TCGA_LUAD_all_tumorPurity.txt"
unpaired <- "./LUAD/unpaired/down_edgeR_TCGA_LUAD_unpaired_tumorPurity.txt"
Venn_Diagram_dataset(paired,all,unpaired,"plots_DEA_comparison_dataset","down_edgeR_TCGA","LUAD")

#---------------------------------------------------------------------------------
# LUAD limma
#--------------------------------------------------------------------------------
#up
paired <- "./LUAD/paired/up_limma_LUAD_paired_tumorPurity.txt"
all <- "./LUAD/all/up_limma_LUAD_all_tss_tumorPurity.txt"
unpaired <- "./LUAD/unpaired/up_limma_LUAD_unpaired_tss_tumorPurity.txt"
Venn_Diagram_dataset(paired,all,unpaired,"plots_DEA_comparison_dataset","up_limma","LUAD")

#down
paired <- "./LUAD/paired/down_limma_LUAD_paired_tumorPurity.txt"
all <- "./LUAD/all/down_limma_LUAD_all_tss_tumorPurity.txt"
unpaired <- "./LUAD/unpaired/down_limma_LUAD_unpaired_tss_tumorPurity.txt"
Venn_Diagram_dataset(paired,all,unpaired,"plots_DEA_comparison_dataset","down_limma","LUAD")

#-------------------------------------------------------------------------------
# LUSC edgeR
#-------------------------------------------------------------------------------
#up
paired <- "./LUSC/paired/up_edgeR_LUSC_paired_tumorPurity.txt"
all <- "./LUSC/all/up_edgeR_LUSC_all_tss_tumorPurity.txt"
unpaired <- "./LUSC/unpaired/up_edgeR_LUSC_unpaired_tss_tumorPurity.txt"
Venn_Diagram_dataset(paired,all,unpaired,"plots_DEA_comparison_dataset","up_limma","LUSC")

#down
paired <- "./LUSC/paired/down_edgeR_LUSC_paired_tumorPurity.txt"
all <- "./LUSC/all/down_edgeR_LUSC_all_tss_tumorPurity.txt"
unpaired <- "./LUSC/unpaired/down_edgeR_LUSC_unpaired_tss_tumorPurity.txt"
Venn_Diagram_dataset(paired,all,unpaired,"plots_DEA_comparison_dataset","down_edgeR","LUSC")

#-------------------------------------------------------------------------------
# LUSC edgeR_TCGA
#-------------------------------------------------------------------------------
#up
paired <- "./LUSC/paired/up_edgeR_TCGA_LUSC_paired_tumorPurity.txt"
all <- "./LUSC/all/up_edgeR_TCGA_LUSC_all_tumorPurity.txt"
unpaired <- "./LUSC/unpaired/up_edgeR_TCGA_LUSC_unpaired_tumorPurity.txt"
Venn_Diagram_dataset(paired,all,unpaired,"plots_DEA_comparison_dataset","up_edgeR_TCGA","LUSC")

#down
paired <- "./LUSC/paired/down_edgeR_TCGA_LUSC_paired_tumorPurity.txt"
all <- "./LUSC/all/down_edgeR_TCGA_LUSC_all_tumorPurity.txt"
unpaired <- "./LUSC/unpaired/down_edgeR_TCGA_LUSC_unpaired_tumorPurity.txt"
Venn_Diagram_dataset(paired,all,unpaired,"plots_DEA_comparison_dataset","down_edgeR_TCGA","LUSC")

#----------------------------------------------------------------------
# LUSC limma
#----------------------------------------------------------------------
#up
paired <- "./LUSC/paired/up_limma_LUSC_paired_tumorPurity.txt"
all <- "./LUSC/all/up_limma_LUSC_all_tss_tumorPurity.txt"
unpaired <- "./LUSC/unpaired/up_limma_LUSC_unpaired_tss_tumorPurity.txt"
Venn_Diagram_dataset(paired,all,unpaired,"plots_DEA_comparison_dataset","up_limma","LUSC")

#down
paired <- "./LUSC/paired/down_limma_LUSC_paired_tumorPurity.txt"
all <- "./LUSC/all/down_limma_LUSC_all_tss_tumorPurity.txt"
unpaired <- "./LUSC/unpaired/down_limma_LUSC_unpaired_tss_tumorPurity.txt"
Venn_Diagram_dataset(paired,all,unpaired,"plots_DEA_comparison_dataset","down_limma","LUSC")
