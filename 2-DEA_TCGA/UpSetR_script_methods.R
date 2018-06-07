
library(UpSetR)


#################################################################################

# LUAD

######### up all  ###########

edgeR <- read.table("./LUAD/all/up_edgeR_LUAD_all_tss_tumorPurity.txt")
edgeRTCGA <- read.table("./LUAD/all/up_edgeR_TCGA_LUAD_all_tumorPurity.txt")
limma<- read.table("./LUAD/all/up_limma_LUAD_all_tss_tumorPurity.txt")

listInput <- list(edgeR$V1,edgeRTCGA$V1,limma$V1)
names(listInput) <- c("edgeR","edgeR-TCGAb","limma")

pdf("./plots_DEA_comparison_methods/overlapMethods_LUAD_all_up.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.8,3,3),
      mainbar.y.label = "Intersection Size LUAD-all", sets.x.label = "DEGs per method",
      line.size = 1)
dev.off()

######### up unpaired ##########

edgeR <- read.table("./LUAD/unpaired/up_edgeR_LUAD_unpaired_tss_tumorPurity.txt")
edgeRTCGA <- read.table("./LUAD/unpaired/up_edgeR_TCGA_LUAD_unpaired_tumorPurity.txt")
limma<- read.table("./LUAD/unpaired/up_limma_LUAD_unpaired_tss_tumorPurity.txt")

listInput <- list(edgeR$V1,edgeRTCGA$V1,limma$V1)
names(listInput) <- c("edgeR","edgeR-TCGAb","limma")

pdf("./plots_DEA_comparison_methods/overlapMethods_LUAD_unpaired_up.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.8,3,3),
      mainbar.y.label = "Intersection Size LUAD-unpaired", sets.x.label = "DEGs per method",
      line.size = 1)
dev.off()

######### up paired #############

edgeR <- read.table("./LUAD/paired/up_edgeR_LUAD_paired_tumorPurity.txt")
edgeRTCGA <- read.table("./LUAD/paired/up_edgeR_TCGA_LUAD_paired_tumorPurity.txt")
limma<- read.table("./LUAD/paired/up_limma_LUAD_paired_tumorPurity.txt")

listInput <- list(edgeR$V1,edgeRTCGA$V1,limma$V1)
names(listInput) <- c("edgeR","edgeR-TCGAb","limma")

pdf("./plots_DEA_comparison_methods/overlapMethods_LUAD_paired_up.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,2,3,3),
      mainbar.y.label = "Intersection Size LUAD-paired", sets.x.label = "DEGs per method",
      line.size = 1)
dev.off()


##############################################################################
# LUSC

# down all

edgeR <- read.table("./LUSC/all/down_edgeR_LUSC_all_tss_tumorPurity.txt")
edgeRTCGA <- read.table("./LUSC/all/down_edgeR_TCGA_LUSC_all_tumorPurity.txt")
limma <- read.table("./LUSC/all/down_limma_LUSC_all_tss_tumorPurity.txt")

listInput <- list(edgeR$V1,edgeRTCGA$V1,limma$V1)
names(listInput) <- c("edgeR","edgeR-TCGAb","limma")

pdf("./plots_DEA_comparison_methods/overlapMethods_LUSC_all_down.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.7,3,3),
      mainbar.y.label = "Intersection Size LUSC-all", sets.x.label = "DEGs per method",
      line.size = 1)
dev.off()

# down unpaired

edgeR <- read.table("./LUSC/unpaired/down_edgeR_LUSC_unpaired_tss_tumorPurity.txt")
edgeRTCGA <- read.table("./LUSC/unpaired/down_edgeR_TCGA_LUSC_unpaired_tumorPurity.txt")
limma <- read.table("./LUSC/unpaired/down_limma_LUSC_unpaired_tss_tumorPurity.txt")

listInput <- list(edgeR$V1,edgeRTCGA$V1,limma$V1)
names(listInput) <- c("edgeR","edgeR-TCGAb","limma")

pdf("./plots_DEA_comparison_methods/overlapMethods_LUSC_unpaired_down.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.7,3,3),
      mainbar.y.label = "Intersection Size LUSC-unpaired", sets.x.label = "DEGs per method",
      line.size = 1)
dev.off()


# down paired

edgeR <- read.table("./LUSC/paired/down_edgeR_LUSC_paired_tumorPurity.txt")
edgeRTCGA <- read.table("./LUSC/paired/down_edgeR_TCGA_LUSC_paired_tumorPurity.txt")
limma <- read.table("./LUSC/paired/down_limma_LUSC_paired_tumorPurity.txt")

listInput <- list(edgeR$V1,edgeRTCGA$V1,limma$V1)
names(listInput) <- c("edgeR","edgeR-TCGAb","limma")

pdf("./plots_DEA_comparison_methods/overlapMethods_LUSC_paired_down.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,2,3,3),
      mainbar.y.label = "Intersection Size LUSC-paired", sets.x.label = "DEGs per method",
      line.size = 1)
dev.off()










