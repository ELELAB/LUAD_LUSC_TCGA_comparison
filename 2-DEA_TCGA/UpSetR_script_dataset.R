library(UpSetR)


#################################################################################

# LUAD

######### up limma  ###########

paired <- read.table("./LUAD/paired/up_limma_LUAD_paired_tumorPurity.txt")
unpaired <- read.table("./LUAD/unpaired/up_limma_LUAD_unpaired_tss_tumorPurity.txt")
all <- read.table("./LUAD/all/up_limma_LUAD_all_tss_tumorPurity.txt")

listInput <- list(paired$V1,unpaired$V1,all$V1)
names(listInput) <- c("paired","unpaired","all")

pdf("./plots_DEA_comparison_dataset/overlapDataset_LUAD_limma_up.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.8,3,3),
      mainbar.y.label = "Intersection Size LUAD-limma-up", sets.x.label = "DEGs per dataset",
      line.size = 1)
dev.off()

######### down limma  ###########

paired <- read.table("./LUAD/paired/down_limma_LUAD_paired_tumorPurity.txt")
unpaired <- read.table("./LUAD/unpaired/down_limma_LUAD_unpaired_tss_tumorPurity.txt")
all <- read.table("./LUAD/all/down_limma_LUAD_all_tss_tumorPurity.txt")

listInput <- list(paired$V1,unpaired$V1,all$V1)
names(listInput) <- c("paired","unpaired","all")

pdf("./plots_DEA_comparison_dataset/overlapDataset_LUAD_limma_down.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.8,3,3),
      mainbar.y.label = "Intersection Size LUAD-limma-down", sets.x.label = "DEGs per dataset",
      line.size = 1)
dev.off()

######### up edgeR  ###########

paired <- read.table("./LUAD/paired/up_edgeR_LUAD_paired_tumorPurity.txt")
unpaired <- read.table("./LUAD/unpaired/up_edgeR_LUAD_unpaired_tss_tumorPurity.txt")
all <- read.table("./LUAD/all/up_edgeR_LUAD_all_tss_tumorPurity.txt")

listInput <- list(paired$V1,unpaired$V1,all$V1)
names(listInput) <- c("paired","unpaired","all")

pdf("./plots_DEA_comparison_dataset/overlapDataset_LUAD_edgeR_up.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.8,3,3),
      mainbar.y.label = "Intersection Size LUAD-edgeR-up", sets.x.label = "DEGs per dataset",
      line.size = 1)
dev.off()

######### down edgeR  ###########

paired <- read.table("./LUAD/paired/down_edgeR_LUAD_paired_tumorPurity.txt")
unpaired <- read.table("./LUAD/unpaired/down_edgeR_LUAD_unpaired_tss_tumorPurity.txt")
all <- read.table("./LUAD/all/down_edgeR_LUAD_all_tss_tumorPurity.txt")

listInput <- list(paired$V1,unpaired$V1,all$V1)
names(listInput) <- c("paired","unpaired","all")

pdf("./plots_DEA_comparison_dataset/overlapDataset_LUAD_edgeR_down.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.8,3,3),
      mainbar.y.label = "Intersection Size LUAD-edgeR-down", sets.x.label = "DEGs per dataset",
      line.size = 1)
dev.off()

######### up edgeR-TCGAb  ###########

paired <- read.table("./LUAD/paired/up_edgeR_TCGA_LUAD_paired_tumorPurity.txt")
unpaired <- read.table("./LUAD/unpaired/up_edgeR_TCGA_LUAD_unpaired_tumorPurity.txt")
all <- read.table("./LUAD/all/up_edgeR_TCGA_LUAD_all_tumorPurity.txt")

listInput <- list(paired$V1,unpaired$V1,all$V1)
names(listInput) <- c("paired","unpaired","all")

pdf("./plots_DEA_comparison_dataset/overlapDataset_LUAD_edgeR-TCGAb_up.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.8,3,3),
      mainbar.y.label = "Intersection Size LUAD-edgeR-TCGAb-up", sets.x.label = "DEGs per dataset",
      line.size = 1)
dev.off()

######### down edgeR-TCGAb  ###########

paired <- read.table("./LUAD/paired/down_edgeR_TCGA_LUAD_paired_tumorPurity.txt")
unpaired <- read.table("./LUAD/unpaired/down_edgeR_TCGA_LUAD_unpaired_tumorPurity.txt")
all <- read.table("./LUAD/all/down_edgeR_TCGA_LUAD_all_tumorPurity.txt")

listInput <- list(paired$V1,unpaired$V1,all$V1)
names(listInput) <- c("paired","unpaired","all")

pdf("./plots_DEA_comparison_dataset/overlapDataset_LUAD_edgeR-TCGAb_down.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.8,3,3),
      mainbar.y.label = "Intersection Size LUAD-edgeR-TCGAb-down", sets.x.label = "DEGs per dataset",
      line.size = 1)
dev.off()


##########################################################3
#LUSC
####################################################

######### up limma  ###########

paired <- read.table("./LUSC/paired/up_limma_LUSC_paired_tumorPurity.txt")
unpaired <- read.table("./LUSC/unpaired/up_limma_LUSC_unpaired_tss_tumorPurity.txt")
all <- read.table("./LUSC/all/up_limma_LUSC_all_tss_tumorPurity.txt")

listInput <- list(paired$V1,unpaired$V1,all$V1)
names(listInput) <- c("paired","unpaired","all")

pdf("./plots_DEA_comparison_dataset/overlapDataset_LUSC_limma_up.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.8,3,3),
      mainbar.y.label = "Intersection Size LUSC-limma-up", sets.x.label = "DEGs per dataset",
      line.size = 1)
dev.off()


######### down limma  ###########

paired <- read.table("./LUSC/paired/down_limma_LUSC_paired_tumorPurity.txt")
unpaired <- read.table("./LUSC/unpaired/down_limma_LUSC_unpaired_tss_tumorPurity.txt")
all <- read.table("./LUSC/all/down_limma_LUSC_all_tss_tumorPurity.txt")

listInput <- list(paired$V1,unpaired$V1,all$V1)
names(listInput) <- c("paired","unpaired","all")

pdf("./plots_DEA_comparison_dataset/overlapDataset_LUSC_limma_down.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.8,3,3),
      mainbar.y.label = "Intersection Size LUSC-limma-down", sets.x.label = "DEGs per dataset",
      line.size = 1)
dev.off()

######### up edgeR  ###########

paired <- read.table("./LUSC/paired/up_edgeR_LUSC_paired_tumorPurity.txt")
unpaired <- read.table("./LUSC/unpaired/up_edgeR_LUSC_unpaired_tss_tumorPurity.txt")
all <- read.table("./LUSC/all/up_edgeR_LUSC_all_tss_tumorPurity.txt")

listInput <- list(paired$V1,unpaired$V1,all$V1)
names(listInput) <- c("paired","unpaired","all")

pdf("./plots_DEA_comparison_dataset/overlapDataset_LUSC_edgeR_up.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.8,3,3),
      mainbar.y.label = "Intersection Size LUSC-edgeR-up", sets.x.label = "DEGs per dataset",
      line.size = 1)
dev.off()

######### down edgeR  ###########

paired <- read.table("../../LUSC/paired/edgeR/down_edgeR_LUSC_paired_tumorPurity.txt")
unpaired <- read.table("../../LUSC/unpaired/edgeR/down_edgeR_LUSC_unpaired_tss_tumorPurity.txt")
all <- read.table("../../LUSC/all/edgeR/down_edgeR_LUSC_all_tss_tumorPurity.txt")

listInput <- list(paired$V1,unpaired$V1,all$V1)
names(listInput) <- c("paired","unpaired","all")

pdf("./plots_DEA_comparison_dataset/overlapDataset_LUSC_edgeR_down.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.8,3,3),
      mainbar.y.label = "Intersection Size LUSC-edgeR-down", sets.x.label = "DEGs per dataset",
      line.size = 1)
dev.off()


######### up edgeR-TCGAb  ###########

paired <- read.table("../../LUSC/paired/edgeR_TCGA/up_edgeR_TCGA_LUSC_paired_tumorPurity.txt")
unpaired <- read.table("../../LUSC/unpaired/edgeR_TCGA/up_edgeR_TCGA_LUSC_unpaired_tumorPurity.txt")
all <- read.table("../../LUSC/all/edgeR_TCGA/up_edgeR_TCGA_LUSC_all_tumorPurity.txt")

listInput <- list(paired$V1,unpaired$V1,all$V1)
names(listInput) <- c("paired","unpaired","all")

pdf("./plots_DEA_comparison_dataset/overlapDataset_LUSC_edgeR-TCGAb_up.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.8,3,3),
      mainbar.y.label = "Intersection Size LUSC-edgeR-TCGAb-up", sets.x.label = "DEGs per dataset",
      line.size = 1)
dev.off()

######### down edgeR-TCGAb  ###########

paired <- read.table("../../LUSC/paired/edgeR_TCGA/down_edgeR_TCGA_LUSC_paired_tumorPurity.txt")
unpaired <- read.table("../../LUSC/unpaired/edgeR_TCGA/down_edgeR_TCGA_LUSC_unpaired_tumorPurity.txt")
all <- read.table("../../LUSC/all/edgeR_TCGA/down_edgeR_TCGA_LUSC_all_tumorPurity.txt")

listInput <- list(paired$V1,unpaired$V1,all$V1)
names(listInput) <- c("paired","unpaired","all")

pdf("./plots_DEA_comparison_dataset/overlapDataset_LUSC_edgeR-TCGAb_down.pdf",width=10, height=7, onefile = FALSE)
upset(fromList(listInput), order.by = "freq",text.scale = c(2,3,2,1.8,3,3),
      mainbar.y.label = "Intersection Size LUSC-edgeR-TCGAb-down", sets.x.label = "DEGs per dataset",
      line.size = 1)
dev.off()