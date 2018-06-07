
## LUAD paired####

dir1 <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/VennDiagram_three_methods/VennDiagram_LUAD-LUSC/final_genes_list/tumor_Purity/"

up_paired_LUAD <- read.table("new_up_paired_LUAD.txt")
up_paired_LUAD1 <- read.table(paste0(dir1,"new_up_paired_LUAD.txt"))
length(up_paired_LUAD$V1)
length(up_paired_LUAD1$V1)
intersect <- intersect(up_paired_LUAD$V1,up_paired_LUAD1$V1)
length(intersect)
genes_paired_LUAD <- up_paired_LUAD$V1[which(!up_paired_LUAD$V1 %in% intersect)]
length(intersect(up_paired_LUAD1$V1,genes_paired_LUAD))
write.table(genes_paired_LUAD,"more_genes_paired_up_LUAD.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
genes_paired_LUAD_less <- up_paired_LUAD1$V1[which(!up_paired_LUAD1$V1 %in% intersect)]
length(genes_paired_LUAD_less)
length(intersect(up_paired_LUAD$V1,genes_paired_LUAD_less))
write.table(genes_paired_LUAD_less,"less_genes_paired_up_LUAD.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)


down_paired_LUAD <- read.table("new_down_paired_LUAD.txt")
down_paired_LUAD1 <- read.table(paste0(dir1,"new_down_paired_LUAD.txt"))
length(down_paired_LUAD$V1)
length(down_paired_LUAD1$V1)
intersect <- intersect(down_paired_LUAD$V1,down_paired_LUAD1$V1)
length(intersect)
genes_paired_LUAD <- down_paired_LUAD$V1[which(!down_paired_LUAD$V1 %in% intersect)]
length(intersect(down_paired_LUAD1$V1,genes_paired_LUAD))
write.table(genes_paired_LUAD,"more_genes_paired_down_LUAD.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
genes_paired_LUAD_less <- down_paired_LUAD1$V1[which(!down_paired_LUAD1$V1 %in% intersect)]
length(genes_paired_LUAD_less)
length(intersect(down_paired_LUAD$V1,genes_paired_LUAD_less))
write.table(genes_paired_LUAD_less,"less_genes_paired_down_LUAD.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)


## LUAD all

up_all_LUAD <- read.table("new_up_all_LUAD.txt")
up_all_LUAD1 <- read.table(paste0(dir1,"new_up_all_LUAD.txt"))
length(up_all_LUAD$V1)
length(up_all_LUAD1$V1)
intersect <- intersect(up_all_LUAD$V1,up_all_LUAD1$V1)
length(intersect)
genes_LUAD <- up_all_LUAD$V1[which(!up_all_LUAD$V1 %in% intersect)]
length(genes_LUAD)
length(intersect(up_all_LUAD1$V1,genes_LUAD))
write.table(genes_LUAD,"more_genes_all_up_LUAD.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
genes_all_LUAD_less <- up_all_LUAD1$V1[which(!up_all_LUAD1$V1 %in% intersect)]
length(genes_all_LUAD_less)
length(intersect(up_all_LUAD$V1,genes_all_LUAD_less))
write.table(genes_all_LUAD_less,"less_genes_all_up_LUAD.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)


down_all_LUAD <- read.table("new_down_all_LUAD.txt")
down_all_LUAD1 <- read.table(paste0(dir1,"new_down_all_LUAD.txt"))
length(down_all_LUAD$V1)
length(down_all_LUAD1$V1)
intersect <- intersect(down_all_LUAD$V1,down_all_LUAD1$V1)
length(intersect)
genes_LUAD <- down_all_LUAD$V1[which(!down_all_LUAD$V1 %in% intersect)]
length(genes_LUAD)
length(intersect(down_all_LUAD1$V1,genes_LUAD))
write.table(genes_LUAD,"more_genes_all_down_LUAD.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
genes_all_LUAD_less <- down_all_LUAD1$V1[which(!down_all_LUAD1$V1 %in% intersect)]
length(genes_all_LUAD_less)
length(intersect(down_all_LUAD$V1,genes_all_LUAD_less))
write.table(genes_all_LUAD_less,"less_genes_all_down_LUAD.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)


## LUSC paired

up_paired_LUSC <- read.table("new_up_paired_LUSC.txt")
up_paired_LUSC1 <- read.table(paste0(dir1,"new_up_paired_LUSC.txt"))
length(up_paired_LUSC$V1)
length(up_paired_LUSC1$V1)
intersect <- intersect(up_paired_LUSC$V1,up_paired_LUSC1$V1)
length(intersect)
genes_LUSC <- up_paired_LUSC$V1[which(!up_paired_LUSC$V1 %in% intersect)]
length(genes_LUSC)
length(intersect(up_paired_LUSC1$V1,genes_LUSC))
write.table(genes_LUSC,"more_genes_paired_up_LUSC.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
genes_paired_LUSC_less <- up_paired_LUSC1$V1[which(!up_paired_LUSC1$V1 %in% intersect)]
length(genes_paired_LUSC_less)
length(intersect(up_paired_LUSC$V1,genes_paired_LUSC_less))
write.table(genes_paired_LUSC_less,"less_genes_paired_up_LUSC.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)


down_paired_LUSC <- read.table("new_down_paired_LUSC.txt")
down_paired_LUSC1 <- read.table(paste0(dir1,"new_down_paired_LUSC.txt"))
length(down_paired_LUSC$V1)
length(down_paired_LUSC1$V1)
intersect <- intersect(down_paired_LUSC$V1,down_paired_LUSC1$V1)
length(intersect)
genes_LUSC <- down_paired_LUSC$V1[which(!down_paired_LUSC$V1 %in% intersect)]
length(genes_LUSC)
length(intersect(down_paired_LUSC1$V1,genes_LUSC))
write.table(genes_LUSC,"more_genes_paired_down_LUSC.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
genes_paired_LUSC_less <- down_paired_LUSC1$V1[which(!down_paired_LUSC1$V1 %in% intersect)]
length(genes_paired_LUSC_less)
length(intersect(down_paired_LUSC$V1,genes_paired_LUSC_less))
write.table(genes_paired_LUSC_less,"less_genes_paired_down_LUSC.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)


### LUSC all

up_all_LUSC <- read.table("new_up_all_LUSC.txt")
up_all_LUSC1 <- read.table(paste0(dir1,"new_up_all_LUSC.txt"))
length(up_all_LUSC$V1)
length(up_all_LUSC1$V1)
intersect <- intersect(up_all_LUSC$V1,up_all_LUSC1$V1)
length(intersect)
genes_LUSC <- up_all_LUSC$V1[which(!up_all_LUSC$V1 %in% intersect)]
length(genes_LUSC)
length(intersect(up_all_LUSC1$V1,genes_LUSC))
write.table(genes_LUSC,"more_genes_all_up_LUSC.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
genes_all_LUSC_less <- up_all_LUSC1$V1[which(!up_all_LUSC1$V1 %in% intersect)]
length(genes_all_LUSC_less)
length(intersect(up_all_LUSC$V1,genes_all_LUSC_less))
write.table(genes_all_LUSC_less,"less_genes_all_up_LUSC.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

down_all_LUSC <- read.table("new_down_all_LUSC.txt")
down_all_LUSC1 <- read.table(paste0(dir1,"new_down_all_LUSC.txt"))
length(down_all_LUSC$V1)
length(down_all_LUSC1$V1)
intersect <- intersect(down_all_LUSC$V1,down_all_LUSC1$V1)
length(intersect)
genes_LUSC <- down_all_LUSC$V1[which(!down_all_LUSC$V1 %in% intersect)]
length(genes_LUSC)
length(intersect(down_all_LUSC1$V1,genes_LUSC))
write.table(genes_LUSC,"more_genes_all_down_LUSC.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
genes_all_LUSC_less <- down_all_LUSC1$V1[which(!down_all_LUSC1$V1 %in% intersect)]
length(genes_all_LUSC_less)
length(intersect(down_all_LUSC$V1,genes_all_LUSC_less))
write.table(genes_all_LUSC_less,"less_genes_all_down_LUSC.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
