

#------------ LUAD-up-paired vs LUSC-up-all-----------

up_paired_LUAD <- read.table("../no_overlap_up_paired_LUAD_only.txt")
up_all_LUSC <- read.table("../no_overlap_up_all_LUSC_only.txt")

genes_to_remove <- intersect(up_paired_LUAD$V1,up_all_LUSC$V1)
new_up_paired_LUAD <- subset(up_paired_LUAD,!up_paired_LUAD$V1 %in% genes_to_remove)
write.table(new_up_paired_LUAD,"new_up_paired_LUAD.txt",quote=FALSE, row.names = FALSE,col.names = FALSE)

new_up_all_LUSC <- subset(up_all_LUSC,!up_all_LUSC$V1 %in% genes_to_remove)
write.table(new_up_all_LUSC,"new_up_all_LUSC.txt",quote=FALSE, row.names = FALSE,col.names = FALSE)
length(intersect(new_up_paired_LUAD$V1,new_up_all_LUSC$V1))

#--------LUAD-down-paired vs LUSC-down-all------------------------

down_paired_LUAD <- read.table("../no_overlap_down_paired_LUAD_only.txt")
down_all_LUSC <- read.table("../no_overlap_down_all_LUSC_only.txt")

genes_to_remove <- intersect(down_paired_LUAD$V1,down_all_LUSC$V1)
new_down_paired_LUAD <- subset(down_paired_LUAD,!down_paired_LUAD$V1 %in% genes_to_remove)
write.table(new_down_paired_LUAD,"new_down_paired_LUAD.txt",quote=FALSE, row.names = FALSE,col.names = FALSE)

new_down_all_LUSC <- subset(down_all_LUSC,!down_all_LUSC$V1 %in% genes_to_remove)
write.table(new_down_all_LUSC,"new_down_all_LUSC.txt",quote=FALSE, row.names = FALSE,col.names = FALSE)
length(intersect(new_down_paired_LUAD$V1,new_down_all_LUSC$V1))

#------LUSC_up_paired vs LUAD_up_all-------------------------------

up_paired_LUSC <- read.table("../no_overlap_up_paired_LUSC_only.txt")
up_all_LUAD <- read.table("../no_overlap_up_all_LUAD_only.txt")

genes_to_remove <- intersect(up_paired_LUSC$V1,up_all_LUAD$V1)
new_up_paired_LUSC <- subset(up_paired_LUSC,!up_paired_LUSC$V1 %in% genes_to_remove)
write.table(new_up_paired_LUSC,"new_up_paired_LUSC.txt",quote=FALSE, row.names = FALSE,col.names = FALSE)

new_up_all_LUAD <- subset(up_all_LUAD,!up_all_LUAD$V1 %in% genes_to_remove)
write.table(new_up_all_LUAD,"new_up_all_LUAD.txt",quote=FALSE, row.names = FALSE,col.names = FALSE)
length(intersect(new_up_paired_LUSC$V1,new_up_all_LUAD$V1))

#--------LUSC_down_paired vs LUAD_down_all-------

down_paired_LUSC <- read.table("../no_overlap_down_paired_LUSC_only.txt")
down_all_LUAD <- read.table("../no_overlap_down_all_LUAD_only.txt")

genes_to_remove <- intersect(down_paired_LUSC$V1,down_all_LUAD$V1)
new_down_paired_LUSC <- subset(down_paired_LUSC,!down_paired_LUSC$V1 %in% genes_to_remove)
write.table(new_down_paired_LUSC,"new_down_paired_LUSC.txt",quote=FALSE, row.names = FALSE,col.names = FALSE)

new_down_all_LUAD <- subset(down_all_LUAD,!down_all_LUAD$V1 %in% genes_to_remove)
write.table(new_down_all_LUAD,"new_down_all_LUAD.txt",quote=FALSE, row.names = FALSE,col.names = FALSE)
length(intersect(new_down_paired_LUSC$V1,new_down_all_LUAD$V1))
