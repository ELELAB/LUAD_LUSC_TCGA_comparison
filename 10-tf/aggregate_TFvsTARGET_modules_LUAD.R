#remember to setwd("~/10-tf")
#extract m1_luad
m1_luad <- read.table("../8-coexpression/CEMiTool_summary/modules/m1_luad.txt")
m1_luad_list <-unlist(m1_luad)
trrust <- read.table("trrust_rawdata.human.tsv", sep="\t")
trrust_d <- as.data.frame.matrix(trrust)
#retain only TF - module 1 luad
trrust_m1luad <-subset(trrust_d, V1  %in% m1_luad_list, select=c(V1,V2,V3,V4))
#secon filter to retain TARGETS - module 1 luad
trrust_m1luad_final <-subset(trrust_m1luad, V2   %in%  m1_luad_list, select=c(V1,V2,V3,V4))
names(trrust_m1luad_final) <- c("TF", "TARGET", "EFFECT", "PMID")
write.table(trrust_m1luad_final,"./TF-TARGET_int/m1_luad_TFvsTARGET.txt", sep="\t", quote=FALSE, row.names = FALSE)

#extract m2_luad
m2_luad <- read.table("../8-coexpression/CEMiTool_summary/modules/m2_luad.txt")
m2_luad_list <-unlist(m2_luad)
#trrust <- read.table("trrust_rawdata.human.tsv", sep="\t")
#trrust_d <- as.data.frame.matrix(trrust)
#retain only TF - module 2 luad
trrust_m2luad <-subset(trrust_d, V1  %in% m2_luad_list, select=c(V1,V2,V3,V4))
#secon filter to retain TARGETS - module 2 luad
trrust_m2luad_final <-subset(trrust_m2luad, V2   %in%  m2_luad_list, select=c(V1,V2,V3,V4))
names(trrust_m2luad_final) <- c("TF", "TARGET", "EFFECT", "PMID")
write.table(trrust_m2luad_final,"./TF-TARGET_int/m2_luad_TFvsTARGET.txt", sep="\t", quote=FALSE, row.names = FALSE)

#extract m3_luad
m3_luad <- read.table("../8-coexpression/CEMiTool_summary/modules/m3_luad.txt")
m3_luad_list <-unlist(m3_luad)
#trrust <- read.table("trrust_rawdata.human.tsv", sep="\t")
#trrust_d <- as.data.frame.matrix(trrust)
#retain only TF - module 3 luad
trrust_m3luad <-subset(trrust_d, V1  %in% m3_luad_list, select=c(V1,V2,V3,V4))
#secon filter to retain TARGETS - module 3 luad
trrust_m3luad_final <-subset(trrust_m3luad, V2   %in%  m3_luad_list, select=c(V1,V2,V3,V4))
names(trrust_m3luad_final) <- c("TF", "TARGET", "EFFECT", "PMID")
write.table(trrust_m3luad_final,"./TF-TARGET_int/m3_luad_TFvsTARGET.txt", sep="\t", quote=FALSE, row.names = FALSE)

#extract m4_luad
m4_luad <- read.table("../8-coexpression/CEMiTool_summary/modules/m4_luad.txt")
m4_luad_list <-unlist(m4_luad)
#trrust <- read.table("trrust_rawdata.human.tsv", sep="\t")
#trrust_d <- as.data.frame.matrix(trrust)
#retain only TF - module 4 luad
trrust_m4luad <-subset(trrust_d, V1  %in% m4_luad_list, select=c(V1,V2,V3,V4))
#secon filter to retain TARGETS - module 4 luad
trrust_m4luad_final <-subset(trrust_m4luad, V2   %in%  m4_luad_list, select=c(V1,V2,V3,V4))
names(trrust_m4luad_final) <- c("TF", "TARGET", "EFFECT", "PMID")
write.table(trrust_m4luad_final,"./TF-TARGET_int/m4_luad_TFvsTARGET.txt", sep="\t", quote=FALSE, row.names = FALSE)


