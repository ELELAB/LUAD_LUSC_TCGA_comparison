#extract m1_lusc
m1_lusc <- read.table("/Users/ele/elena/CEMiTool_summary/modules/m1_lusc.txt")
m1_lusc_list <-unlist(m1_lusc)
trrust <- read.table("/Users/ele/elena/tf/trrust_rawdata.human.tsv", sep="\t")
trrust_d <- as.data.frame.matrix(trrust)
#retain only TF - module 1 lusc
trrust_m1lusc <-subset(trrust_d, V1  %in% m1_lusc_list, select=c(V1,V2,V3,V4))
#secon filter to retain TARGETS - module 1 lusc
trrust_m1lusc_final <-subset(trrust_m1lusc, V2   %in%  m1_lusc_list, select=c(V1,V2,V3,V4))
names(trrust_m1lusc_final) <- c("TF", "TARGET", "EFFECT", "PMID")
write.table(trrust_m1lusc_final,"/Users/ele/elena/tf/TF-TARGET_int/m1_lusc_TFvsTARGET.txt", sep="\t", quote=FALSE, row.names = FALSE)

#extract m2_lusc
m2_lusc <- read.table("/Users/ele/elena/CEMiTool_summary/modules/m2_lusc.txt")
m2_lusc_list <-unlist(m2_lusc)
#trrust <- read.table("/Users/ele/elena/tf/trrust_rawdata.human.tsv", sep="\t")
#trrust_d <- as.data.frame.matrix(trrust)
#retain only TF - module 2 lusc
trrust_m2lusc <-subset(trrust_d, V1  %in% m2_lusc_list, select=c(V1,V2,V3,V4))
#secon filter to retain TARGETS - module 2 lusc
trrust_m2lusc_final <-subset(trrust_m2lusc, V2   %in%  m2_lusc_list, select=c(V1,V2,V3,V4))
names(trrust_m2lusc_final) <- c("TF", "TARGET", "EFFECT", "PMID")
write.table(trrust_m2lusc_final,"/Users/ele/elena/tf/TF-TARGET_int/m2_lusc_TFvsTARGET.txt", sep="\t", quote=FALSE, row.names = FALSE)

#extract m3_lusc
m3_lusc <- read.table("/Users/ele/elena/CEMiTool_summary/modules/m3_lusc.txt")
m3_lusc_list <-unlist(m3_lusc)
#trrust <- read.table("/Users/ele/elena/tf/trrust_rawdata.human.tsv", sep="\t")
#trrust_d <- as.data.frame.matrix(trrust)
#retain only TF - module 3 lusc
trrust_m3lusc <-subset(trrust_d, V1  %in% m3_lusc_list, select=c(V1,V2,V3,V4))
#secon filter to retain TARGETS - module 3 lusc
trrust_m3lusc_final <-subset(trrust_m3lusc, V2   %in%  m3_lusc_list, select=c(V1,V2,V3,V4))
names(trrust_m3lusc_final) <- c("TF", "TARGET", "EFFECT", "PMID")
write.table(trrust_m3lusc_final,"/Users/ele/elena/tf/TF-TARGET_int/m3_lusc_TFvsTARGET.txt", sep="\t", quote=FALSE, row.names = FALSE)

#extract m4_lusc
m4_lusc <- read.table("/Users/ele/elena/CEMiTool_summary/modules/m4_lusc.txt")
m4_lusc_list <-unlist(m4_lusc)
#trrust <- read.table("/Users/ele/elena/tf/trrust_rawdata.human.tsv", sep="\t")
#trrust_d <- as.data.frame.matrix(trrust)
#retain only TF - module 4 lusc
trrust_m4lusc <-subset(trrust_d, V1  %in% m4_lusc_list, select=c(V1,V2,V3,V4))
#secon filter to retain TARGETS - module 4 lusc
trrust_m4lusc_final <-subset(trrust_m4lusc, V2   %in%  m4_lusc_list, select=c(V1,V2,V3,V4))
names(trrust_m4lusc_final) <- c("TF", "TARGET", "EFFECT", "PMID")
write.table(trrust_m4lusc_final,"/Users/ele/elena/tf/TF-TARGET_int/m4_lusc_TFvsTARGET.txt", sep="\t", quote=FALSE, row.names = FALSE)

#extract m5_lusc
m5_lusc <- read.table("/Users/ele/elena/CEMiTool_summary/modules/m5_lusc.txt")
m5_lusc_list <-unlist(m5_lusc)
#trrust <- read.table("/Users/ele/elena/tf/trrust_rawdata.human.tsv", sep="\t")
#trrust_d <- as.data.frame.matrix(trrust)
#retain only TF - module 5 lusc
trrust_m5lusc <-subset(trrust_d, V1  %in% m5_lusc_list, select=c(V1,V2,V3,V4))
#secon filter to retain TARGETS - module 5 lusc
trrust_m5lusc_final <-subset(trrust_m5lusc, V2   %in%  m5_lusc_list, select=c(V1,V2,V3,V4))
names(trrust_m5lusc_final) <- c("TF", "TARGET", "EFFECT", "PMID")
write.table(trrust_m5lusc_final,"/Users/ele/elena/tf/TF-TARGET_int/m5_lusc_TFvsTARGET.txt", sep="\t", quote=FALSE, row.names = FALSE)

#extract m6_lusc
m6_lusc <- read.table("/Users/ele/elena/CEMiTool_summary/modules/m6_lusc.txt")
m6_lusc_list <-unlist(m6_lusc)
#trrust <- read.table("/Users/ele/elena/tf/trrust_rawdata.human.tsv", sep="\t")
#trrust_d <- as.data.frame.matrix(trrust)
#retain only TF - module 6 lusc
trrust_m6lusc <-subset(trrust_d, V1  %in% m6_lusc_list, select=c(V1,V2,V3,V4))
#secon filter to retain TARGETS - module 6 lusc
trrust_m6lusc_final <-subset(trrust_m6lusc, V2   %in%  m6_lusc_list, select=c(V1,V2,V3,V4))
names(trrust_m6lusc_final) <- c("TF", "TARGET", "EFFECT", "PMID")
write.table(trrust_m6lusc_final,"/Users/ele/elena/tf/TF-TARGET_int/m6_lusc_TFvsTARGET.txt", sep="\t", quote=FALSE, row.names = FALSE)
