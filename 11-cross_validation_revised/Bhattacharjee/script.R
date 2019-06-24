#set the working directory

setwd("~/cross-validation/Bhattacharjee")

#######
# gconvert

# convert probe set to gene names and collapse them 
# (multiple probe set match same gene)

library(gProfileR)
library(WGCNA)

data <- read.csv("DatasetA_12600gene.csv")
#remove "gene description" column
data <- data[,-2]
rownames(data) <- data$probe.set
data <- data[,-1]
rowID <- rownames(data)
dim(data)

conversion <- gconvert(rowID, organism = "hsapiens", target="AFFY_HG_U133_PLUS_2")
conversion <- read.csv("conversion_file.csv",row.names = 1)


gene_symbol <- conversion$name[match(rowID,tolower(conversion$alias))]
length(gene_symbol)
# remove all the non-converted probes
conversion_data <- na.omit(cbind(rowID,as.character(gene_symbol)))
new_data <- subset(data,rownames(data) %in% conversion_data[,"rowID"])
dim(new_data)

rowGroup <- conversion_data[,2]
rowID <- rownames(new_data)

# collapse
collapse.object=collapseRows(datET=new_data, rowGroup=rowGroup, rowID=rowID)

dat1Collapsed=data.frame(collapse.object$group2row, collapse.object$datETcollapsed)
# check if the gene names are unique
nrow(dat1Collapsed)
length(unique(dat1Collapsed$group))
length(unique(dat1Collapsed$selectedRowID))

# remove group and selectedRowID columns
dat1Collapsed <- dat1Collapsed[,- c(1,2)]

# retain just adenocarcinoma and squamous samples
dat1Collapsed <- dat1Collapsed[,union(grep("AD",colnames(dat1Collapsed)),grep("SQ",colnames(dat1Collapsed)))]
write.csv(dat1Collapsed,"collapsed_data.csv")

#####################################################################

dat1Collapsed <- read.csv("collapsed_data.csv", row.names = 1)

# retain only the genes of interest
gene_list <- c("MUC5B","HABP2","MUC21","KCNK5","ICA1","CSTA","P2RY1","ANXA8","FZD7","ITGA6",
               "CHST7","ACOX2","ALDOC","AQP5","ARSE","FABP5","SLC2A9", "NRCAM","AGR2", "SPDEF")

dat1Collapsed <- dat1Collapsed[which(rownames(dat1Collapsed) %in% gene_list),]

dat1Collapsed_scaled <- scale(dat1Collapsed)


################
# clustering
#################
library(gplots)

mydist=function(c) {dist(c,method="euclidean")}
myclust=function(c) {hclust(c,method="complete")}
# plot heatmap
clab <- rep("magenta",ncol(dat1Collapsed))
clab[grep("SQ",colnames(dat1Collapsed))] <- "cyan"

png("heatmap_Bhattacharjee_allGenes.png")
heatmap.2(dat1Collapsed_scaled, scale = "none", distfun = mydist,
          hclustfun = myclust, labCol = FALSE, ColSideColors = clab, 
          trace = "none", density.info = "none", margins = c(12,10))
par(lend = 1) 
legend("topright", legend = c("LUAD","LUSC"), col = c("magenta","cyan"),
       lty = 1,lwd = 10, border = FALSE, bty="n", y.intersp = 0.7, cex=0.9)
dev.off()

pdf("heatmap_Bhattacharjee_allGenes.pdf")
heatmap.2(dat1Collapsed_scaled, scale = "none", distfun = mydist,
          hclustfun = myclust, labCol = FALSE, ColSideColors = clab, 
          trace = "none", density.info = "none", margins = c(12,10))
par(lend = 1) 
legend("topright", legend = c("LUAD","LUSC"), col = c("magenta","cyan"),
       lty = 1,lwd = 10, border = FALSE, bty="n", y.intersp = 0.7, cex=0.9)
dev.off()
###-------------------
#retain the genes of interest
dat1Collapsed_scaled <- dat1Collapsed_scaled[c('ALDOC','ARSE','ANXA8',"CSTA"),]

png("heatmap_Bhattacharjee.png")
heatmap.2(dat1Collapsed_scaled, scale = "none", distfun = mydist,
          hclustfun = myclust, labCol = FALSE, ColSideColors = clab, 
          trace = "none", density.info = "none", margins = c(10,7))
par(lend = 1) 
legend("topright", legend = c("LUAD","LUSC"), col = c("magenta","cyan"),
       lty = 1,lwd = 10, border = FALSE, bty="n", y.intersp = 0.7, cex=0.9)
dev.off()


margins = c(10,7)

pdf("heatmap_Bhattacharjee.pdf")
heatmap.2(dat1Collapsed_scaled, scale = "none", distfun = mydist,
          hclustfun = myclust, labCol = FALSE, ColSideColors = clab, 
          trace = "none", density.info = "none", margins = c(10,7))
par(lend = 1) 
legend("topright", legend = c("LUAD","LUSC"), col = c("magenta","cyan"),
       lty = 1,lwd = 10, border = FALSE, bty="n", y.intersp = 0.7, cex=0.9)
dev.off()


margins = c(10,7)
