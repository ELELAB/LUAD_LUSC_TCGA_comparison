#set the working directory

setwd("~/cross-validation/GSE33532")

library(gProfileR)
library(WGCNA)

data <- read.csv("final_matrix.csv",row.names = 1)
rowID <- rownames(data)
dim(data)

conversion <- gconvert(rowID, organism = "hsapiens", target="AFFY_HG_U133_PLUS_2")
write.csv <- write.csv(conversion,file = "conversion_file.csv")
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
write.csv(dat1Collapsed,file="collapsed_data.csv")


#########################################################################

dat1Collapsed <- read.csv("collapsed_data.csv",row.names = 1)
# retain only the genes of interest

gene_list <- c("KCNK5","ICA1","CSTA","P2RY1","FZD7","CHST7","ACOX2","ARSE","SLC2A9","SPDEF", "MUC5B","HABP2")

dat1Collapsed <- dat1Collapsed[which(rownames(dat1Collapsed) %in% gene_list),]

################
# clustering
#################
library(gplots)

mydist=function(c) {dist(c,method="euclidean")}
myclust=function(c) {hclust(c,method="complete")}
# plot heatmap
clab <- rep("magenta",ncol(dat1Collapsed))
clab[grep("LUSC",colnames(dat1Collapsed))] <- "cyan"


pdf("heatmap_GSE33532.pdf")
heatmap.2(as.matrix(dat1Collapsed), scale = "none", distfun = mydist,
          hclustfun = myclust, labCol = FALSE, ColSideColors = clab, 
          trace = "none", density.info = "none", margins = c(10,9))
par(lend = 1) 
legend('topright', legend = c("LUAD","LUSC"), col = c("magenta","cyan"),
       lty = 1,lwd = 10, border = FALSE, bty="n", y.intersp = 0.7, cex=0.9)
dev.off()
