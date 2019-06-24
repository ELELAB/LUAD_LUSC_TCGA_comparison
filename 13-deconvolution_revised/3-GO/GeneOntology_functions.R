# '''Thilde Bagger Terkelsen, thilde@cancer.dk, DCRC 15-03-2017'''

library(IRanges)
library(DBI)
library(topGO)


# --------------------------------------------------------------------------------------
# Download Human Ontology Terms and Associated Genes
# --------------------------------------------------------------------------------------

#file <- basename("ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz")
#download.file("ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz", file)


# --------------------------------------------------------------------------------------
# Convert to GO Background
# --------------------------------------------------------------------------------------

#HUGO <- read.delim('goa_human.gaf.gz', skip = 12, sep = "\t", header = FALSE)
#HUGO <- data.frame(HUGO$V5, HUGO$V3)
#colnames(HUGO) <- c("GOID", "gene")
#HUGO <- HUGO[order(HUGO$GOID),]
#HUGO <- split(HUGO$gene, HUGO$GOID)

#save(HUGO, file = "HUGO.RData")

#load("HUGO.RData")

# --------------------------------------------------------------------------------------
# Functions which creates GOobject. 
# The function takes arguments: ont = the ontology to enrich for, univ = a univers of genes, 
# intgenes = a set of genes of interest and my.GO.background = GO-background (dataframe where genes are coupled to GO-terms)
# --------------------------------------------------------------------------------------

GOobject <-function(ont, univ, intgenes, my.GO.background ){
  if(is.null(univ)){
    univ=(unique(unlist(my.GO.background)))
  }
  # factorise genelist
  geneList = factor(as.integer(univ%in%intgenes))
  names(geneList) = univ
  #new topgo object
  GOdata = new("topGOdata", ontology = as.character(ont), allGenes = geneList, annot = annFUN.GO2genes, GO2genes=my.GO.background)
  return(GOdata)
}  


# --------------------------------------------------------------------------------------
# Function which creates GOobject and identifies significant GO-terms. 
# The function takes the arguments: ont = the ontology to enrich for, univ = a univers of genes, 
# intgenes = a set of genes of interest, my.GO.background = GO-background (list object where genes are coupled to GO-terms) 
# and nterms = number of GO hits to display.
# --------------------------------------------------------------------------------------


TOPGO <- function(ont, univ, intgenes, my.GO.background, nterms, my.name) {
  
  GOdata <- GOobject(ont, univ, intgenes, my.GO.background)
  test.stat = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher = getSigGroups(GOdata, test.stat)
  test.stat = new("elimCount", testStatistic = GOFisherTest, name = "Fisher test", cutOff = 0.01)
  resultElim = getSigGroups(GOdata, test.stat)

  res <- GenTable(GOdata, fisher = resultFisher, elim = resultElim, orderBy = "elim", ranksOf = "elim", topNodes = nterms)
  
  # pvalue correction
  FDR_fisher <- data.frame(p.adjust(res$fisher, method = "fdr", n = nrow(res)))
  colnames(FDR_fisher) <- "FDR_fisher"
  FDR_elim <- data.frame(p.adjust(res$elim, method = "fdr", n = nrow(res)))
  colnames(FDR_elim) <- "FDR_elim"
  
  res=cbind(res, FDR_fisher, FDR_elim)
  save(res, file=paste0(my.name, ".RData"))
  return(res)
}



  

