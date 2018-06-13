library(WGCNA)
m=60
# number of rows (e.g. variables or probes on a microarray) 
n=500
# seed module eigenvector for the simulateModule function
MEtrue=rnorm(m)
# numeric data frame of n rows and m columns
datNumeric=data.frame(t(simulateModule(MEtrue,n)))
RowIdentifier=paste("Probe", 1:n, sep="")
ColumnName=paste("Sample",1:m, sep="")
dimnames(datNumeric)[[2]]=ColumnName
# Let us now generate a data frame whose first column contains the rowID
dat1=data.frame(RowIdentifier, datNumeric)
#we simulate a vector with n/5 group labels, i.e. each row group corresponds to 5 rows
rowGroup=rep(  paste("Group",1:(n/5),  sep=""), 5 )

# Typical Input Data 
# Since the first column of dat1 contains the RowIdentifier, we use the following code
datET=dat1[,-1]
rowID=dat1[,1]

# assign row names according to the RowIdentifier 
dimnames(datET)[[1]]=rowID
# run the function and save it in an object

collapse.object=collapseRows(datET=datET, rowGroup=rowGroup, rowID=rowID)

# this creates the collapsed data where 
# the first column contains the group name
# the second column reports the corresponding selected row name (the representative)
# and the remaining columns report the values of the representative row
dat1Collapsed=data.frame( collapse.object$group2row, collapse.object$datETcollapsed)
dat1Collapsed[1:5,1:5]
