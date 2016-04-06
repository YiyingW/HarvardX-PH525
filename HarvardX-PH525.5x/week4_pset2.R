# the variance of the average of a set N of normal variables with mean 0 and variance 1
# is 1/N. This can be demonstrated in R with a simulation of 10000 averages for N=10.
var(rowMeans(matrix(rnorm(10000 * 10, 0, 1), ncol=10)))

# suppose we have 10 Normal variables with mean 0 and variance 1, but each randome variable
# has a correlation of 0.7 with the other variables 
# create such random variables in R using multivariate normal distribution
library(MASS)
Sigma = matrix(.7, ncol=10, nrow=10)
diag(Sigma) = 1
averages <- apply(mvrnorm(n=10000,mu=rep(0,10),Sigma=Sigma), 1, mean)



## Gene Set Analysis with the roast algorithm 
library(GEOquery)
e <- getGEO("GSE34313")[[1]]
# the groups are defined in the characteristics_ch1.2 variable
e$condition <- e$characteristics_ch1.2
# rename the conditions
levels(e$condition) <- c("dex24", "dex4", "control")
# e$condition is equivalent to pData(e)$condition

# this data includes GO terms for every gene. 
names(fData(e))
fData(e)$GO_ID[1:4]

# compare the control samples to those treated with dexamethasone for 4 hours
lvls <- c("control", "dex4")
es <- e[,e$condition %in% lvls] # a subset of e containing control and dex4
es$condition <- factor(es$condition, levels = lvls)

# run the linear model in limma
# the top genes are common immune-response genes.
library(limma)
library(qvalue)
design <- model.matrix(~es$condition)
fit <- lmFit(es, design=design)
fit <- eBayes(fit)
topTable(fit)[, c(6,7,18,22)]

# use the ROAST method for gene set testing. we can test a single gene set by looking up the genes
# which contain a certain GO ID and providing this to the roast function

# Immune response
set.seed(1)
idx <- grep("GO:0006955", fData(es)$GO_ID)
length(idx)
r1 <- roast(es, idx, design)
r1

# question
set.seed(1)
idx <- grep("GO:0045454", fData(es)$GO_ID)
length(idx)
r2 <- roast(es, idx, design)
r2


## Inference for multiple gene sets
library(org.Hs.eg.db)  # package used to gather the gene set information
org.Hs.egGO2EG
go2eg <- as.list(org.Hs.egGO2EG)

# unlist the list, gets matches for each Entrez gene ID to the index in the ExpressionSet.
govector <- unlist(go2eg)
golengths <- sapply(go2eg, length)
head(fData(es)$GENE)

idxvector <- match(govector, fData(es)$GENE)
table(is.na(idxvector))
idx <- split(idxvector, rep(names(go2eg), golengths))
go2eg[[1]]
fData(es)$GENE[idx[[1]]]

# clean the list to remove NA and remove gene sets which have less than 10 genes
idxclean <- lapply(idx, function(x) x[!is.na(x)])
idxlengths <- sapply(idxclean, length)
idxsub <- idxclean[idxlengths > 10]
length(idxsub)

set.seed(1)
r3 <- mroast(es, idxsub, design)
head(r3)
r3[which.max(r3$PropUp),]

# rerun mroast to only use gene sets that have 50 or more than 50 genes
idxclean <- lapply(idx, function(x) x[!is.na(x)])
idxlengths <- sapply(idxclean, length)
idxsub <- idxclean[idxlengths >= 50]

set.seed(1)
r4 <- mroast(es, idxsub, design)
r4[which.max(r4$PropUp),]

## Annotation check
library(GO.db)
library(AnnotationDbi)
select(GO.db, keys="GO:0000776", columns=c("GOID","TERM"), keytype="GOID")

library(GO.db)
columns(GO.db)
keytypes(GO.db)
select(GO.db, keys="GO:0000776",columns="TERM") 









