library(devtools)
install_github("genomicsclass/maPooling")
library(matrixStats)
library(Biobase)
library(maPooling)
library(genefilter)
library(qvalue)
data(maPooling)
# the rows represent samples and the columns represent mice
pd = pData(maPooling)
pooled = which(rowSums(pd)==12)
individuals=which(rowSums(pd)==1)
## remove replicates
individuals=individuals[-grep("tr", names(individuals))]

# create two measurement matrices repreenting technical replicates
pool = exprs(maPooling)[,pooled]
indiv = exprs(maPooling)[,individuals]
strain = ifelse(grepl("a", rownames(pData(maPooling))),0,1)
g_pool = strain[pooled] # technical replicates
g_indiv = strain[individuals] # biological replicates

## 1 comparing technical and biological variation genome-wide
biologicalsd <- rowSds(indiv[,g_indiv==1])
technicalsd <- rowSds(pool[,g_pool==1])
sum(biologicalsd>technicalsd)/length(biologicalsd)

## 2 two-sample tests, genome wide, with FDR
res <- rowttests(pool, as.factor(g_pool))
qs <- qvalue(res$p.value)$qvalues
sum(qs<0.05)

## 3 can the claims based on pooled data be confirmed
index <- which(qs<0.05)
res2 <- rowttests(indiv[index, ], as.factor(g_indiv))
sum(res2$p.value>0.05)/length(res2$p.value)
## 4 application of the moderated t-test
library(limma)
fit <- lmFit(indiv[index,], design=model.matrix(~as.factor(g_indiv)))
fit <- eBayes(fit)
tt <- topTable(fit, coef=2, number=Inf)
qs3 <- qvalue(tt$P.Value)$qvalues
sum(qs3<0.05)/length(qs3)

