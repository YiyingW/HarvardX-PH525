tmpfile <- tempfile()
tmpdir <- tempdir()
download.file("http://seanlahman.com/files/database/lahman-csv_2014-02-14.zip",tmpfile)
## this shows us files
filenames <- unzip(tmpfile, list = T)
players <- read.csv(unzip(tmpfile, files="Batting.csv", exdir=tmpdir), as.is=T)
unlink(tmpdir)
file.remove(tmpfile)

# use dplyr to obtain the necessary information to perform a hierarchical model
# this dplyr commands gives the batting averages (AVG) for players with more than 500 AB in 2012
# filter(players, yearID==2012) %>% mutate(AVG=H/AB) %>% filter(AB>=500) %>% select(AVG)
library(dplyr)
# obtain all the batting averages from 2010, 2011, 2012 and removing rows with AB < 500
b_avg <- filter(players, yearID==2012|yearID==2011|yearID==2010) %>% mutate(AVG=H/AB) %>% filter(AB>=500) %>% select(AVG)
mean(b_avg[[1]])
sd(b_avg[[1]])
mypar()
qqnorm(b_avg[[1]])
qqline(b_avg[[1]])
# check if it is approximately normally distributed 

p <- 0.45

sqrt(p*(1-p)/20)

B <- 0.11^2/(0.11^2 + 0.027^2)
mu <- 0.275
E <- mu + (1-B)*(0.45 - mu)

## to install:
library(rafalib)
#install_bioc("SpikeInSubset")
library(Biobase)
library(SpikeInSubset)
data(rma95)
y <- exprs(rma95)


## RNA was obtained from the same background pool to create six replicate samples
## RNA from 16 genes were artificially added in different quantities to each sample
pData(rma95)
g <- factor(rep(0:1, each=3))  # define two groups
# create an index of which rows are associated with the artificially added genes
spike <- rownames(y) %in% colnames(pData(rma95))

# note that only these 16 genes are differentially expressed since these six samples differ only due to random sampling
# perform a t-test on each gene using the rowttests function in the genefilter package
# what proportion of genes with a p-value < 0.01 are not part of the artificially added (false positive)?

library(genefilter)
sig_genes <- rownames(y[rowttests(y, g)$p.value < 0.01,])
1-sum(sig_genes%in%colnames(pData(rma95)))/length(sig_genes)


## compute the within group sample SD for each gene (can use group 1). Based on the p-value<0.01 cutoff, split the genes
## into true positive, false positive, true negatives and false negatives. create a boxplot comparing the sampe SDs

sds <- apply(y[,4:6], 1, sd)
positives <- which(rowttests(y, g)$p.value < 0.01)  # the positives predicted by t-test
real <- which(spike)  # the real 
TP <- intersect(positives, real)
FP <- setdiff(positives, real)
negatives <- which(rowttests(y, g)$p.value >= 0.01)
TN <- setdiff(negatives, real)
FN <- intersect(negatives, real)
boxplot(sds[TP], sds[FP], sds[TN], sds[FN])

## the random variability associated with the sample standard deviation leads to
## t-statistics that are large by chance
## the sample standard deviation we use in the t-test is an estimate and that
## with just a pair of triplicate samples, the variability associated with the 
## denominator in the t-test can be large
## perform the basic limma analysis, the eBayes step uses a hierarchical model
## that provides a new estimate of the gene specific standard error
library(limma)
fit <- lmFit(y, design=model.matrix(~g))
colnames(coef(fit))
fit <- eBayes(fit)
sampleSD = fit$sigma
posteriorSD = sqrt(fit$s2.post)
plot(sampleSD, posteriorSD)

## use these new estimates of SD in the denominator of the t-test and compute p-values
## second coefficient relates to differences between group
pval <- fit$p.value[,2]
sig_genes <- rownames(y[pval < 0.01,])
1-sum(sig_genes%in%colnames(pData(rma95)))/length(sig_genes)


## exploratory data analysis
data(mas133)
e <- exprs(mas133)
plot(e[,1],e[,2],main=paste0("corr=",signif(cor(e[,1],e[,2]),3)),cex=0.5)
k <- 3000
b <- 1000 #a buffer
polygon(c(-b,k,k,-b),c(-b,-b,k,k),col="red",density=0,border="red")

total <- nrow(e)
nrow(e[e[,1]<k & e[,2]<k,])/total

## make the sample plot with log
plot(log2(e[,1]),log2(e[,2]),main=paste0("corr=",signif(cor(log2(e[,1]),log2(e[,2])),2)),cex=0.5)
k <- log2(3000)
b <- log2(0.5)
polygon(c(b,k,k,b),c(b,b,k,k),col="red",density=0,border="red")


## make an MA-plot
e <- log2(exprs(mas133))
plot((e[,1]+e[,2])/2,e[,2]-e[,1],cex=0.5)

sd(e[,2]-e[,1])
