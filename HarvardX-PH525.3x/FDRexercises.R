# define error controlling procedures for experimental data
# make list of genes based on q-values
library(GSE5859Subset)
data(GSE5859Subset)
library(genefilter)

# exercise 1, how many genes have p-values smaller than 0.05
sum(rowttests(geneExpression, as.factor(sampleInfo$group))$p.value<0.05)

# exercise 2, apply the Bonferroni correction to the pvalues above to achieve
# a FWER of 0.05. How many genes are called significant under this procedure?

sum(rowttests(geneExpression, as.factor(sampleInfo$group))$p.value<0.05/nrow(geneExpression))
# we are quite certain had no false positives


# exercise 3
# to define the q-value, order features tested by p-value then compute the FDRs for a list with the
# most significant, the two most significant, the three most significant, etc...
# the FDR of the list with the m most significant tests is defined as the q-value of the m-th most significant feature
# the q-value of a feature is the FDR of the biggest list that includes that gene

pvals <- rowttests(geneExpression, as.factor(sampleInfo$group))$p.value
qvals <-p.adjust(pvals, method='fdr')
sum(qvals < 0.05)
# Here we are saying that we think this list of genes has about 5% false positives



# exercise 4
library(qvalue)
qs <- qvalue(pvals)$qvalue
sum(qs < 0.05)

# exercise 5
qvalue(pvals)$pi0

# exercise 7



set.seed(1)

FP_Bon <- replicate(1000, {
  # create a Monte Carlo Simulation, simulate measurements from 8793 genes for 24 samples
  n <- 24
  m <- 8793
  mat <- matrix(rnorm(n*m), m, n)
  
  # for 500 genes, there is a difference of 2 between cases and controls
  delta <- 2
  positives <- 500
  mat[1:positives, 1:(n/2)] <- mat[1:positives, 1:(n/2)] + delta
  # the null hypothesis is true for 8793-500 genes. 
  # m = 8793, m0 = 8293, m1 = 500
  
  # set the seed at 1, run this experiment 1000 times with a Monte Carlo simulation. For each instance,
  # compute p-values using a t-test and create three lists of genes using:
  # 1. Bonferroni correction to achieve an FWER of 0.05,
  # 2. p-adjust estimates of FDR to achieve an FDR of 0.05,
  # 3. qvalue estimates of FDR to achieve an FDR of 0.05.
  # for each of these three lists, compute the number of false positives in the list and the number of
  # false negatives: genes not in the list that should have been because the null hypothesis is not true.
  m0 <- 8293
  m1 <- 500

  # Bonferroni
  pval_Bon <- rowttests(mat)$p.value
  sig_Bon <- pval_Bon<0.05/m
  # number of false positives, aka, TRUE in sig_Bon but index >= 501
  # FP/FP+TN
  FP <- sum(sig_Bon[501:m])
  TN <- m0
  FP/m0

})

mean(FP_Bon)


# what is the false negative rate if we use Bonferroni

set.seed(1)

FN_Bon <- replicate(1000, {
  # create a Monte Carlo Simulation, simulate measurements from 8793 genes for 24 samples
  n <- 24
  m <- 8793
  mat <- matrix(rnorm(n*m), m, n)
  
  # for 500 genes, there is a difference of 2 between cases and controls
  delta <- 2
  positives <- 500
  mat[1:positives, 1:(n/2)] <- mat[1:positives, 1:(n/2)] + delta
  # the null hypothesis is true for 8793-500 genes. 
  # m = 8793, m0 = 8293, m1 = 500
  
  # set the seed at 1, run this experiment 1000 times with a Monte Carlo simulation. For each instance,
  # compute p-values using a t-test and create three lists of genes using:
  # 1. Bonferroni correction to achieve an FWER of 0.05,
  # 2. p-adjust estimates of FDR to achieve an FDR of 0.05,
  # 3. qvalue estimates of FDR to achieve an FDR of 0.05.
  # for each of these three lists, compute the number of false positives in the list and the number of
  # false negatives: genes not in the list that should have been because the null hypothesis is not true.
  m0 <- 8293
  m1 <- 500
  
  # Bonferroni
  pval_Bon <- rowttests(mat)$p.value
  sig_Bon <- pval_Bon<0.05/m
  # number of false negatives, aka, FALSE in sig_Bon but index <=500
  # FN/TP+FN
  FN <- m1-sum(sig_Bon[1:m1])
  TP <- m1
  FN/(m0)
  
})


set.seed(1)
library(qvalue)
library(genefilter)
n <- 24
m <- 8793
B <- 1000
delta <-2
positives <- 500
g <- factor(rep(c(0,1),each=12))
result <- replicate(B,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.val
  ## q-values from p.adjust
  qvals <-p.adjust(pvals, method='fdr')
  FP1 <- sum(qvals[-(1:positives)]<=0.05)  
  FP1
})

mean(result/(m-positives))


set.seed(1)
library(qvalue)
library(genefilter)
n <- 24
m <- 8793
B <- 1000
delta <-2
positives <- 500
g <- factor(rep(c(0,1),each=12))
result <- replicate(B,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.val
  ## q-values from p.adjust
  qvals <-p.adjust(pvals, method='fdr')
  FP1 <- sum(qvals[-(1:positives)]<=0.05)  
  FN1 <- sum(qvals[1:positives]>0.05)
  c(FP1,FN1)
})
mean(result[2,]/(positives))


set.seed(1)
library(qvalue)
library(genefilter)
n <- 24
m <- 8793
B <- 1000
delta <-2
positives <- 500
g <- factor(rep(c(0,1),each=12))
result <- replicate(B,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.val
  ## q-values from qvalue function
  qvals <- qvalue(pvals)$qvalue
  FP1 <- sum(qvals[-(1:positives)]<=0.05)  
  FP1
})

mean(result/(m-positives))



set.seed(1)
library(qvalue)
library(genefilter)
n <- 24
m <- 8793
B <- 1000
delta <-2
positives <- 500
g <- factor(rep(c(0,1),each=12))
result <- replicate(B,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.val
  ## q-values from qvalue function
  qvals <- qvalue(pvals)$qvalue
  FP1 <- sum(qvals[-(1:positives)]<=0.05)  
  FN1 <- sum(qvals[1:positives]>0.05)
  c(FP1,FN1)
})
mean(result[2,]/(positives))
