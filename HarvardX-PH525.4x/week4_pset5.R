# Surrogate Variable Analysis
library(sva)
library(Biobase)
library(GSE5859Subset)
data(GSE5859Subset)

s <- svd(geneExpression-rowMeans(geneExpression))
cor(sampleInfo$group, s$v[,1]) # the first factor was correlated with outcome of interest

sex <- sampleInfo$group
mod <- model.matrix(~sex)
svafit <- sva(geneExpression, mod)
head(svafit$sv) # estimated factors 

# fit a linear model to estimate the difference between males and females for each gene
# include the factors estimated by sva in the model
p_values <- sapply(1:nrow(geneExpression), function(i){
  X = model.matrix(~sex+svafit$sv)
  y = geneExpression[i,]
  fit = lm(y~X-1)
  summary(fit)$coef[2,4]
})
q_values <- qvalue(p_values)$qvalues
sum(q_values<0.1)

index <- geneAnnotation$CHR[q_values < 0.1]%in%c("chrX", "chrY")
mean(index)







