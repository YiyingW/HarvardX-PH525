# Factor Analysis and PCA excercise 

library(Biobase)
library(GSE5859Subset)
library(rafalib)
library(RColorBrewer)
data(GSE5859Subset)

y = geneExpression - rowMeans(geneExpression)

## Compute and plot an image of the correlation for each sample. 
## Make two image plots of these correlations. In the first one, plot the correlation as image. 
## In the second, order the samples by date and then plot the image of the correlation. 
## The only difference in these plots is the order in which the samples are plotted.

## 1
# the first plot
cors <- cor(y)
n <- ncol(geneExpression)
cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
mypar()
image(1:n, 1:n, cors, col=cols, xlab="", ylab="", zlim=c(-1,1))

# the second plot
# order the sample by date use a new index order
index <- order(sampleInfo$date)
cors <- cor(y[,index])
n <- ncol(geneExpression)
cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
mypar()
image(1:n, 1:n, cors, col=cols, xlab="", ylab="", zlim=c(-1,1))

## 2
# there are at least two hidden factors. Use PCA estimate these two factors. 
pcs <- svd(y)$v[,1:2]


## 3
# plot each of the estimated factor ordered by date. Use color to denote month. The first factor is clearly
# related to date. 

o <- order(sampleInfo$date)
month <- as.factor(format(sampleInfo$date, "%m"))
cols <- as.numeric(month)
plot(s$v[o,1], pch=21, cex = 1.25, bg=cols[o], ylab="First PC", xlab="")
legend("topleft", c("Month 06", "Month 10"), col = 1:2, pch=16, box.lwd=0)


## 4
s <- svd(y)
plot(s$d^2/sum(s$d^2), ylab="variance explained", xlab="PC")

## 5
# which PC most correlates (negative or positive correlation) with month
correlations <- sapply(1: ncol(s$v), function(i){
  cor(s$v[,i], as.numeric(month))
  })
which.max(abs(correlations))
correlations[1]

## 6
# which PC most correlates with sex
sex <- sampleInfo$group
correlations <- sapply(1: ncol(s$v), function(i){
  cor(s$v[,i], as.numeric(sex))
})
which.max(abs(correlations))
correlations[1]


## 7
# instead of using month, add the two estimated factors in Factor Analysis
# to the linear model and apply this model to each gene, compute q-values for
# the sex difference
p_values <- sapply(1:nrow(geneExpression), function(i){
  X = model.matrix(~sex+s$v[,1:2])
  y = geneExpression[i,]
  fit = lm(y~X-1)
  summary(fit)$coef[2,4]
})
q_values <- qvalue(p_values)$qvalues
sum(q_values<0.1)

index <- geneAnnotation$CHR[q_values < 0.1]%in%c("chrX", "chrY")
mean(index)






