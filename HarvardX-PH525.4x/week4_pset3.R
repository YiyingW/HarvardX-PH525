library(GSE5859Subset)
data(GSE5859Subset)

library(rafalib)
library(genefilter)
library(qvalue)

sex <- sampleInfo$group
month <- factor(format(sampleInfo$date, "%m"))
table(sampleInfo$group, month)

res <- rowttests(geneExpression, as.factor(sex))
qvalues <- qvalue(res$p.value)$qvalues
sum(qvalues<0.1)

ind <- which(qvalues < 0.1)
p_small <- geneAnnotation[ind,]
nrow(p_small[which(p_small$CHR%in%c("chrX", "chrY")),])
20/59



###2
index <- geneAnnotation$CHR[qvalues < 0.1]%in%c("chrX", "chrY")
mean(index)

### 3
ind <- which(qvalues < 0.1)
out <- which(qvalues < 0.1 & geneAnnotation$CHR%in%c("chrX","chrY"))
index <- setdiff(ind, out)
res3 <- rowttests(geneExpression[index,], month)
p_value3 <- res3$p.value
sum(p_value3<0.05)/length(p_value3)


## 5
p_values <- sapply(1:nrow(geneExpression), function(i){
       X = model.matrix(~sex+month)
       y = geneExpression[i,]
       fit = lm(y~X-1)
       summary(fit)$coef[2,4]
       })
q_values <- qvalue(p_values)$qvalues
sum(q_values<0.1)

## 6
index <- geneAnnotation$CHR[q_values < 0.1]%in%c("chrX", "chrY")
mean(index)


##7
p_values <- sapply(1:nrow(geneExpression), function(i){
  X = model.matrix(~sex+month)
  y = geneExpression[i,]
  fit = lm(y~X-1)
  summary(fit)$coef[3,4]
})
q_values <- qvalue(p_values)$qvalues
sum(q_values<0.1)








