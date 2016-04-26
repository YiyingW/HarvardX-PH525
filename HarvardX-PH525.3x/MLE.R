## MLE exercises

# question: is there a section of the human cytomegalovirus genome in which the rate of palindromes is higher than expected?
library(dagdata)
data(hcmv)

library(rafalib)
mypar()
# these are the locations of palindromes on the genome of this virus
plot(locations,rep(1,length(locations)),ylab="",yaxt="n")

breaks = seq(0, 4000*round(max(locations)/4000), 4000)
tmp=cut(locations,breaks)
counts=as.numeric(table(tmp))
hist(counts)

probs <- dpois(counts, 4)
likelihood <- prod(probs)
likelihood

logprobs <- dpois(counts, 4, log=T)
loglikelihood <- sum(logprobs)
loglikelihood

## write a function that takes lambda and the vector of counts as input, and returns the 
## log-likelihood. Compute this log-likelihood for lambdas = seq(0, 15, len=300) and make a plot

compute_likelihood <- function(lambda, counts){
  logprobs <- dpois(counts, lambda, log=T)
  loglikelihood <- sum(logprobs)
}

lambdas <- seq(0, 15, len=300)
log_likelihood <- sapply(lambdas, compute_likelihood, counts)

plot(lambdas, log_likelihood)

## solution:
loglikelihood = function(lambda, x){
  sum(dpois(x, lambda, log=T))
  
}

lambdas <- seq(1, 15, len=300)
l <- sapply(lambdas, function(lambda) loglikelihood(lambda, counts))
plot(lambdas, l)
mle = lambdas[which.max(l)]
abline(v=mle)
print(mle)

mean(counts)


## create a plot and see the counts per location
binLocation=(breaks[-1]+breaks[-length(breaks)])/2
plot(binLocation, counts, type="l", xlab=)
binLocation[which.max(counts)]
max(counts)

## we have identified the lcoation with the largest palindrome count, now we want to know
## if by chance we could see a value this big.

lambda = mean(counts[-which.max(counts)])
# what is the probability of seeing a count of 14 or more
1-ppois(14, lambda)
1-ppois(13, lambda)

## create a qq-plot to see if our poisson model is a good fit
ps <- (seq(along=counts) - 0.5) / length(counts)
lambda <- mean(counts[-which.max(counts)])
poisq <- qpois(ps, lambda)
qqplot(poisq, counts)
abline(0, 1)







