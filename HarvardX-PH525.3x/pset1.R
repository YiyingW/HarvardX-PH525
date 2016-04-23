set.seed(1)
filename <- 'femaleControlsPopulation.csv'
population <- read.csv(filename)
pvals <- replicate(1000, {
  control = sample(population[,1], 12)
  treatment = sample(population[,1], 12)
  t.test(treatment, control)$p.val
  
})
head(pvals)
hist(pvals)


# exercises #3

set.seed(100)
pvals <- replicate(20, {
  cases = rnorm(10, 30, 2)
  controls = rnorm(10, 30, 2)
  t.test(cases, controls)$p.val
  
})
sum(pvals < 0.05)

# exercises #4

set.seed(100)
p_less_0.05 <- replicate(1000, {
  pvals <- replicate(20, {
  cases = rnorm(10, 30, 2)
  controls = rnorm(10, 30, 2)
  t.test(cases, controls)$p.val
    
  })
  sum(pvals < 0.05)
  

  }
)
mean(p_less_0.05)

# exercises #5
sum(p_less_0.05 >= 1)/1000


# bonferroni correction exercises (monte carlo simulation)
set.seed(1)
false_num <- replicate(10000, {
  pvals <- runif(8793, 0, 1)  # to simulate the p-value results of 8793 t-tests for which the null is true
  m <- 8793
  cutoff <- 0.05/m
  sum(pvals < cutoff)
  
})

mean(false_num)

# a matrix version of above code
set.seed(1)
B <- 10000  # number of times the simulation run
m <- 8793  # number of genes/t-tests
alpha <- 0.05  # controlled error rate
pvals <- matrix(runif(B*m, 0,1), B, m)
k <- alpha/m  # adjusted cutoff for calling significant
mistakes <- rowSums(pvals < k)
mean(mistakes > 0)


# use sidak's cutoff and report the FWER 
set.seed(1)
B <- 10000  # number of times the simulation run
m <- 8793  # number of genes/t-tests
alpha <- 0.05  # controlled error rate
pvals <- matrix(runif(B*m, 0,1), B, m)
k <- 1-(1-alpha)^(1/m)  # adjusted cutoff for calling significant
mistakes <- rowSums(pvals < k)
mean(mistakes > 0)










