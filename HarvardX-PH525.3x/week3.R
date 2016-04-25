## the probability of conceiving a girl is 0.49, what is the probability that
## a family with 4 children has 2 girls and 2 boys

p = 0.49
dbinom(2, 4, p)

## what is the probability that a family with 10 children has 4 girls and 6 boys
pbinom(4, 10, p) - pbinom(3, 10, p)
# or
dbinom(4, 10, p)


## the genome has 3 billion bases, 20% C, 20% G, 30% A, 30% T. you take a random interval of 20 bases,
## what is the probability that the GC-content is strictly above 0.5 in this interval. 
p <- 0.4
size <- 20
1-pbinom(10, size, p)


## The probability of winning the lottery is 1 in 175,223,510. 
## If 189,000,000 randomly generated (with replacement) tickets are sold, 
## what is the probability that at least one winning tickets is sold? (give your answer as a proportion not percentage)

p <- 1/175223510
size <- 189000000
1-pbinom(0, size, p)

## Using the information from the previous question, what is the probability that two or more winning tickets are sold?
1-pbinom(1, size, p)


## the binomial distribution is approximately normal with N is large and p is not too close to 0 or 1
## The genome has 3 billion bases. About 20% are C, 20% are G, 30% are T and 30% are A. 
## Suppose you take a random interval of 20 bases, what is the exact probability that the GC-content (proportion of Gs of Cs) is greater than 0.35 and smaller or equal to 0.45 in this interval? HINT: use the binomial distribution.


p <- 0.4
size <- 20
pbinom(0.45*size, size, p) - pbinom(0.35*size, size, p)


## what is the normal approximation to the probability?
pnorm(0.45*20, mean=8, sd=sqrt(4.8)) - pnorm(0.35*20, mean=8, sd=sqrt(4.8))

## repeat statistical models before, but using an interval of 1000 bases. 
## what is the difference (in absolute value) between the normal approximation and the
## exact probability of the GC-content being greater than 0.35 and lesser or
## equal to 0.45?

p <- 0.4
size <- 1000
bino <- pbinom(0.45*size, size, p) - pbinom(0.35*size, size, p)
norl <- pnorm(0.45*size, mean=p*size, sd=sqrt(p*(1-p)*size)) - pnorm(0.35*size, mean=p*size, sd=sqrt(p*(1-p)*size))
bino - norl


## use possion approximation when p is very small 
## use the rate lamda = Np representing the number of tickets per N that win the lottery
N <- 189000000
p <- 1/175223510
dpois(2, N*p)

## assumptions for possion to work: N is very very large and Np is not 0
1- dpois(0, N*p) - dpois(1, N*p)
