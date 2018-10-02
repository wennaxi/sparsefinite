
###Example Code from Online
# summary statistics of sample
n    <- 30
ybar <- 15
s2   <- 3

# sample from the joint posterior (mu, tau | data)
mu     <- rep(NA, 11000)
tau    <- rep(NA, 11000)
T      <- 1000    # burnin
tau[1] <- 1  # initialisation
for(i in 2:11000) {   
  mu[i]  <- rnorm(n = 1, mean = ybar, sd = sqrt(1 / (n * tau[i - 1])))    
  tau[i] <- rgamma(n = 1, shape = n / 2, scale = 2 / ((n - 1) * s2 + n * (mu[i] - ybar)^2))
}
mu  <- mu[-(1:T)]   # remove burnin
tau <- tau[-(1:T)] # remove burnin
hist(mu)
hist(tau)





###My Actual Code Starts Here
for(i in 2:11000) {  
  #a
  for (k in 1:K) {
    ifelse(N[k] != 0, , runif(n = 1, min = 0, max = 1))
  }
  #b
  for (k in 1:(K - 1)) {
    nu[k, i]  = rbeta(n = 1, shape1 = e0 + N[k], shape2 = (K - k) * e0 + sum(N[(k + 1):K]))   
  }
  nu[K, i] = 1
  tau[i] = rgamma(n = 1, shape = n / 2, scale = 2 / ((n - 1) * s2 + n * (mu[i] - ybar)^2))
}