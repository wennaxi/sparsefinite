



rm(list=ls())


#Parameters for MCMC
nburn = 1000
niter = 11000


#Parameters for model
N.sum = 1000
K = 10



#Parameters in model theta[k] S[i] N[k] nu[k] eta[k] e0
theta = matrix(NA, nrow = niter, ncol = K)
S = matrix(NA, nrow = niter, ncol = N.sum)
N = matrix(NA, nrow = niter, ncol = K)
nu = matrix(NA, nrow = niter, ncol = K)
eta = matrix(NA, nrow = niter, ncol = K)
e0 = rep(NA, niter)
ae = 1
be = 10
#Initial Values
theta[1, ] = rep(.5, K)
S[1, ] = c(sample(1:K, N.sum, replace = T))
for (k in 1:K) {
  N[1, k] = sum(S[1, ] == k)
}
for (k in 1:K) {
  nu[1, k] = 1 / (1 + K - k) #inital value set to mean of its prior
  eta[1, k] = ifelse(k == 1, nu[1, k], nu[1, k] * prod(1 - nu[1, 1:(k - 1)])) #transfrom back to eta's
}
e0[1] = ae / be #inital value set to mean of its prior


#Data y[i]
y = rep(NA, N.sum)
y[1:200] = rbinom(n = 200, size = 1, p = .1)
y[201:400] = rbinom(n = 200, size = 1, p = .3)
y[401:600] = rbinom(n = 200, size = 1, p = .5)
y[601:800] = rbinom(n = 200, size = 1, p = .7)
y[801:1000] = rbinom(n = 200, size = 1, p = .9)


###MCMC
for(i in 2:niter) {  
  #a
  for (k in 1:K) {
    theta[i, k] = ifelse(N[i - 1, k] != 0, 
                           rbeta(n = 1, shape1 = sum((S[i - 1, ] == k) * y) + 2, shape2 = sum((S[i - 1, ] == k) * (1 - y)) + 1), 
                           runif(n = 1, min = 0, max = 1))
  }
  #b
  for (k in 1:(K - 1)) {
    nu[i, k]  = rbeta(n = 1, shape1 = e0 + N[i - 1, k], shape2 = (K - k) * e0 + sum(N[i - 1, (k + 1):K]))   
    eta[i, k] = ifelse(k == 1, nu[i, k], nu[i, k] * prod(1 - nu[i, 1:(k - 1)])) #transfrom back to eta's
  }
  nu[i, K] = 1
  eta[i, K] = nu[i, K] * prod(1 - nu[i, 1:(K - 1)])
  #c
  for (n in 1:N.sum) {
    pr = rep(NA, K)
    for (k in 1:K) {
      pr[k] = eta[i, k] * theta[i, k]^y[n] * (1 - theta[i, k])^(1 - y[n])
    }
    S[i, n] = which(rmultinom(n = 1, size = 1, prob = pr) == 1)
  }
  for (k in 1:K) {
    N[i, k] = sum(S[i, ] == k) #update N[[i]]
  }
  #d
  e0.proposed = rnorm(n = 1, mean = e0[i - 1], sd = 1)
  if (e0.proposed <= 0) {
    e0[i] = e0[i - 1]
  } else {
    log_f = function(e0_value) {
      # partition_conditional = gamma(K * e0_value) / gamma(N.sum + K * e0_value) * prod(gamma(N[i, ] + e0_value) / gamma(e0_value))
      # e0_prior = e0_value^(ae - 1) * exp(-be * e0_value)
      # return(partition_conditional * e0_prior)
      log_partition_conditional = lgamma(K * e0_value) - lgamma(N.sum + K * e0_value) + sum(lgamma(N[i, ] + e0_value) - lgamma(e0_value))
      log_e0_prior = (ae - 1) * log(e0_value) - be * e0_value
      return(log_partition_conditional + log_e0_prior)
    }
    #alpha = f(e0.proposed) / f(e0[i - 1])
    alpha = exp(log_f(e0.proposed) - log_f(e0[i - 1]))
    u = runif(n = 1, min = 0, max = 1)
    e0[i] = ifelse(u <= alpha, e0.proposed, e0[i - 1])
  }
}

dim(theta)
dim(S)
dim(N)
dim(nu)
dim(eta)
length(e0)


# remove burnin theta[k] S[i] N[k] nu[k] eta[k] e0
theta = theta[-(1:nburn), ]   
S = S[-(1:nburn), ] 
N = N[-(1:nburn), ]
nu = nu[-(1:nburn), ]
eta = eta[-(1:nburn), ]
e0 = e0[-(1:nburn)]




dim(theta)
dim(S)
dim(N)
dim(nu)
dim(eta)
length(e0)

hist(theta[, 1])
for (k in 1:K) {
  hist(theta[, k])
}



colMeans(theta)
#colMeans(S)
colMeans(N)
colMeans(nu)
colMeans(eta)
mean(e0)


dim(N)
N[9901:10000, ]
