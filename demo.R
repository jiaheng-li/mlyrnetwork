rm(list = ls())

library("devtools")
devtools::document()

set.seed(123456)

## Parameter settings
burnin <- 100
k <- 4  # number of layers
H <- 2  # highest order of interactions
mdim <- 0
for(i in 1:H){
  mdim <- mdim + choose(k,i)
}
mterm <- 'SBM'  # generate stochastic block model basis
intv <- 3
iter_max <- 30000
# SBM parameters
B <- 3
p_within <- .25
p_between <- .01 

# generate theta for layers
params <- matrix(runif(mdim,-1,1),1,mdim)
params[,c(3,6)] = 0
theta <- params[1,]

N <- 100
iter_max <- 30000
seed <- sample(1:9999999,1,replace = FALSE)

# set basis network arguments for SBM basis
basis_arguments <- c(B, p_within ,p_between, seed)

## Sample a k-layer multilayer network
mlnet <- samp_ml(theta,N = N, k = k, H = H ,mdim = mdim, seed = seed, mterm = 'SBM', basis_arguments = basis_arguments)
draw_mlnet(mlnet$net,N)


## Sample and estimate a multilayer network, by default it is 
estimates <- est_ml(NetMat = mlnet$net, theta=theta, N = N, k = k, H = H, mdim = mdim, mterm = mterm,
                     seed = seed, basis_arguments = basis_arguments) 

print(estimates)

