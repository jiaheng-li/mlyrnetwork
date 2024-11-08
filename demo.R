rm(list = ls())
setwd("C:/Users/ljhhe/OneDrive - Florida State University/GitHub/mlyrnetwork")
library("devtools")
devtools::document()


set.seed(123456)
seed <- sample(1:9999999,1,replace = FALSE)

## Parameter settings
burnin <- 100
k <- 4  # number of layers, okay to change
H <- 2  # highest order of interactions, okay to change (<= k)
mdim <- 0
for(i in 1:H){
  mdim <- mdim + choose(k,i)
}

intv <- 3 # sampling iterval
iter_max <- 30000 # maximum iterations for estimation


# generate stochastic block model basis as an example
# can be set to 'Ber', 'LSM' and 'other'
mterm <- 'SBM'  
# set SBM parameters
B <- 3
p_within <- .25
p_between <- .01 
# set basis network arguments for SBM basis
basis_arguments <- c(B, p_within ,p_between, seed)

# generate theta for layers
params <- matrix(runif(mdim,-1,1),1,mdim)
params[,c(3,6)] = 0
theta <- params[1,]

N <- 100 # number of nodes

## Sample a k-layer multilayer network
mlnet <- samp_ml(theta,N = N, k = k, H = H ,mdim = mdim, seed = seed, mterm = 'SBM', basis_arguments = basis_arguments)
summary_sim(mlnet)

## plot the k-layer multilayer network sampled above
draw_mlnet(mlnet$net,N)


## Sample and estimate the multilayer network sampled above
estimates <- est_ml(NetMat = mlnet$net, N = N, k = k, H = H, mdim = mdim, mterm = mterm,
                     seed = seed, basis_arguments = basis_arguments) 

### Add standard errors ###
summary_est(estimates) 

## Examples of invalid input
est_ml(NetMat = mlnet$net, N = N, k = k, H = k+1, mdim = mdim, mterm = mterm,
       seed = seed, basis_arguments = basis_arguments) 


## Examples of data analysis using the Lazega lawyer network data set
obs <- compute_suffstats_Lazega() # compute the observed sufficient statistics from the Lazega network.

# compute the basis network of Lazega
unique_dyads <- unique(Lazega_lawyer_network[,1:2])
unique_dyads_vec <- as.vector(t(unique_dyads))

# estimate ERGM model parameters using MPLE.
data_res <- est_ml(NetMat = Lazega_lawyer_network, N = 71, k = 3, H = 2, mdim = 6, mterm = 'other',
       seed = seed, basis_arguments = unique_dyads_vec)  

# obtain the estimated theta for Lazega network
data_theta <- data_res$theta_est
summary_est(data_res)

# m indicates the number of repetitions
simulated_suff <- simulate_suffstats_Lazega(data = Lazega_lawyer_network, m = 20 , theta = data_theta) # reproduce the Lazega network using the estimated parameters and the basis network induced from the observed network.
draw_box_plot_Lazega(simulated_suff, obs)

############################################################
### test of sufficient statistic functions for each dyad ###
############################################################
rm(list = ls())
setwd("C:/Users/ljhhe/OneDrive - Florida State University/GitHub/mlyrnetwork")

library("devtools")
devtools::document()

set.seed(123456)
seed <- sample(1:9999999,1,replace = FALSE)
N <- 1000 # number of nodes
## Parameter settings
burnin <- 100
k <- 4  # number of layers, okay to change
H <- 2  # highest order of interactions, okay to change (<= k)
mdim <- 0
for(i in 1:H){
  mdim <- mdim + choose(k,i)
}

intv <- 3 # sampling iterval
iter_max <- 30000 # maximum iterations for estimation
dyads <- 1000

# generate stochastic block model basis as an example
# can be set to 'Ber', 'LSM' and 'other'
mterm <- 'BER'  
gy <- dyads / choose(N,2) 
# set basis network arguments for SBM basis
basis_arguments <- c(gy, seed)

# generate theta for layers
params <- matrix(runif(mdim,-1,1),1,mdim)
params[,c(3,6)] = 0
theta <- params[1,]


mlnet <- samp_ml(theta,N = N, k = k, H = H ,mdim = mdim, seed = seed, mterm = mterm, basis_arguments = basis_arguments)

a = rcpp_compute_dyad_suffstats(mlnet$net, 1, 100, intv, mdim, mterm, N, k, H, seed, basis_arguments)
colSums(a$dyad_suffstats)
apply(a$dyad_suffstats,2,sd)
mlnet$suff_stats
