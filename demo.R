rm(list = ls())
setwd("C:/Users/ljhhe/OneDrive - Florida State University/GitHub/mlyrnetwork")
#usethis::use_package("expm")
library("devtools")
devtools::document()
#devtools::load_all()
#?mlyrnetwork::comp_inner_prod

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
mlnet <- samp_ml(theta, N = N, k = k, H = H ,mdim = mdim, seed = seed, mterm = 'SBM', basis_arguments = basis_arguments)
summary_sim(mlnet)

## plot the k-layer multilayer network sampled above
draw_mlnet(mlnet$net,N)


## Sample and estimate the multilayer network sampled above
estimates <- est_ml(NetMat = mlnet$net, N = N, k = k, H = H, mdim = mdim, mterm = mterm,
                     seed = seed, basis_arguments = basis_arguments) 

## Add standard errors ##
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
N <- 50 # number of nodes
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



############################################################
### test of standard error  ################################
############################################################
rm(list = ls())
setwd("C:/Users/ljhhe/OneDrive - Florida State University/GitHub/mlyrnetwork")

library("devtools")
devtools::document()

set.seed(123456)
seed <- sample(1:9999999,1,replace = FALSE)
N <- 50 # number of nodes
## Parameter settings
burnin <- 100
k <- 3  # number of layers, okay to change
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
summary_sim(mlnet)
## Sample and estimate the multilayer network sampled above
estimates <- est_ml(NetMat = mlnet$net, N = N, k = k, H = H, mdim = mdim, mterm = mterm,
                    seed = seed, basis_arguments = basis_arguments) 
summary_est(estimates) 

sim_est_res <- sim_est(theta,N = N, k = k, H = H, mdim = mdim, mterm = mterm,
                    seed = seed,basis_arguments =basis_arguments)

summary_sim(sim_est_res)
summary_est(sim_est_res)


########################
### Sample iid dyads ###
########################

rm(list = ls())
N <- 71
iid_dyad <- matrix(sample(N, 2*0.2*choose(N,2), replace = TRUE),ncol = 2)
for(i in 1:length(iid_dyad[,1])){
  if(iid_dyad[i,1] > iid_dyad[i,2]){
    t = iid_dyad[i,1]
    iid_dyad[i,1] = iid_dyad[i,2]
    iid_dyad[i,2] <- t
  }

}
iid_dyad <- iid_dyad[-which(iid_dyad[,1] == iid_dyad[,2]),]
iid_dyad <- unique(iid_dyad[,1:2])
iid_dyad <- cbind(iid_dyad,sample(3,length(iid_dyad[,1]),replace = TRUE))
for(i in 1:length(iid_dyad[,1])){
  if(runif(1) < 0.3){
    e <- sample(2,1)
    a <- c(1,2,3)
    a <- a[-which(a == iid_dyad[i,3])]
    iter<-0
    e2 <- sample(a,e)
    for(j in e2){
      iter <- iter + 1
      if(iter > 2){print(iter)}
      iid_dyad <- rbind(iid_dyad,c(iid_dyad[i,1],iid_dyad[i,2],j))
    }
    
  }
}
iid_dyad <- iid_dyad[order(iid_dyad[,1],iid_dyad[,2]),]
dyad_ip <- comp_inner_prod(iid_dyad)
x1 <- dyad_ip$neighboring_inner_prod
x2 <- dyad_ip$non_neighboring_inner_prod
ks.test(x1,x2)
plot(ecdf(x1), 
     xlim = range(c(x1, x2)), 
     col = "blue")
plot(ecdf(x2), 
     add = TRUE, 
     lty = "dashed",
     col = "red")

lazega_ip <- comp_inner_prod()
y1 <- lazega_ip$neighboring_inner_prod
y2 <- lazega_ip$non_neighboring_inner_prod
ks.test(x2,y2)
plot(ecdf(x2), 
     xlim = range(c(x2, y2)), 
     col = "blue")
plot(ecdf(y2), 
     add = TRUE, 
     lty = "dashed",
     col = "red")


#######################################
### Independence test for each dyad ###
#######################################
B <- 1000
x <- comp_indep_test_network(B = B)
quant <- which(sort(x) == x[1])
min_pos <- quant[1]
max_pos <- quant[length(quant)]
#hist(x)
data_frame <- data.frame(x)
ggplot( data_frame, aes( x= x ) ) + 
  geom_histogram( bins=ceiling(range(x)[2]-range(x)[1]+1) )+
  geom_vline(xintercept = x[1],col = "red")
print(min_pos/B)
print(max_pos/B)




