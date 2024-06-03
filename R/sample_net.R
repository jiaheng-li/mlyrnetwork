
#' Title
#'
#' @param theta 
#' @param N 
#' @param samp_num 
#' @param burnin 
#' @param k 
#' @param H 
#' @param mdim 
#' @param mterm 
#' @param intv 
#' @param seed 
#' @param iter_max 
#' @param basis_arguments 
#'
#' @return
#' @export
#'
#' @examples
est_ml <- function(NetMat,N,samp_num = 1,burnin = 100,k = 3, H = 2, mdim, mterm = "BER",intv = 3,
                    seed = 0, iter_max = 30000, basis_arguments = c(1,0)){
  
  if(N <= 1){
    stop("N must be greater than one")
  }
  if(k < H){
    stop("H must not be greater than k")
  }
  tic <- Sys.time()
  result <- rcpp_estimate_model_ml_Hway(NetMat, samp_num,  burnin, intv, mdim, mterm, N, k,H, seed,basis_arguments, iter_max) 
  dist_mat <- data.frame(result)
  values <- rep(0,mdim)
  for (i in c(1:mdim)){
    values[i] <- dist_mat[i,1]
  }
  
  
  toc <- Sys.time()
  
  # Prepare result list 
  res <- list(theta_est = values, 
              time = as.numeric(difftime(toc, tic, units = "secs")), seed = seed, net_size = N, num_layers = k, mod_dim = mdim, interaction_order = H)
  
  
  class(res) <- 'estimate_class'
  return(res)
  
}



#' Title
#'
#' @param theta model-generating parameter
#' @param N network size
#' @param samp_num 
#' @param burnin 
#' @param k 
#' @param mdim 
#' @param mterm 
#' @param intv 
#' @param iter_max 
#' @param seed 
#' @param gy 
#'
#' @return a list
#' @export
#'
#' @examples
sim_est <- function(theta,N,samp_num = 1,burnin = 100,k = 3, H = 2, mdim, mterm = "BER",intv = 3,
                    seed = 0, iter_max = 30000, basis_arguments = c(1,0)){
  
  
  tic <- Sys.time()
  RNETWORK <- rcpp_simulate_ml_Hway(samp_num, burnin, intv, mdim, mterm, N, k, theta, H, seed, basis_arguments)
  NetMat <- RNETWORK$elist
  result <- rcpp_estimate_model_ml_Hway(NetMat, samp_num,  burnin, intv, mdim, mterm, N, k,H, seed,basis_arguments, iter_max) 
  
  
  # calculate covariance matrix of mple by sufficient statistics
  a <- rcpp_compute_dyad_suffstats(NetMat, samp_num, burnin, intv, mdim, mterm, N, k, H, seed,basis_arguments)
  suff_matrix <- a$dyad_suffstats
  df <- data.frame(suff_matrix)
  cond_suff_mat <- df[rowSums(df[,3:mdim]) > 0,]
  mple_cov <- cov(cond_suff_mat[,3:(mdim+2)])
  num_of_dyad <- length(cond_suff_mat[,1])
  
  dist_mat <- data.frame(result)
  values <- rep(0,mdim)
  for (i in c(1:mdim)){
    values[i] <- dist_mat[i,1]
  }
  
  RL2err <- norm(values - theta, "2")/norm(theta,"2")
  toc <- Sys.time()
  
  
  
  # Prepare result list 
  res <- list(theta = theta,theta_est = values, RL2err = RL2err, suff_matrix = suff_matrix,
              time = as.numeric(difftime(toc, tic, units = "secs")), seed = seed, net_size = N, num_of_layers = k, highest_order = H,
              num_of_dyad = num_of_dyad, mple_cov = mple_cov, basis_net = basis_arguments)
  
  
  
  return(res)
  
}





#' Title
#'
#' @param res_list 
#'
#' @return
#' @export
#'
#' @examples
sim_est_save <- function(res_list, sim_id){
  
  file_ <- paste0("results/res_", sim_id, ".rda")
  save(res_list, file = file_)
  sink()
  
  
}


#' Title
#'
#' @param theta 
#' @param N 
#' @param samp_num 
#' @param burnin 
#' @param k 
#' @param H 
#' @param mdim 
#' @param mterm 
#' @param intv 
#' @param seed 
#' @param iter_max 
#' @param basis_arguments 
#'
#' @return
#' @export
#'
#' @examples
sim_est_wrapper <- function(sim_id, thetas, N, samp_num = 1, burnin = 100, k, H, mdim, mterm, intv = 3,
                            seeds, iter_max = 30000, basis_arguments){
  
  file_ <- paste0("out/out_", sim_id, ".out")
  sink(file_)
  theta <- thetas[sim_id,]
  cat(paste0("\nNetwork size: ", N, "\n"))
  #cat(paste0("\ntheta: ", theta, "\n"))
  cat(paste0("\nStarting replication: ", sim_id, "\n"))
  cat(paste0("\nModel term is ", mterm, "\n"))
  cat(paste0("\nThe number of layers is k = ", k, "\n"))
  cat(paste0("\nThe highest order of interaction is H = ",  H, "\n"))
  cat(paste0("\nThe number of activated dyads is: ", length(basis_arguments)/2, "\n"))
  seed <- seeds[sim_id]
  
  res_list <- sim_est(theta, N, samp_num, burnin, k, H, mdim, mterm, intv = 3,
                      seed, iter_max = 30000, basis_arguments)
  
  cat(paste0("\nTime elapsed: ", round(res_list$time/60), " minutes. Seed: ", seed))
  cat(paste0("\nThe relative L2 error is : ", res_list$RL2err))

  sim_est_save(res_list, sim_id)
  return(res_list)
  
  
}


#' Title
#'
#' @param theta model-generating parameter vector by lexicographic order
#' @param N number of nodes
#' @param samp_num number of simulated samples
#' @param burnin number of burn-ins for the MCMC algorithm
#' @param k number of layers
#' @param mdim number of parameter dimensions
#' @param mterm model terms
#' @param intv length of the update interval for the MCMC algorithm
#' @param seed random seed
#' @param gy edge probability for the Bernoulli basis network
#'
#' @return return a list of k-layer network with a Bernoulli basis network and ...
#' @export
#'
#' @examples 
samp_ml <- function(theta,N = 10,samp_num = 1,burnin = 100,k = 3,mdim,
                       mterm = 'BER',intv = 3, H = 2,
                       seed = 0, basis_arguments = c(0.5,0)){
  
  
  tic <- Sys.time()
  RNETWORK <- rcpp_simulate_ml_Hway(samp_num, burnin, intv, mdim, mterm, N, k, theta, H, seed,basis_arguments)
  NetMat <- RNETWORK$elist
  suff_stats <-RNETWORK$suff_list
  basis <- RNETWORK$basis_list
  
  
  
  toc <- Sys.time()
  
  # Prepare result list 
  res <- list(theta = theta, net = NetMat, suff_stats = suff_stats, basis = basis,
              time = as.numeric(difftime(toc, tic, units = "secs")), seed = seed, net_size = N)
  
  
  class(res) <- 'simulate_class'
  return(res)
  
}




sim_sample_save <- function(){
  
}




