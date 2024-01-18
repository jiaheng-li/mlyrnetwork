
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
er.simulate <- function(theta,N,samp_num = 1,burnin = 100,k = 3,mdim = 6,
                       mterm = rep("ml_order2",mdim),intv = 3,
                       seed = 0, gy = 1){
  
  
  tic <- Sys.time()
  RNETWORK <- rcpp_simulate_ml(samp_num, burnin, intv, mdim, mterm, N, k, theta,seed,gy)
  NetMat <- RNETWORK$elist
  suff_stats <-RNETWORK$suff_list
  basis <- RNETWORK$basis_list
  
  
  
  toc <- Sys.time()
  
  # Prepare result list 
  res <- list(theta = theta, net = NetMat, suff_stats = suff_stats, basis = basis,
              time = as.numeric(difftime(toc, tic, units = "secs")), seed = seed, net_size = N, gy = gy)
  
  
  
  return(res)
  
}



sbm.simulate <- function(theta,N,samp_num = 1,burnin = 100,k = 3,mdim = 6,
                        mterm = rep("ml_order2",mdim),intv = 3,
                        seed = 0, gy = 1,block_num = 5, p_within = 0.5, p_between = 0.05){
  
  
  tic <- Sys.time()
  RNETWORK <- rcpp_simulate_ml_SBM(samp_num, burnin, intv, mdim, mterm, N, k, theta,seed,gy,block_num,p_within,p_between)
  NetMat <- RNETWORK$elist
  suff_stats <-RNETWORK$suff_list
  basis <- RNETWORK$basis_list
  
  
  
  toc <- Sys.time()
  
  # Prepare result list 
  res <- list(theta = theta, net = NetMat, suff_stats = suff_stats, basis = basis,
              time = as.numeric(difftime(toc, tic, units = "secs")), seed = seed, net_size = N, block = block_num, p_w = p_within, p_b = p_between)
  
  
  
  return(res)
  
}


lsm.simulate <- function(theta,N,samp_num = 1,burnin = 100,k = 3,mdim = 6,
                         mterm = rep("ml_order2",mdim),intv = 3,
                         seed = 0, gy = 1, fixed_effect = 0.4, mu = 0, std = 1){
  
  
  tic <- Sys.time()
  RNETWORK <- rcpp_simulate_ml_LSM(samp_num, burnin, intv, mdim, mterm, N, k, theta,seed,gy,fixed_effect, mu, std)
  NetMat <- RNETWORK$elist
  suff_stats <-RNETWORK$suff_list
  basis <- RNETWORK$basis_list
  
  
  
  toc <- Sys.time()
  
  # Prepare result list 
  res <- list(theta = theta, net = NetMat, suff_stats = suff_stats, basis = basis,
              time = as.numeric(difftime(toc, tic, units = "secs")), seed = seed, net_size = N, fe = fixed_effect, mu = mu, std = std)
  
  
  
  return(res)
  
}


sim_sample_save <- function(){
  
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
sim_est <- function(theta,N,samp_num,burnin,k,mdim,mterm,intv,iter_max,
                            seed, gy){
  
  
  tic <- Sys.time()
  RNETWORK <- rcpp_simulate_ml(samp_num, burnin, intv, mdim, mterm, N, k, theta,seed,gy)
  NetMat <- RNETWORK$elist
  result <- rcpp_estimate_model_ml(NetMat, iter_max, samp_num,  burnin, intv, mdim, mterm, N, k, TRUE, seed,gy) 
  dist_mat <- data.frame(result)
  values <- rep(0,mdim)
  for (i in c(1:mdim)){
    values[i] <- dist_mat[i,1]
  }
  
  
  toc <- Sys.time()
  
  # Prepare result list 
  res <- list(theta = theta,theta_est = values, 
              time = as.numeric(difftime(toc, tic, units = "secs")), seed = seed, net_size = N, gy = gy)
  
  
  
  return(res)
  
}





sim_est_save <- function(res_list){
  N <- res_list$net_size
  theta <- res_list$theta
  theta_est <- res_list$theta_est
  time <- res_list$time 
  seed <- res_list$seed
  gy <- res_list$gy
  theta_est <- res_list$theta_est
  
  file_ <- paste0("out/out_", ".out")
  sink(file_)
  
  cat(paste0("\nNetwork size: ", N, "\n"))
  cat(paste0("theta: ", theta, ","))
  cat("\n\n")
  cat(paste0("theta_est: ", theta_est, ","))
  cat("\n\n")
  cat(paste0("\n\nTime elapsed: ", round(time/60), " minutes. Seed: ", seed, " \nP(Y_ij = 1) = ",gy ))
  cat("\n\n")

  file_ <- paste0("results/res_", ".rda")
  save(res_list, file = file_)
  sink()
  closeAllConnections()
  
}

sim_est_wrapper <- function(){
  
}

