
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
  dist_mat <- data.frame(result)
  values <- rep(0,mdim)
  for (i in c(1:mdim)){
    values[i] <- dist_mat[i,1]
  }
  
  
  toc <- Sys.time()
  
  # Prepare result list 
  res <- list(theta = theta,theta_est = values, 
              time = as.numeric(difftime(toc, tic, units = "secs")), seed = seed, net_size = N)
  
  
  
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




