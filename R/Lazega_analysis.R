
est_data <- function(net_data, samp_num = 1,  burnin = 100, intv = 3, mdim = 6, mterm = 'other', N , k = 3, H = 2, seed = 0, iter_max = 30000 ){
  unique_dyads <- unique(net_data[,1:2])
  unique_dyads_vec <- as.vector(t(unique_dyads))
  
  tic <- Sys.time()
  result <- rcpp_estimate_model_ml_Hway(net_data, samp_num,  burnin, intv, mdim, mterm, N , k, H, seed,  unique_dyads_vec, iter_max) 
  dist_mat <- data.frame(result)
  values <- rep(0,mdim)
  for (i in c(1:mdim)){
    values[i] <- dist_mat[i,1]
  }
  toc <- Sys.time()
  
  # Prepare result list 
  res <- list( theta_est = values, 
              time = as.numeric(difftime(toc, tic, units = "secs")), seed = seed, net_size = N)
  
  
  
  return(res)
}
