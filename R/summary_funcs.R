#' Title
#'
#' @param res_list 
#'
#' @return
#' @export
#'
#' @examples
summary_est <- function(res_list){
  theta_est <- res_list$theta_est
  dim = res_list$mod_dim
  K = res_list$num_layers
  H = res_list$interaction_order
  theta_len = length(theta_est)
  df <- data.frame(matrix(theta_est, ncol = theta_len))
  names <- rep('na',dim)
  iter <- 0 
  for(i in 1:H){
    index_mat = utils::combn(1:K,i)
    for(j in 1:length(index_mat[1,])){
      iter <- iter + 1
      names[iter] <- paste('theta', paste0(index_mat[,j],collapse = ''))
    }
  }
  
  colnames(df) <- names
  table <- insight::export_table(insight::format_table(df))
  return(table)
  
  
  
}


#' Title
#'
#' @param res_list 
#'
#' @return
#' @export
#'
#' @examples
summary_sim <- function(res_list){
  suff_mat <- res_list$suff_stats
  df <- data.frame(t(round(suff_mat[,3])))
  names <- c('s1','s2','s3','s12','s13','s23','s123')

  colnames(df) <- names
  table <- insight::export_table(insight::format_table(df))
  return(table)
  
  
  
}



# Function to generate unique pairs
#' Title
#'
#' @param pool 
#' @param num_pairs 
#'
#' @return
#' @export
#'
#' @examples
generate_basis <- function(n, pool_size = 100) {
  pairs <- list()
  while (length(pairs) < n) {
    pair <- sort(sample(1:pool_size, 2))
    if (!any(sapply(pairs, function(x) all(x == pair)))) {
      pairs <- c(pairs, list(pair))
    }
  }
  return(as.vector(t(do.call(rbind, pairs))))
}
