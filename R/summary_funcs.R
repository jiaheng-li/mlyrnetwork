


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
  H = res_list$highest_order
  theta_len = length(theta_est)
  mple_sd <- res_list$mple_sd
  df <- data.frame(matrix(c(theta_est,mple_sd), ncol = theta_len, byrow = TRUE))
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
  suff_stats <- res_list$suff_stats
  dim = res_list$mod_dim
  K = res_list$num_layers
  H = res_list$highest_order
  suff_len = length(suff_stats)
  suff_stats_sd <- res_list$suff_stats_sd
    df <- data.frame(matrix(c(suff_stats,suff_stats_sd), ncol = suff_len, byrow = TRUE))
  names <- rep('na',dim)
  iter <- 0 
  for(i in 1:H){
    index_mat = utils::combn(1:K,i)
    for(j in 1:length(index_mat[1,])){
      iter <- iter + 1
      names[iter] <- paste('S', paste0(index_mat[,j],collapse = ''))
    }
  }
  
  colnames(df) <- names
  table <- insight::export_table(insight::format_table(df))
  return(table)
  
  
  
  
}

