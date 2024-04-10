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
    index_mat = utils::combn(1:k,i)
    for(j in 1:length(index_mat[1,])){
      iter <- iter + 1
      names[iter] <- paste('theta', paste0(index_mat[,j],collapse = ''))
    }
  }
  
  colnames(df) <- names
  table <- insight::export_table(insight::format_table(df))
  return(table)
  
  
  
}