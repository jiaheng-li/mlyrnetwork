#' Title
#'
#' @param N 
#' @param node_membs 
#' @param p_within 
#' @param p_between 
#'
#' @return
#' @export
#'
#' @examples
gen_sbm <- function(N, node_membs, p_within, p_between) { 
  net <- matrix(0, nrow = N, ncol = N) 
  for (i in 1:(N-1)) { 
    for (j in (i+1):N) {
      net[i, j] <- rbinom(1, 1, ifelse(node_membs[i] == node_membs[j], p_within, p_between))
      net[j, i] <- net[i, j]
    }
  }
  #net <- as.network(net, directed = FALSE)
  return(net)
}

#' Title
#'
#' @param N 
#' @param num_of_block 
#' @param p_within 
#' @param p_between 
#' @param M 
#'
#' @return
#' @export
#'
#' @examples
generate_SBM_basis <- function(N=3*10, num_of_block = 3, p_within=0.5, p_between=0.1,M=0){
  node_membs <- rep(c(1:num_of_block),each = N/num_of_block)
  net <- gen_sbm(N, node_membs, p_within, p_between)
  dyads = matrix(0,1,2)
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      if(net[i,j] == 1){
        dyads <- rbind(dyads,c(i,j))
      }
    }
  }
  dyads <- dyads[-1,]
  
  # if we need to sample a subgraph of M activated dyads
  if(M){
    num_of_dyad <- length(dyads[,1])
    if(num_of_dyad < M){
      stop("The number of activated dyads required is more than the number of dyads in the network.\n
         Try reducing the number of activated dyads")
    }
    dyads <- dyads[sample(1:num_of_dyad,M,replace = FALSE),]
  }
  
  
  return(as.vector(t(dyads)))
}




# Function to generate Bernoulli random graphs
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