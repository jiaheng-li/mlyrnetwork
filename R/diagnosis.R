#=================================================================
# This file contains functions for conditional independence test
# and link tracing sampling of separable multilayer networks.
#            
#=================================================================

#' Title
#'
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
comp_dyad_mat <- function(data = Lazega_lawyer_network){
  k <- max(data[,3])
  N <- max(data)
  unique_dyads <- unique(data[,1:2])
  n2 <- length(unique_dyads[,1])
  dyad_vec <- matrix(0,n2,k+2)
  n1 <- length(data[,1])
  ind1 <- data[1,1]
  ind2 <- data[1,2]
  dyad_ind <- 1
  # Form a matrix where each row is a dyad vector
  for(i in 1:n1){
    if(ind1 == data[i,1] & ind2 == data[i,2]){
      dyad_vec[dyad_ind,1] <- ind1
      dyad_vec[dyad_ind,2] <- ind2
      dyad_vec[dyad_ind,data[i,3]+2] <- 1
    }
    else{
      dyad_ind <- dyad_ind + 1
      ind1 <- data[i,1]
      ind2 <- data[i,2]
      dyad_vec[dyad_ind,1] <- ind1
      dyad_vec[dyad_ind,2] <- ind2
      #dyad_ind <- data[i,2] - data[i,1] + (N-1 + (N-data[i,1]+1)* (data[i,1]!=2)) * (data[i,1]-1) / (2/((data[i,1]==2)+1))  * (data[i,1]!=1)
      dyad_vec[dyad_ind,data[i,3]+2] <- 1
    }
  }
  return(dyad_vec)
}


#' Compute the independence test statistic for the whole network
#' 
#' Construct the network-wise test statistic by summing over each dyad's inner product of its neighboring dyads
#'
#' @param data n by 3 matrix with each column indicating index1, index2, and layer. The dyad index must be lexicographically ordered.
#' @param B 
#'
#' @return
#' @export
#'
#' @examples
comp_indep_test_network <- function(data = Lazega_lawyer_network, B = 100){
  dyad_vec <- comp_dyad_mat(data)
  n2 <- length(dyad_vec[,1])
  
  
  TC <- rep(0,B+1) # list of test statistics
  
  ## test statistic from the observation
  neighboring_inner_prod <- rep(0,n2)
  # sum of neighboring inner products for each activated dyads
  for(i in 1:n2){
    node1 <- dyad_vec[i,1]
    node2 <- dyad_vec[i,2]
    #non_neighboring_inner_prod <- rep(0,B)
    neib_ind <- which(dyad_vec[,1] == node1 | dyad_vec[,1] == node2 | dyad_vec[,2] == node1 | dyad_vec[,2] == node2)
    #non_neib_ind <- which(!(dyad_vec[,1] == node1 | dyad_vec[,1] == node2 | dyad_vec[,2] == node1 | dyad_vec[,2] == node2))
    for(j in neib_ind){
      if(j == i){next}
      neighboring_inner_prod[i] <- neighboring_inner_prod[i] + sum(dyad_vec[i,3:5] * dyad_vec[j,3:5])
    }
  }
  TC[1] <- mean(neighboring_inner_prod)
  
  ## test statistic calculated by permutation
  for(b in 1:B){
    neighboring_inner_prod <- rep(0,n2)
    permuted_dyad <- cbind(dyad_vec[,1:2],dyad_vec[sample(nrow(dyad_vec)),3:5])
    # sum of neighboring inner products for each activated dyads
    for(i in 1:n2){
      node1 <- permuted_dyad[i,1]
      node2 <- permuted_dyad[i,2]
      neib_ind <- which(permuted_dyad[,1] == node1 | permuted_dyad[,1] == node2 | permuted_dyad[,2] == node1 | permuted_dyad[,2] == node2)
      for(j in neib_ind){
        if(j == i){next}
        neighboring_inner_prod[i] <- neighboring_inner_prod[i] + sum(permuted_dyad[i,3:5] * permuted_dyad[j,3:5])
      }
    }
    TC[b+1] <- mean(neighboring_inner_prod)
  }
  
  
  return(TC)
  
}


#' Title
#'
#' @param data 
#' @param init number of initial nodes
#' @param w number of waves
#' @param N 
#'
#' @return
#' @export
#'
#' @examples
link_tracing_samp <- function(data = Lazega_lawyer_network,init = 3, w = 2, N = 0){
  dyad_vec <- comp_dyad_mat(data)
  if(!N){
    N <- max(data[,1:2])
  }
  wave <- list()
  n2 <- length(dyad_vec[,1])
  wave[[1]] <- sample(N,init)
  for(j in 1:w){
    dyad_ind <- c()
    for(i in wave[[j]]){
      dyad_ind <- append(dyad_ind, which(dyad_vec[,1] == i | dyad_vec[,2] == i))
    }
    wave[[j+1]] <- setdiff(unique(c(dyad_vec[dyad_ind,1:2])),unlist(wave))
  }
  node_set <- unlist(wave)
  lt_samp <- matrix(0,1,3)
  for(i in 1:n2){
    if(dyad_vec[i,1] %in% node_set & dyad_vec[i,2] %in% node_set){
      for(j in 3:length(dyad_vec[1,])){
        if(dyad_vec[i,j] == 1){
          lt_samp <- rbind(lt_samp,c(dyad_vec[i,1],dyad_vec[i,2],j-2))
        }
        
      }
    }
  }
  lt_samp <- lt_samp[-1,]
  return(lt_samp)
}


#' Compute the independence test statistic for each activated dyad
#' 
#' Construct the test statistic by computing each dyad's inner product of its neighboring dyads
#' 
#'
#' @param data 
#'
#' @return a list of inner products. Each element in the list corresponds to a dyad.
#'  The first element of each list is the inner product of all neighboring dyads, and the rest are inner products after random permutation.
#' @export
#'
#' @examples
comp_indep_test_dyad <- function(data = Lazega_lawyer_network, B = 100){
  k <- max(data[,3])
  N <- max(data)
  unique_dyads <- unique(data[,1:2])
  n2 <- length(unique_dyads[,1])
  dyad_vec <- matrix(0,n2,k+2)
  n1 <- length(data[,1])
  ind1 <- data[1,1]
  ind2 <- data[1,2]
  dyad_ind <- 1
  # Form a matrix where each row is a dyad vector
  for(i in 1:n1){
    if(ind1 == data[i,1] & ind2 == data[i,2]){
      dyad_vec[dyad_ind,1] <- ind1
      dyad_vec[dyad_ind,2] <- ind2
      dyad_vec[dyad_ind,data[i,3]+2] <- 1
    }
    else{
      dyad_ind <- dyad_ind + 1
      ind1 <- data[i,1]
      ind2 <- data[i,2]
      dyad_vec[dyad_ind,1] <- ind1
      dyad_vec[dyad_ind,2] <- ind2
      #dyad_ind <- data[i,2] - data[i,1] + (N-1 + (N-data[i,1]+1)* (data[i,1]!=2)) * (data[i,1]-1) / (2/((data[i,1]==2)+1))  * (data[i,1]!=1)
      dyad_vec[dyad_ind,data[i,3]+2] <- 1
    }
    
  }
  
  ip_dist <- list()
  for(i in 1:n2){
    node1 <- dyad_vec[i,1]
    node2 <- dyad_vec[i,2]
    neighboring_inner_prod <- 0
    non_neighboring_inner_prod <- rep(0,B)
    neib_ind <- which(dyad_vec[,1] == node1 | dyad_vec[,1] == node2 | dyad_vec[,2] == node1 | dyad_vec[,2] == node2)
    non_neib_ind <- which(!(dyad_vec[,1] == node1 | dyad_vec[,1] == node2 | dyad_vec[,2] == node1 | dyad_vec[,2] == node2))
    for(j in neib_ind){
      if(j == i){next}
      neighboring_inner_prod <- neighboring_inner_prod + sum(dyad_vec[i,3:5] * dyad_vec[j,3:5])
    }
    n3 <- length(neib_ind)
    for(b in 1:B){
      n4 <- sample(n3, 1)
      replaced_ind <- sample(neib_ind,n4)
      replacing_ind <- sample(non_neib_ind, n4)
      new_ind <- append(neib_ind[-replaced_ind],replacing_ind)
      for(j in new_ind){
        non_neighboring_inner_prod[b] <- non_neighboring_inner_prod[b] + sum(dyad_vec[i,3:5] * dyad_vec[j,3:5])
      }
    }
    
    ip_dist[[i]] <- c(neighboring_inner_prod, non_neighboring_inner_prod)
  }
  
  
  res <- list(ip_dist = ip_dist)
  
  
  class(res) <- 'inner_prod_list'
  return(res)
  
}



#' Compute the average neighboring and non-neighboring inner product of each activated dyad.
#' 
#' This function computes the neighboring and non-neighboring inner product mean of each activated dyad.
#'
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
comp_inner_prod_avg <- function(data = Lazega_lawyer_network){
  k <- max(data[,3])
  N <- max(data)
  unique_dyads <- unique(data[,1:2])
  n2 <- length(unique_dyads[,1])
  dyad_vec <- matrix(0,n2,k+2)
  n1 <- length(data[,1])
  ind1 <- data[1,1]
  ind2 <- data[1,2]
  dyad_ind <- 1
  # Form a matrix where each row is a dyad vector
  for(i in 1:n1){
    if(ind1 == data[i,1] & ind2 == data[i,2]){
      dyad_vec[dyad_ind,1] <- ind1
      dyad_vec[dyad_ind,2] <- ind2
      dyad_vec[dyad_ind,data[i,3]+2] <- 1
    }
    else{
      dyad_ind <- dyad_ind + 1
      ind1 <- data[i,1]
      ind2 <- data[i,2]
      dyad_vec[dyad_ind,1] <- ind1
      dyad_vec[dyad_ind,2] <- ind2
      #dyad_ind <- data[i,2] - data[i,1] + (N-1 + (N-data[i,1]+1)* (data[i,1]!=2)) * (data[i,1]-1) / (2/((data[i,1]==2)+1))  * (data[i,1]!=1)
      dyad_vec[dyad_ind,data[i,3]+2] <- 1
    }
    
  }
  neighboring_inner_prod <- rep(0,n2)
  non_neighboring_inner_prod <- rep(0,n2)
  for(i in 1:n2){
    node1 <- dyad_vec[i,1]
    node2 <- dyad_vec[i,2]
    iter1 <- 0
    iter2 <- 0
    for(j in 1:n2){
      if(j==i){next}
      else if(dyad_vec[j,1] == node1 | dyad_vec[j,1] == node2 | dyad_vec[j,2] == node1 | dyad_vec[j,2] == node2){
        neighboring_inner_prod[i] <- neighboring_inner_prod[i] + sum(dyad_vec[i,3:5] * dyad_vec[j,3:5])
        iter1 <- iter1 + 1
      }
      else{
        non_neighboring_inner_prod[i] <- non_neighboring_inner_prod[i] + sum(dyad_vec[i,3:5] * dyad_vec[j,3:5])
        iter2 <- iter2 + 1
      }
      
    }
    neighboring_inner_prod[i] <- neighboring_inner_prod[i] / iter1
    non_neighboring_inner_prod[i] <- non_neighboring_inner_prod[i] / iter2
  }
  
  
  res <- list(neighboring_inner_prod = neighboring_inner_prod, non_neighboring_inner_prod = non_neighboring_inner_prod)
  
  
  class(res) <- 'inner_prod_class'
  return(res)
  
}





