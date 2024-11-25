


#' Title
#'
#' @param net_data 
#' @param N 
#'
#' @return
#' @export
#'
#' @examples
compute_suffstats_Lazega <- function(net_data = Lazega_lawyer_network, N = 71){
  work_sim <- matrix(0,N,N)
  adv_sim <- matrix(0,N,N)
  friend_sim <- matrix(0,N,N)
  for(r in c(1:length(net_data[,1]))){
    row <- net_data[r,1]
    col <- net_data[r,2]
    layer <- net_data[r,3]
    
    if(layer == 1){
      work_sim[row,col] <- 1
    }
    if(layer == 2){
      adv_sim[row,col] <- 1
    }
    if(layer == 3){
      friend_sim[row,col] <- 1
    }
  }
  
  
  suff_mat <- matrix(0,choose(N,2),7) ## matrix of sufficient statistics for each dyad: number of edges \times dimension
  for(i in c(1:(N-1))){
    for(j in c((i + 1):N)){
      m <- ((N-1)+(N-i+1))*(i-1)/2 + j-i
      if (work_sim[i,j] == 1){
        suff_mat[m, 1] <- 1 
        if(adv_sim[i,j] == 1){
          suff_mat[m, 4] <- 1 
          if(friend_sim[i,j] == 1){
            suff_mat[m, 7] <- 1 
          }
        }
        
      }
      
      
      if (adv_sim[i,j] == 1){
        suff_mat[m, 2] <- 1
        if (friend_sim[i,j] == 1 ){
          suff_mat[m, 6] <- 1
        }
        
      }
      if (friend_sim[i,j] == 1 ){
        suff_mat[m, 3] <- 1
        if (work_sim[i,j] == 1){
          suff_mat[m, 5] <- 1 
        }
        
      }
      
    }
  }
  suff_mat <- suff_mat[rowSums(suff_mat[])>0,]
  return(colSums(suff_mat))

  
}



#' Title
#'
#' @param m 
#' @param theta 
#'
#' @return
#' @export
#'
#' @examples
simulate_suffstats_Lazega <- function(data = Lazega_lawyer_network, m = 10, theta = c(-1.4497737, -3.3336582, -2.6945731,  1.8011888,  0.2176545,  2.4582079)){
  # the default theta is the MPLE of the Lazega network data
  reproduced_suff <- matrix(0,m,7)
  seeds <- sample(1:9999999,m,replace = FALSE)
  sim_suff <- matrix(0,m, 7)
  unique_dyads <- unique(data[,1:2])
  unique_dyads_vec <- as.vector(t(unique_dyads))

  for(iter in c(1:m)){
    seed <- seeds[iter]
    reproduced <- samp_ml(theta, N = 71, mdim = 6,
                          mterm = 'other',intv = 3, H = 2, seed = seeds[iter],
                          basis_arguments = unique_dyads_vec)
    
    simsamp <- reproduced$suff_stats
    reproduced_suff[iter,] <- simsamp[,3]
    
  
    
  }
  return(reproduced_suff) 
}


#' Title
#'
#' @param reproduced_suff 
#' @param obs 
#'
#' @return
#' @export
#'
#' @examples
draw_box_plot_Lazega <- function(reproduced_suff, obs){
  suff_y <- c(reproduced_suff)
  suff_mean <- colMeans(reproduced_suff)
  m = length(reproduced_suff[,1])
  dim_x <- c("Coworker", "Advice", "Friendship", "CxA","CxF","AxF","CxAxF")
  df <- data.frame(suff_y,name = as.factor(rep(dim_x,each = m)))
  df$name <- factor(df$name,levels = dim_x)
  
  fig <- ggplot2::ggplot(df,ggplot2::aes(x = name, y = suff_y)) +
    ggplot2::geom_boxplot(fill = "grey") + 
    ggplot2::geom_point( size = 3,ggplot2::aes(x = name,y=rep(obs,each = m)),color="red") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = expression(paste("Box-plot of the reproduced sufficient statistic")) , x = "Layer interaction", y ="Sufficient statistic") + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size = 15),
          axis.text = ggplot2::element_text(size=10),
          axis.title = ggplot2::element_text(size=15),
          panel.border = ggplot2::element_rect(colour = "black", fill=NA,linewidth = 1),
          #legend.title = element_text( size = 30),
          #legend.text = element_text( size = 30),
          #legend.position = "right",
          #legend.key.width = unit(5, 'cm'),
          #legend.background = element_blank(),
          legend.box.background = ggplot2::element_rect(colour = "black",linewidth = 1)) +
    ggplot2::scale_color_manual(values = c("red")) 
  
  return(fig)
  
  
}

#' Compute the inner product of activated dyads.
#' 
#' This function computes inner product of activated dyads.
#'
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
comp_inner_prod <- function(data = Lazega_lawyer_network){
  k <- max(data[,3])
  N <- max(data)
  unique_dyads <- unique(data[,1:2])
  n2 <- length(unique_dyads[,1])
  dyad_vec <- matrix(0,n2,k+2)
  n1 <- length(data[,1])
  ind1 <- data[1,1]
  ind2 <- data[1,2]
  dyad_ind <- 1
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


#' Compute the independence test statistic
#' 
#' Construct the test statistic by computing each dyad's inner product of its neighboring dyads
#' 
#'
#' @param data 
#'
#' @return a list of inner products
#' @export
#'
#' @examples
comp_indep_test <- function(data = Lazega_lawyer_network, B = 100){
  k <- max(data[,3])
  N <- max(data)
  unique_dyads <- unique(data[,1:2])
  n2 <- length(unique_dyads[,1])
  dyad_vec <- matrix(0,n2,k+2)
  n1 <- length(data[,1])
  ind1 <- data[1,1]
  ind2 <- data[1,2]
  dyad_ind <- 1
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
      if(j == i) {next}
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


