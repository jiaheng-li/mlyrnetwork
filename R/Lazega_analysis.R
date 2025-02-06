

#' Title
#'
#' @param adj 
#' @param K 
#'
#' @return
#' @export
#'
#' @examples
transform_adjmatrix <- function(adj, K){
  N <- length(adj[[1]][1,])
  mlnetwork <- matrix(0,1,3)
  
  ### Create symmetric adjacent matrix
  for(i in c(1:(N-1))){
    for(j in c((i + 1):N)){
      for(k in 1:K)
        if(adj[[k]][i,j] == 1 && adj[[k]][j,i] == 1){
          mlnetwork <- rbind(mlnetwork,c(i,j,k))
  
        }
      
    }
  }
  mlnetwork_vec <- mlnetwork[-1,]
  return(mlnetwork_vec)
  
  
}

#' Title
#' Only calculate sufficient stats for models without interactions
#' @param adj 
#' @param H 
#'
#' @return
#' @export
#'
#' @examples
compute_suffstats <- function(adj, H = 1){
  K <- length(adj)
  N <- length(adj[[1]][1,])
  mdim <- 0
  for(h in 1:H){
    mdim <- mdim + choose(K,h)
  }
  
  suff_mat <- matrix(0,choose(N,2),mdim+1) ## matrix of sufficient statistics for each dyad: number of edges \times dimension
  for(i in c(1:(N-1))){
    for(j in c((i + 1):N)){
      m <- ((N-1)+(N-i+1))*(i-1)/2 + j-i
      #iter1 <- 1
      f <- 0
      for(d in 1:mdim){
        # if(d > choose(K,iter)){
        #   p <- d - K
        #   
        # }
        if(adj[[d]][i,j] == 1){
          suff_mat[m,d] <- 1
          f <- f + 1
        }
        if(f == K){suff_mat[m,K+1] <- 1}
      }
      
    }
  }
  suff_mat <- suff_mat[rowSums(suff_mat[])>0,]
  return(colSums(suff_mat))
}


#' Title
#'
#' @param net_data 
#' @param N 
#'
#' @return
#' @export
#'
#' @examples
compute_suffstats_Lazega <- function(net_data = Lazega_lawyer_network,mdim = 6, N = 71){
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
  
  
  suff_mat <- matrix(0,choose(N,2),mdim) ## matrix of sufficient statistics for each dyad: number of edges \times dimension
  for(i in c(1:(N-1))){
    for(j in c((i + 1):N)){
      m <- ((N-1)+(N-i+1))*(i-1)/2 + j-i
      if (work_sim[i,j] == 1){
        suff_mat[m, 1] <- 1 
        if(adv_sim[i,j] == 1){
          suff_mat[m, 4] <- 1 
          # if(friend_sim[i,j] == 1){
          #   suff_mat[m, 7] <- 1 
          # }
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
#' @param data 
#' @param N 
#' @param H 
#' @param mdim 
#'
#' @return
#' @export
#'
#' @examples
simulate_suffstats_Lazega <- function(data = Lazega_lawyer_network, N = 71, H = 2, mdim = 6, m = 10, theta = c(-1.4497737, -3.3336582, -2.6945731,  1.8011888,  0.2176545,  2.4582079)){
  # the default theta is the MPLE of the Lazega network data
  reproduced_suff <- matrix(0,m,mdim)
  seeds <- sample(1:9999999,m,replace = FALSE)
  sim_suff <- matrix(0,m, mdim)
  unique_dyads <- unique(data[,1:2])
  unique_dyads_vec <- as.vector(t(unique_dyads))
  K <- max(data[,3])
  for(iter in c(1:m)){
    seed <- seeds[iter]
    reproduced <- samp_ml(theta, N = N, k = K, mdim = mdim,
                          mterm = 'other',intv = 3, H = H, seed = seeds[iter],
                          basis_arguments = unique_dyads_vec)
    
    simsamp <- reproduced$suff_stats
    reproduced_suff[iter,] <- simsamp
    
  
    
  }
  return(reproduced_suff) 
}


#' Title
#'
#' @param reproduced_suff 
#' @param obs 
#' @param dimx_name 
#'
#' @return
#' @export
#'
#' @examples
draw_box_plot_Lazega <- function(reproduced_suff, obs, dimx_name = c("Coworker", "Advice", "Friendship", "CxA","CxF","AxF")){
  suff_y <- c(reproduced_suff)
  suff_mean <- colMeans(reproduced_suff)
  m = length(reproduced_suff[,1])
  dim_x <- dimx_name #,"CxAxF")
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






