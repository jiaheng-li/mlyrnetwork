
draw_mlnet <- function(mlnet,N){
  
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  k <- max(unique(mlnet[,3]))
  node_colors <- sample(col_vector, k)
  fig <- vector("list",length = k)
  net <- list()
  for(j in 1:k){
    layer <- mlnet[mlnet[,3] == j,]
    n <- length(layer[,1])
    net <- append(net,list(matrix(0,N,N)))
    for(i in 1:n){
      ind1 <- layer[i,1]
      ind2 <- layer[i,2]
      net[[j]][ind1,ind2] <- 1
      net[[j]][ind2, ind1] <-  net[[j]][ind1,ind2]
      net[[j]] <- as.network(net[[j]], directed = FALSE)
    }
    
    if(j == 1){
      x = sna::gplot.layout.fruchtermanreingold(net[[j]], NULL)
    }
         
    net[[j]] %v% "x" = x[, 1]
    net[[j]] %v% "y" = x[, 2]
    fig[[j]] <- GGally::ggnet2(net[[j]], node.size = 3, node.color = node_colors[j], mode = c("x", "y") ) + 
      ggplot2::labs(title = paste0("Layer ", as.character(j))) +  ## Change the graph title here
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size = 10))   ## Change the font size here  
        
  }
  
  ## Arrange all networks in one plot
  ggpubr::ggarrange(plotlist = fig, ncol = k, nrow = 1)
  

}



