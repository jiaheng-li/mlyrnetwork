## calculate network satistics

count_edges <- function(){
  
}

count_cooccurrences <- function(){
  
}

count_triangle <- function(){
  
}

# function to compute the inverse square root of a matrix
fnMatSqrtInverse = function(mA) {
  ei = eigen(mA)
  d = ei$values
  d = (d+abs(d))/2
  d2 = 1/sqrt(d)
  d2[d == 0] = 0
  return(ei$vectors %*% diag(d2) %*% t(ei$vectors))
}
