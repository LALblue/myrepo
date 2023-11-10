sqrtmm <-function(A){
  a.eig <- eigen(A)
  return(a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors))
}

#int x to infty
incomplete_gammafunction <- function(x, p){
  lgamma(p/2) + pgamma(x, p/2, 1, lower = FALSE, log = TRUE)
}

LOGDET <- function(x){
  return(log(det(x)))
}

# Compute the adj(A) (transpose of cofactor matrix, adjoint of matrix)
compute_adjugate <- function(A) {
  # Ensure A is a square matrix
  if (nrow(A) != ncol(A)) {
    stop("Matrix must be square")
  }
  
  n <- nrow(A)
  
  # Initialize cofactor matrix
  C <- matrix(0, n, n)
  
  # Compute each cofactor
  for (i in 1:n) {
    for (j in 1:n) {
      minor_matrix <- A[-i, -j]
      C[i, j] <- (-1)^(i + j) * det(minor_matrix)
    }
  }
  
  # Return the adjugate (transpose of the cofactor matrix)
  return(t(C))
}
# 
# # Test
# A <- matrix(c(1, 2, 3, 0, 4, 5, 1, 0, 6), nrow=3)
# print(compute_adjugate(A))


sqrtmm <-function(A){
  a.eig <- eigen(A)
  return(a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors))
}
# LoRaD general version
# Only need to input an unconstrain matrix T by p MCMC parameter samples and a p by 1 vector of loglikelihood
LoRaD <- function(MCMCmatrix, LOGLIKELIHOOD, rq = 0.2){
  eta1 <- MCMCmatrix[1:(nrow(MCMCmatrix)/2),]
  eta2 <- MCMCmatrix[-(1:(nrow(MCMCmatrix)/2)),]
  loglikelihood <- LOGLIKELIHOOD[-(1:(nrow(MCMCmatrix)/2))]
  p <- ncol(MCMCmatrix)
  MCMCsize <- nrow(MCMCmatrix)
  MCMCsize2 <- MCMCsize/2
  #First half: Cal r_max
  covp <- cov(eta1)
  eta1_bar_matrix <- matrix(nrow = MCMCsize2, ncol = p)
  eta1_bar <- colMeans(eta1)
  for(i in 1:MCMCsize2){
    eta1_bar_matrix[i,]<- eta1_bar
  }
  psi <- solve(sqrtmm(covp)) %*% t(eta1 - eta1_bar_matrix)
  theta <- as.matrix(psi)
  theta <- t(theta) 
  r_sq <- rowSums(theta^2)
  r_sq_max <- quantile(r_sq, rq)
  
  #Second half: Cal MCMC kernel
  covp <- cov(eta2)
  eta2_bar_matrix <- matrix(nrow = MCMCsize2, ncol = p)
  eta2_bar <- colMeans(eta2)
  for(i in 1:MCMCsize2){
    eta2_bar_matrix[i,] <- eta2_bar
  }
  psi <- solve(sqrtmm(covp)) %*% t(eta2 - eta2_bar_matrix)
  log_J_3 <- sum(log(diag(chol(covp))))
  q_star <-rep(0,MCMCsize - MCMCsize2)
  for(i in 1:(MCMCsize - MCMCsize2)){
    q_star[i] <- -0.5 * p * log(2*pi) - 0.5 * sum(psi[,i]^2)
  }
  ########### log likelihood
  theta <-as.matrix(psi)
  theta <-t(theta) 
  r_sq <- rowSums(theta^2)
  q <- loglikelihood + log_J_3
  ratio_star <- q_star - q 
  #r_max_store <- c(r_max_store, sqrt(r_sq_max) )
  keep <- which(r_sq <= r_sq_max)
  #no_points <- c(no_points, length(keep) )
  theta <- theta[keep, ]
  r_sq <- r_sq[keep]
  ratio_star <- ratio_star[keep]
  #q: posterior kernel of the standardized MCMC sample
  #q_star: proposed posterior kernel based on the standard multivariate normal
  all <- cbind(theta, sqrt(r_sq), ratio_star, r_sq)
  temp <- all[, (p + 2)]
  d_bot_star <- log(sum(exp(temp - max(temp) ) ) ) + max(temp) - log(MCMCsize2) #q^*
  d_top_star <- log(pgamma(r_sq_max/2,p/2) ) 
  result_aPWK <- d_top_star - d_bot_star
  return(result_aPWK)
}






