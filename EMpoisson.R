
EM_poismix <- function(y, k = 3, lambda_hat = NULL, p_hat = NULL, tol = 10^-4, nit = 1, nitmax = 1000, criterio = 1){
  # Startings
  if(is.null(p_hat))  p_hat <- rep(1/k, k)
  if (is.null(lambda_hat)) {
    cuts <- stats::quantile(y, probs = seq(0, 1, length.out = k + 1))
    grp <- cut(y, breaks = unique(cuts), include.lowest = TRUE)
    lambda_hat <- tapply(y, grp, mean)
    }

  logLik_mix <- function(y, lambda, p){
    mat <- sapply(lambda, function(lambda) dpois(y, lambda))
    colnames(mat) <- paste0("pois_", 1:length(lambda))
    matty <- t(t(mat) * p)
    loglike <- sum(log(matty))
  }
  poismix <- function(y, lambda, p){
    mat <- sapply(lambda, function(lambda) dpois(y, lambda))
    colnames(mat) <- paste0("pois_", 1:length(lambda))
    matty <- t(t(mat) * p)
  }
  loglik.old <- logLik_mix(y, lambda_hat, p_hat)

  # 
  lambda_hat.new <- lambda_hat
  p_hat.new <- p_hat
  
  while(criterio > tol){
  # Expectation

  numeratore <- poismix(y, lambda_hat, p_hat)
  denominatore <- rowSums(numeratore)
  
  probs <- numeratore / denominatore
  p_hat.new <- colMeans(probs)
  
  # Maximization
  weights <- probs * y
  lambda_hat.new <- colSums(weights)/colSums(probs)
  
  loglik <- abs(logLik_mix(y, lambda_hat.new, p_hat.new))
  
  # Stop 
  criterio <- abs(loglik - loglik.old)
  nit <- nit + 1
  loglik.old <- loglik
  
  estimates <- list(lambda = lambda_hat.new, mixing = p_hat.new, loglik = loglik, nit = nit, difference = criterio)
  if(nitmax == nit) break
  }
  return(estimates)
 }



