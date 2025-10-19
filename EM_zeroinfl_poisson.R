EM_zeroinfl_poisson <- function(y, k = 3, lambda_hat = NULL, p_hat = NULL, theta_hat = NULL,
                                tol = 1e-6, nitmax = 1000, criterio=1, nit=1) {

  if (is.null(lambda_hat)) {
    cuts <- stats::quantile(y, probs = seq(0, 1, length.out = k + 1))
    grp <- cut(y, breaks = unique(cuts), include.lowest = TRUE)
    lambda_hat <- tapply(y, grp, mean, na.rm = TRUE)
  }
  if (is.null(theta_hat)) theta_hat <- 0.5
  if (is.null(p_hat)) p_hat <- rep(1/k, k)
  
  logLik_mix <- function(y, lambda, theta, p) {
    mat <- sapply(lambda, function(l) dpois(y, l))
    colnames(mat) <- paste0("pois_", 1:length(lambda))
    mix <- rowSums(t(t(mat) * p))
    matty <- theta*(y == 0) + (1 - theta)*mix
    matty[matty <= 0] <- .Machine$double.eps
    sum(log(matty))
  }
  
  poismix <- function(y, lambda, p) {
    mat <- sapply(lambda, function(l) dpois(y, l))
    colnames(mat) <- paste0("pois_", 1:length(lambda))
    t(t(mat) * p)
  }

  loglik.old <- logLik_mix(y, lambda_hat, theta_hat, p_hat)

  
  repeat {
    # Expectation
    val <- rowSums(poismix(y, lambda_hat, p_hat))
    probEsti <- numeric(length(y))
    probEsti[y == 0] <- theta_hat / (theta_hat + (1 - theta_hat)*val[y == 0])
    val <- poismix(y, lambda_hat, p_hat)
    probs <- (1 - probEsti) * (val / rowSums(val))
    
    # Maximization
    p_hat <- colMeans(probs)
    theta_hat <- mean(probEsti)
    lambda_hat <- colSums(probs * y) / colSums(probs)
    
    # Stop
    loglik <- logLik_mix(y, lambda_hat, theta_hat, p_hat)
    criterio <- abs(loglik - loglik.old)
    loglik.old <- loglik
    nit <- nit + 1
    
    if (criterio < tol || nit >= nitmax) break
  }
  
  list(lambda = lambda_hat,
       mixing = p_hat,
       theta = theta_hat,
       loglik = loglik.old,
       nit = nit,
       difference = criterio)
}
