#' Expit function
#' 
#' @param x Numeric input
#' @return Numeric output
expit <- function(x) {1/(1+exp(-x))}



#' Take a posterior sample of parameter vector theta
#' 
#' @param fit Stan model fit returned by fit_stan()
#' @param size Number of samples to take (m); one per imputation
#' @return A list of parameter samples
posterior_param_sample <- function(fit, size) {
  
  n.iter <- nrow(fit[[1]])
  chains <- sample(c(1,2), size=size, replace=TRUE)
  rows <- sample(c(1:n.iter), size=size, replace=FALSE)
  alpha0 <- psi1 <- psi2 <- c()
  for (i in 1:C$m) {
    alpha0 <- c(alpha0, fit[[chains[i]]][[rows[i],"alpha0"]])
    psi1 <- c(psi1, fit[[chains[i]]][[rows[i],"psi1"]])
    psi2 <- c(psi2, fit[[chains[i]]][[rows[i],"psi2"]])
  }
  
  return(list(
    alpha0=alpha0, psi1=psi1, psi2=psi2
  ))
  
}
