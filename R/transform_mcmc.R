#' Transform dataset into JAGS format
#'
#' @param dat_baseline A dataset returned by generate_data_baseline()
#' @param dat_events A dataset returned by generate_data_events()
#' @return A dataset in JAGS format

transform_mcmc <- function(dat_baseline, dat_events) {
  
  I <- attr(dat_baseline, "num_patients")
  
  # Create empty data structures
  mtx <- cbind(rep(0,I), matrix(NA, nrow=I, ncol=12*(C$end_year-C$start_year)))
  dat_mcmc <- list(
    I = I,
    T_i = c(),
    max_T_i = NULL,
    sex = c(),
    b_age = c(),
    v = mtx,
    x = mtx,
    y = mtx,
    z = mtx
  )
  
  # Fill in dat_mcmc
  for (i in 1:I) {
    dat_mcmc$v[i,] <- dat_events[[i]]$v
    dat_mcmc$x[i,] <- dat_events[[i]]$x
    dat_mcmc$y[i,] <- dat_events[[i]]$y
    dat_mcmc$z[i,] <- dat_events[[i]]$z
    dat_mcmc$T_i <- c(dat_mcmc$T_i, dat_events[[i]]$T_i)
  }
  dat_mcmc$b_age <- dat_baseline$b_age
  dat_mcmc$sex <- dat_baseline$sex
  dat_mcmc$max_T_i <- max(dat_mcmc$T_i)
  
  return(dat_mcmc)
  
}
