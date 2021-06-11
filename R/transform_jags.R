#' Transform dataset into JAGS format
#'
#' @param dat_baseline A dataset returned by generate_data_baseline()
#' @param dat_events A dataset returned by generate_data_events()
#' @return A dataset in JAGS format
#' @notes
#'     - The X and Z matrices receive a "baseline" column

transform_jags <- function(dat_baseline, dat_events) {
  
  I <- attr(dat_baseline, "num_patients")
  
  # Create empty data structures
  mtx <- cbind(rep(0,I), matrix(NA, nrow=I, ncol=12*(C$end_year-C$start_year)))
  dat_jags <- list(
    I = I,
    J = c(),
    max_J = NULL,
    sex = c(),
    b_age = c(),
    v = mtx,
    x = mtx,
    y = mtx,
    z = mtx
  )
  
  # Fill in dat_jags
  for (i in 1:I) {
    dat_jags$v[i,] <- c(0,dat_events[[i]]$v)
    dat_jags$x[i,] <- c(0,dat_events[[i]]$x)
    dat_jags$y[i,] <- c(0,dat_events[[i]]$y)
    dat_jags$z[i,] <- c(0,dat_events[[i]]$z)
    dat_jags$J <- c(dat_jags$J, dat_events[[i]]$J)
  }
  dat_jags$b_age <- dat_baseline$b_age
  dat_jags$sex <- dat_baseline$sex
  dat$max_J <- max(dat$J)
  
  # Recode NA values
  dat$v %<>% replace_na(999)
  dat$x %<>% replace_na(999)
  dat$y %<>% replace_na(999)
  dat$z %<>% replace_na(999)
  
  return(dat_jags)
  
}
