#' Generate dataset
#'
#' @param n Number of patients in cohort
#' @param max_time Number of time points in the study
#' @param params List of parameters governing the distribution.
#' @return A dataframe, one row per patient, containing the following fields:
#'     - id: patient ID variable
#' @notes
#'     - TO DO

generate_data <- function(n, max_time, params) {
  
  # Set baseline hazard functions
  # Note: these have temporarily been set to constants, passed in through the
  #     `params` argument
  # a_x <- function(t) { log(0.005) }
  # a_y <- function(t) { log(0.003) }
  # a_v <- function(t) { log(0.1) }
  
  # Generate dataframe to hold results
  dat <- data.frame(
    "id" = integer(),
    "t_start" = integer(),
    "t_end" = integer(),
    "w_sex" = integer(),
    "w_age" = double(),
    "x" = integer(),
    "z" = integer(),
    "y" = integer(),
    "v" = integer(),
    "d" = integer(),
    "u" = integer()
  )
  
  # Generate baseline covariates
  id <- c(1:n)
  w_sex <- sample(c(0,1), size=n, replace=T)
  w_age <- sample(c(1:80), size=n, replace=T)
  
  # !!!!! Temp: scale age variable
  w_age <- w_age/100
  
  # Loop through individuals/time to generate events
  p <- params
  i_T_minus <- i_T_plus <- i_case <- c()
  for (i in c(1:n)) {
    
    x <- y <- v <- z <- c()
    event <- x_prev <- z_prev <- 0
    known_pos <- known_pos_prev <- 0
    j <- 0
    w_sex_ <- w_sex[i]
    w_age_ <- w_age[i]
    
    while (!event && j<=max_time) {
      
      # Increment time
      j <- round(j+1)
      
      # Sample serostatus (x)
      p_sero <- x_prev + (1-x_prev) * exp(
        p$a_x + p$g_x[1]*w_sex_ + p$g_x[2]*w_age_
      )
      x[j] <- x_prev <- rbinom(n=1, size=1, prob=p_sero)
      
      # Sample testing
      p_test <- (1-known_pos_prev) * exp(
        p$a_v + p$g_v[1]*w_sex_ + p$g_v[2]*w_age_
      )
      v[j] <- rbinom(n=1, size=1, prob=p_test)
      
      # Calculate "known positive" variable
      known_pos[j] <- known_pos_prev + (1-known_pos_prev)*v[j]*x[j]
      known_pos_prev <- known_pos[j]
      
      # Sample ART status
      if (known_pos[j]==0) {
        p_art <- 0
      } else if (z_prev==1) {
        p_art <- 1
      } else {
        p_art <- exp(
          p$a_z + p$g_z[1]*w_sex_ + p$g_z[2]*w_age_
        )
      }
      z[j] <- z_prev <- rbinom(n=1, size=1, prob=p_art)
      
      # Sample events
      p_event <- exp(
        p$a_y + p$g_y[1]*w_sex_ + p$g_y[2]*w_age_ +
          p$beta_x*x[j] + p$beta_z*z[j]
      )
      event <- rbinom(n=1, size=1, prob=p_event)
      y[j] <- event
      
    }

    # Calculate case indicators
    case_i <- case(x,v)
    
    # Calculate time(s) of most recent negative test and/or positive test
    T_pm <- T_plusminus(case_i, j, x, v)
    
    # Calculate Delta
    d <- g_delta(case_i, j, T_pm$T_minus, T_pm$T_plus)
    
    # Add results to dataframe
    dat <- rbind(dat, list(
      id=rep(i,j), t_start=c(0:(j-1)), t_end=c(1:j), w_sex=rep(w_sex_,j),
      w_age=rep(w_age_,j), x=x, z=z, y=y, v=v, d=d, u=x*d
    ))
    
    # Store additional vectors
    i_T_minus[i] <- T_pm$T_minus
    i_T_plus[i] <- T_pm$T_plus
    i_case[i] <- case_i
    
  }
  
  attr(dat, "n") <- n
  attr(dat, "max_time") <- max_time
  attr(dat, "params") <- params
  attr(dat, "T_minus") <- i_T_minus
  attr(dat, "T_plus") <- i_T_plus
  attr(dat, "case") <- i_case
  
  return (dat)
  
}
