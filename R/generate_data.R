#' Generate dataset
#'
#' @param n Number of patients in cohort
#' @param max_time Number of time points in the study
#' @param params List of parameters governing the distribution; see levels.R
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
    "w_1" = integer(),
    "w_2" = double(),
    "x" = integer(),
    "z" = integer(),
    "y" = integer(),
    "v" = integer(),
    "d" = integer(),
    "u" = integer()
  )
  
  # Generate baseline covariates
  id <- c(1:n)
  w_1 <- sample(c(0,1), size=n, replace=T)
  w_2 <- sample(c(13:80), size=n, replace=T)
  
  # Scale age variable
  w_2 <- w_2/100
  
  # Sample start time
  s_i <- sample(c(1:100), size=n, replace=T)
  
  # Loop through individuals/time to generate events
  p <- params
  vec_T_minus <- vec_T_plus <- vec_case <- vec_s_i <- vec_t_i <- c()
  for (i in c(1:n)) {
    
    # Initial values
    x <- y <- v <- z <- c()
    event <- z_prev <- 0
    known_pos <- known_pos_prev <- 0
    j <- 1
    w_1_ <- w_1[i]
    w_2_ <- w_2[i]
    cal_time <- s_i_ <- s_i[i] # cal_time currently unused
    
    while (!event && j<=max_time) {
      
      # !!!!! Need to increment age each loop
      
      # Sample serostatus (x)
      # !!!!! Add calendar time trend
      if (j==1) {
        # Sample baseline serostatus
        # !!!!! Add calendar time trend: + p$t_s*cal_time
        p_sero <- expit(p$a_s + sum(p$g_s*c(w_1_,w_2_)))
      } else {
        p_sero <- x_prev + (1-x_prev) * exp(
          p$a_x + p$g_x[1]*w_1_ + p$g_x[2]*w_2_
        )
      }
      x[j] <- x_prev <- rbinom(n=1, size=1, prob=p_sero)
      
      # Sample testing
      p_test <- (1-known_pos_prev) * exp(
        p$a_v + p$g_v[1]*w_1_ + p$g_v[2]*w_2_
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
          p$a_z + p$g_z[1]*w_1_ + p$g_z[2]*w_2_
        )
      }
      z[j] <- z_prev <- rbinom(n=1, size=1, prob=p_art)
      
      # Sample events
      # !!!!! Add calendar time trend
      p_event <- exp(
        p$a_y + p$g_y[1]*w_1_ + p$g_y[2]*w_2_ +
          p$beta_x*x[j] + p$beta_z*z[j]
      )
      event <- rbinom(n=1, size=1, prob=p_event)
      y[j] <- event
      
      # Increment time
      j <- round(j+1)
      cal_time <- cal_time+1 # Currently unused; also account for scaling
      
    }
    
    # Need to subtract one from j since it is incremented at the end of the loop
    j <- j-1

    # Calculate case indicators
    case_i <- case(x,v)
    
    # Calculate time(s) of most recent negative test and/or positive test
    T_pm <- T_plusminus(case=case_i, s_i=s_i_, t_i=j+s_i_-1, x=x, v=v)

    # Calculate Delta
    d <- g_delta(case=case_i, s_i=s_i_, t_i=s_i_+j-1, T_minus=T_pm$T_minus,
                 T_plus=T_pm$T_plus)
    
    # Add results to dataframe
    dat <- rbind(dat, list(
      id=rep(i,j), t_start=c((s_i_-1):(s_i_+j-2)), t_end=c(s_i_:(s_i_+j-1)),
      w_1=rep(w_1_,j), w_2=rep(w_2_,j), x=x, z=z, y=y, v=v, d=d, u=x*d
    ))
    
    # Store additional vectors
    vec_T_minus[i] <- T_pm$T_minus
    vec_T_plus[i] <- T_pm$T_plus
    vec_case[i] <- case_i
    vec_s_i[i] <- s_i_
    vec_t_i[i] <- j+s_i_-1
    
  }
  
  attr(dat, "n") <- n
  attr(dat, "max_time") <- max_time
  attr(dat, "params") <- params
  attr(dat, "T_minus") <- vec_T_minus
  attr(dat, "T_plus") <- vec_T_plus
  attr(dat, "case") <- vec_case
  attr(dat, "s_i") <- vec_s_i
  attr(dat, "t_i") <- vec_t_i
  
  return (dat)
  
}
