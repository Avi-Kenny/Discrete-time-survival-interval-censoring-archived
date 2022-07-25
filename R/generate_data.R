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
  a_x <- function(t) { log(0.005) }
  a_y <- function(t) { log(0.003) }
  a_v <- function(t) { log(0.01) }
  # b_hazard_x <- function(t) { log(0.003 + 0.003*(t/100)) }
  # b_hazard_y <- function(t) { log(0.002 + 0.002*(t/100)) }
  
  # Generate dataframe to hold results
  dat <- data.frame(
    "id" = integer(),
    "t_start" = integer(),
    "t_end" = integer(),
    "w_sex" = integer(),
    "w_age" = integer(),
    "x" = integer(),
    "y" = integer(),
    "v" = integer(),
    "u" = integer(),
    "d" = integer(),
    "xs" = integer()
  )
  
  # Generate baseline covariates
  id <- c(1:n)
  w_sex <- sample(c(0,1), size=n, replace=T)
  w_age <- sample(c(1:80), size=n, replace=T)
  
  # Loop through individuals/time to generate events
  p <- params
  i_T_minus <- i_T_plus <- i_case <- c()
  for (i in c(1:n)) {
    
    x <- y <- v <- u <- c() # !!!!! check if all of these are needed
    event <- u_prev <- x_prev <- 0
    j <- 0
    w_sex_ <- w_sex[i]
    w_age_ <- w_age[i]
    
    while (!event && j<=max_time) {
      
      # Increment time
      j <- round(j+1)
      
      # Sample seroconversion (x)
      p_sero <- x_prev + (1-x_prev) * exp(
        a_x(j) + p$g_x[1]*w_sex_ + p$g_x[2]*w_age_
      )
      x[j] <- x_prev <- rbinom(n=1, size=1, prob=p_sero)
      
      # Sample events
      p_event <- exp(a_y(j) + p$g_y[1]*w_sex_ + p$g_y[2]*w_age_ + p$b*x[j])
      event <- rbinom(n=1, size=1, prob=p_event)
      y[j] <- event
      
      # Sample testing
      p_test <- (1-u_prev) * exp(a_v(j) + p$g_v[1]*w_sex_ + p$g_v[2]*w_age_)
      v[j] <- rbinom(n=1, size=1, prob=p_test)
      
      # Calculate additional variables
      u[j] <- u_prev + (1-u_prev)*v[j]*x[j]
      u_prev <- u[j]
      
    }
    J <- j
    
    # Calculate case indicators
    case_i <- case(x,v)
    
    # Calculate time(s) of most recent negative test and/or positive test
    T_pm <- T_plusminus(case_i, J, x, v)
    
    # Calculate Delta
    d <- g_delta(case_i, J, T_pm$T_minus, T_pm$T_plus)
    
    # Add results to dataframe
    dat <- rbind(dat, list(
      id=rep(i,J), t_start=c(0:(J-1)), t_end=c(1:J), w_sex=rep(w_sex_,J),
      w_age=rep(w_age_,J), x=x, y=y, v=v, u=u, d=d, xs=x*d
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
