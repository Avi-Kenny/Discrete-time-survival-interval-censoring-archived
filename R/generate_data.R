#' Generate dataset
#'
#' @param n Number of patients in cohort
#' @param max_time Number of time points in the study
#' @param params List of parameters governing the distribution. Should be of the
#'     following form: list(a=c(1,2), b=c(1,2,3))
#' @return A dataframe, one row per patient, containing the following fields:
#'     - id: patient ID variable
#' @notes
#'     - TO DO

generate_data <- function(n, max_time, params) {
  
  # !!!!!
  if (F) {
    params <- list(
      g_x = c(log(1.3),log(1.002)),
      g_y = c(log(1.2),log(1.001)),
      beta = log(1.5)
    )
    dat <- generate_data(n=100, max_time=100, params=params)
    sum(dplyr::summarize(group_by(dat,id), x=max(x))$x) # Number seroconverted
    sum(dat$y) # Number of events
    sum(dat$x*dat$y) # Number of events among HIV+
    
  }
  
  # Set baseline hazard functions
  b_hazard_x <- function(t) { 0.005 } # !!!!! alpha_j
  b_hazard_y <- function(t) { 0.003 } # !!!!! beta_j
  # b_hazard_x <- function(t) { 0.003 + 0.003*(t/100) }
  # b_hazard_y <- function(t) { 0.002 + 0.002*(t/100) }
  
  # Generate dataframe to hold results
  dat <- data.frame(
    "id" = integer(),
    "t_start" = integer(),
    "t_end" = integer(),
    "w_sex" = integer(),
    "w_age" = integer(),
    "x" = integer(),
    "y" = integer()
  )
  
  # Generate baseline covariates
  id <- c(1:n)
  w_sex <- sample(c(0,1), size=n, replace=T)
  w_age <- sample(c(1:80), size=n, replace=T)
  
  # Loop through individuals/time to generate events
  p <- params
  for (i in c(1:n)) {
    
    # Create variables for serostatus (x) and outcome indicator (y)
    x <- y <- x_prev <- c()
    event <- 0
    t <- 1
    w_sex_ <- w_sex[i]
    w_age_ <- w_age[i]
    
    while (!event && t<=max_time) {
      
      # Sample seroconversion (x)
      if (t==1) {
        p_sero <- b_hazard_x(t) * exp(p$g_x[1]*w_sex_ + p$g_x[2]*w_age_)
        x[t] <- rbinom(n=1, size=1, prob=p_sero)
        x_prev <- 0
      } else {
        x_prev[t] <- x[round(t-1)]
        if (x_prev[t]==1) { x[t] <- 1 } else {
          p_sero <- b_hazard_x(t) * exp(p$g_x[1]*w_sex_ + p$g_x[2]*w_age_)
          x[t] <- rbinom(n=1, size=1, prob=p_sero)
        }
      }
      
      # Sample events
      p_event <- b_hazard_y(t) *
        exp(p$g_y[1]*w_sex_ + p$g_y[2]*w_age_ + p$b*x[t])
      event <- rbinom(n=1, size=1, prob=p_event)
      y[t] <- event
      
      # Increment time
      t <- round(t+1)
      
    }
    
    # Add results to dataframe
    len <- length(x)
    x_prev <- 
    dat <- rbind(dat, list(
      id=rep(i,len), t_start=c(0:(len-1)), t_end=c(1:len),
      w_sex=rep(w_sex_,len), w_age=rep(w_age_,len), x=x, y=y, x_prev=x_prev
    ))
    
  }
  
  attr(dat, "n") <- n
  attr(dat, "max_time") <- max_time
  attr(dat, "params") <- params
  
  return (dat)
  
}
