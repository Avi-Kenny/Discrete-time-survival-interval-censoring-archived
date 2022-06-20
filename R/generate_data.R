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
    params <- list(a=c(log(1.3),log(1.002)), b=c(log(1.2),log(1.001),log(1.5)))
    dat2 <- generate_data(n=100, max_time=100, params=params)
    sum(summarize(group_by(dat2,id), x=max(x))$x) # Number seroconverted
    sum(dat2$d) # Number of events
  }
  
  # Set baseline hazard functions
  b_hazard_x <- function(t) { 0.003 + 0.003*(t/100) }
  b_hazard_y <- function(t) { 0.002 + 0.002*(t/100) }
  
  # Generate dataframe to hold results
  dat <- data.frame(
    "id" = integer(),
    "t_start" = integer(),
    "t_end" = integer(),
    "z_sex" = integer(),
    "z_age" = integer(),
    "x" = integer(),
    "d" = integer()
  )
  
  # Generate baseline covariates
  id <- c(1:n)
  z_sex <- sample(c(0,1), size=n, replace=T)
  z_age <- sample(c(1:80), size=n, replace=T)
  
  # Loop through individuals/time to generate events
  p <- params
  for (i in c(1:n)) {
    
    # Create variables for serostatus (x) and event indicator (d)
    x <- d <- c()
    event <- 0
    t <- 1
    z_sex_ <- z_sex[i]
    z_age_ <- z_age[i]
    
    while (!event && t<=max_time) {
      
      # Sample seroconversion (x)
      if (t==1) {
        p_sero <- b_hazard_x(t) * exp(p$a[1]*z_sex_) * exp(p$a[2]*z_age_)
        x[t] <- rbinom(n=1, size=1, prob=p_sero)
      } else {
        if (x[round(t-1)]==1) { x[t] <- 1 } else {
          p_sero <- b_hazard_x(t) * exp(p$a[1]*z_sex_) * exp(p$a[2]*z_age_)
          x[t] <- rbinom(n=1, size=1, prob=p_sero)
        }
      }
      
      # Sample events
      p_event <- b_hazard_y(t) * exp(p$b[1]*z_sex_) * exp(p$b[2]*z_age_) *
        exp(p$b[3]*x[t])
      event <- rbinom(n=1, size=1, prob=p_event)
      d[t] <- event
      
      # Increment time
      t <- round(t+1)
      
    }
    
    # Add results to dataframe
    len <- length(x)
    dat <- rbind(dat, list(
      id=rep(i,len), t_start=c(0:(len-1)), t_end=c(1:len),
      z_sex=rep(z_sex_,len), z_age=rep(z_age_,len), x=x, d=d
    ))
    
  }
  
  attr(dat, "n") <- n
  attr(dat, "max_time") <- max_time
  attr(dat, "params") <- params
  
  return (dat)
  
}
