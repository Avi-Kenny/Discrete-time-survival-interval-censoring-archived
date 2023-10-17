#' Negative log-likelihood across individuals and time
#' @param dat A dataset returned by generate_dataset()
#' @param par Vector of parameters governing the distribution.
#' @return Numeric likelihood
#' @notes This corresponds to the missing data structure
construct_negloglik_miss <- function(dat) {
  
  # Construct a data object specific to each individual
  n <- attr(dat, "n")
  dat_objs <- lapply(c(1:n), function(i) {
    
    d <- list()
    d$dat_i <- dat[dat$id==i,] # Does dat_i need to be returned with r?
    
    # Start and end times
    d$s_i <- min(d$dat_i$t_end)
    d$t_i <- max(d$dat_i$t_end)
    
    # Data vectors
    d$w <- subset(d$dat_i, select=c(w_sex,w_age)) # !!!!! Subset command is temporary
    d$y <- d$dat_i$y
    d$z <- d$dat_i$z
    d$v <- d$dat_i$v
    d$d <- d$dat_i$d
    d$u <- d$dat_i$u
    d$cal_time <- d$dat_i$t_end / 100 # Rescaling is done to prevent optimization issues
    
    # Calculate the set X_i to sum over
    d$X_i_set <- list()
    for (j in c(1:(d$t_i-d$s_i+2))) {
      x_ <- c(rep(0,d$t_i-d$s_i-j+2), rep(1,j-1))
      case_i_ <- case(x_,d$v)
      T_pm <- T_plusminus(case=case_i_, s_i=d$s_i, t_i=d$t_i, x=x_, v=d$v)
      if (case_i_!=9 &&
          all(d$u==x_*d$d) &&
          all(d$d==g_delta(case=case_i_, s_i=d$s_i, t_i=d$t_i,
                           T_minus=T_pm$T_minus, T_plus=T_pm$T_plus))
      ) {
        d$X_i_set <- c(d$X_i_set, list(x_))
      }
    }
    
    return(d)
    
  })
  
  negloglik_miss <- function(par) {
    
    # Convert parameter vector to a named list
    p <- as.numeric(par)
    params <- list(a_x=p[1], g_x=c(p[2],p[3]), a_y=p[4], g_y=c(p[5],p[6]),
                   beta_x=p[7], beta_z=p[8], t_x=p[9], t_y=p[10],
                   a_s=p[11], t_s=p[12], g_s=c(p[13],p[14]))
    
    # Compute the negative likelihood across individuals
    -1 * sum(log(unlist(lapply(c(1:n), function(i) {
      
      # Extract data for individual i
      d <- dat_objs[[i]]

      # Compute the likelihood for individual i
      w <- t(d$w)
      f2 <- sum(unlist(lapply(d$X_i_set, function(x) {
        x_prev <- c(0,x[1:(length(x)-1)])
        prod(unlist(lapply(c(1:(d$t_i-d$s_i+1)), function(j) {
          # Note: here, j is the index of time within an individual rather than
          #       a calendar time index
          w_ij <- as.numeric(w[,j])
          return(
            f_x(x=x[j], x_prev=x_prev[j], w=w_ij, j=d$cal_time[j],
                s=In(d$cal_time[j]==d$s_i), params=params) *
              f_y(y=d$y[j], x=x[j], w=w_ij, z=d$z[j], j=d$cal_time[j],
                  params=params)
          )
        })))
      })))
      if (f2<=0) {
        f2 <- 1e-10
        # warning("Likelihood of zero")
      }
      
      return(f2)
      
    }))))
    
  }
  
  return(negloglik_miss)
  
}



#' Calculate likelihood component f_x
#'
#' @param x Seroconversion indicator (time j)
#' @param x_prev Seroconversion indicator (time j-1)
#' @param w Vector of covariates (time j)
#' @param params Named list of parameters
#' @return Numeric likelihood
f_x <- function(x, x_prev, w, j, s, params) {
  if (s==0) {
    if (x==1) {
      if (x_prev==1) {
        return(1)
      } else {
        return(exp(params$a_x + params$t_x*j + sum(params$g_x*w)))
      }
    } else {
      if (x_prev==1) {
        return(0)
      } else {
        return(1 - exp(params$a_x + params$t_x*j + sum(params$g_x*w)))
      }
    }
  } else {
    expitlin <- expit(params$a_s + params$t_s*j + sum(params$g_s*w)) # !!!!! Need to add these params
    if (x==1) {
      return(expitlin)
    } else {
      return(1-expitlin)
    }
  }
}



#' Calculate likelihood component f_y
#'
#' @param y Event indicator (time j)
#' @param x Seroconversion indicator (time j)
#' @param w Vector of covariates (time j)
#' @param z ART indicator (time j)
#' @param params Named list of parameters
#' @return Numeric likelihood
f_y <- function(y, x, w, z, j, params) {
  explin <- exp(params$a_y + params$t_y*j + sum(params$g_y*w) + params$beta_x*x +
                  params$beta_z*z)
  if (y==1) { return(explin) } else { return(1-explin) }
}
