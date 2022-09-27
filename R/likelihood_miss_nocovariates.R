#' Negative log-likelihood across individuals and time
#' @param dat A dataset returned by generate_dataset()
#' @param par Vector of parameters governing the distribution.
#' @return Numeric likelihood
#' @notes This corresponds to the missing data structure
negloglik_miss_nocovariates <- function(dat, par) {
  
  # Convert parameter vector to a named list
  p <- par
  params <- list(a_x=p[1], a_y=p[2], beta=p[3], a_v=p[4])
  
  # Compute the negative likelihood across individuals
  n <- attr(dat, "n")
  -1 * sum(log(unlist(lapply(c(1:n), function(i) {
    
    dat_i <- filter(dat, id==i)
    J <- nrow(dat_i)
    
    # Calculate vectors for patient i
    y <- dat_i$y
    v <- dat_i$v
    u <- dat_i$u
    d <- dat_i$d
    xs <- dat_i$xs
    u_prev <- c(0, dat_i$u[c(1:(J-1))])
    
    # Calculate the set X_i to sum over
    X_i_set <- list()
    for (j in c(1:(J+1))) {
      x_ <- c(rep(0,J-j+1), rep(1,j-1))
      case_i_ <- case(x_,v)
      T_pm <- T_plusminus(case_i_, J, x_, v)
      if (case_i_!=9 &&
          all(xs==x_*d) &&
          all(d==g_delta(case_i_, J, T_pm$T_minus, T_pm$T_plus)) &&
          prod(u[1:J]==u_prev[1:J]+(1-u_prev[1:J])*v[1:J]*x_[1:J])
      ) {
        X_i_set <- c(X_i_set, list(x_))
      }
    }
    
    # Compute the likelihood for individual i
    f2 <- sum(unlist(lapply(X_i_set, function(x) {
      x_prev <- c(0,x[1:(length(x)-1)])
      prod(unlist(lapply(c(1:J), function(j) {
        f_x2(x=x[j], x_prev=x_prev[j], params=params) *
          f_y2(y=y[j], x=x[j], params=params) *
          f_v2(v=v[j], u_prev=u_prev[j], params=params)
      })))
    })))
    if (f2<=0) {
      f2 <- 1e-10
      # warning("Likelihood of zero")
    }
    
    return(f2)
    
  }))))
  
}



#' Calculate likelihood component f_x2
#'
#' @param x Seroconversion indicator (time j)
#' @param x_prev Seroconversion indicator (time j-1)
#' @param params Named list of parameters
#' @return Numeric likelihood
f_x2 <- function(x, x_prev, params) {
  p <- params
  if (x==1) {
    if (x_prev==1) {
      return(1)
    } else {
      return(min(exp(p$a_x),0.99999))
    }
  } else {
    if (x_prev==1) {
      return(0)
    } else {
      return(1 - min(exp(p$a_x),0.99999))
    }
  }
}



#' Calculate likelihood component f_y2
#'
#' @param y Event indicator (time j)
#' @param x Seroconversion indicator (time j)
#' @param params Named list of parameters
#' @return Numeric likelihood
f_y2 <- function(y, x, params) {
  p <- params
  explin <- min(exp(p$a_y + p$beta*x),0.99999)
  if (y==1) { return(explin) } else { return(1-explin) }
}



#' Calculate likelihood component f_v2
#'
#' @param v Testing indicator (time j)
#' @param u_prev Variable U indicator (time j-1)
#' @param params Named list of parameters
#' @return Numeric likelihood
f_v2 <- function(v, u_prev, params) {
  p <- params
  if (v==1) {
    if (u_prev==1) {
      warning("v==1 and u_prev==1")
      return(0.00001)
    } else {
      return(min(exp(p$a_v),0.99999))
    }
  } else {
    if (u_prev==1) {
      return(1)
    } else {
      return(1-min(exp(p$a_v),0.99999))
    }
  }
}
