#' Negative log-likelihood across individuals and time
#' @param dat A dataset returned by generate_dataset()
#' @param par Vector of parameters governing the distribution.
#' @return Numeric likelihood
#' @notes This corresponds to the missing data structure
negloglik_miss <- function(dat, par) {
  
  # Convert parameter vector to a named list
  p <- as.numeric(par)
  params <- list(a_x=p[1], g_x=c(p[2],p[3]), a_y=p[4], g_y=c(p[5],p[6]),
                 beta_x=p[7], beta_z=p[8])
  
  # Compute the negative likelihood across individuals
  n <- attr(dat, "n")
  -1 * sum(log(unlist(lapply(c(1:n), function(i) {
    
    dat_i <- filter(dat, id==i)
    J <- nrow(dat_i)
    
    # Calculate vectors for patient i
    w <- subset(dat_i, select=c(w_sex,w_age)) # !!!!! TEMP
    y <- dat_i$y
    z <- dat_i$z
    v <- dat_i$v
    d <- dat_i$d
    u <- dat_i$u

    # Calculate the set X_i to sum over
    X_i_set <- list()
    for (j in c(1:(J+1))) {
      x_ <- c(rep(0,J-j+1), rep(1,j-1))
      case_i_ <- case(x_,v)
      T_pm <- T_plusminus(case_i_, J, x_, v)
      if (case_i_!=9 &&
          all(u==x_*d) &&
          all(d==g_delta(case_i_, J, T_pm$T_minus, T_pm$T_plus))
      ) {
        X_i_set <- c(X_i_set, list(x_))
      }
    }
    
    # Compute the likelihood for individual i
    w <- t(w)
    f2 <- sum(unlist(lapply(X_i_set, function(x) {
      x_prev <- c(0,x[1:(length(x)-1)])
      prod(unlist(lapply(c(1:J), function(j) {
        w_ij <- as.numeric(w[,j])
        f_x(x=x[j], x_prev=x_prev[j], w=w_ij, params=params) *
          f_y(y=y[j], x=x[j], w=w_ij, z=z[j], params=params)
      })))
    })))
    if (f2<=0) {
      f2 <- 1e-10
      # warning("Likelihood of zero")
    }
    
    return(f2)
    
  }))))
  
}



#' Negative log-likelihood across individuals and time
#' @param dat A dataset returned by generate_dataset()
#' @param par Vector of parameters governing the distribution.
#' @return Numeric likelihood
#' @notes This corresponds to the missing data structure
negloglik_miss_old <- function(dat, par) {
  
  # Convert parameter vector to a named list
  p <- as.numeric(par)
  params <- list(a_x=p[1], g_x=c(p[2],p[3]), a_y=p[4], g_y=c(p[5],p[6]),
                 beta_x=p[7], beta_z=p[8])
  
  # Compute the negative likelihood across individuals
  n <- attr(dat, "n")
  -1 * sum(log(unlist(lapply(c(1:n), function(i) {
    
    dat_i <- filter(dat, id==i)
    J <- nrow(dat_i)
    
    # Calculate vectors for patient i
    w <- subset(dat_i, select=c(w_sex,w_age)) # !!!!! TEMP
    y <- dat_i$y
    z <- dat_i$z
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
    w <- t(w)
    f2 <- sum(unlist(lapply(X_i_set, function(x) {
      x_prev <- c(0,x[1:(length(x)-1)])
      prod(unlist(lapply(c(1:J), function(j) {
        w_ij <- as.numeric(w[,j])
        f_x(x=x[j], x_prev=x_prev[j], w=w_ij, params=params) *
          f_y(y=y[j], x=x[j], w=w_ij, z=z[j], params=params)
      })))
    })))
    if (f2<=0) {
      f2 <- 1e-10
      # warning("Likelihood of zero")
    }
    
    return(f2)
    
  }))))
  
}



#' Calculate likelihood component f_x
#'
#' @param x Seroconversion indicator (time j)
#' @param x_prev Seroconversion indicator (time j-1)
#' @param w Vector of covariates (time j)
#' @param params Named list of parameters
#' @return Numeric likelihood
f_x <- function(x, x_prev, w, params) {
  p <- params
  if (x==1) {
    if (x_prev==1) {
      return(1)
    } else {
      return(exp(p$a_x + sum(p$g_x*w)))
    }
  } else {
    if (x_prev==1) {
      return(0)
    } else {
      return(1 - exp(p$a_x + sum(p$g_x*w)))
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
f_y <- function(y, x, w, z, params) {
  p <- params
  explin <- exp(p$a_y + sum(p$g_y*w) + p$beta_x*x + p$beta_z*z)
  if (y==1) { return(explin) } else { return(1-explin) }
}
