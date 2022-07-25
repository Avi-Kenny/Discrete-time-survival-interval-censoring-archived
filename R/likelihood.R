#' Calculate overall likelihood
#'
#' @param dat A dataset returned by generate_dataset()
#' @param params List of parameters governing the distribution.
#' @return Numeric likelihood
lik <- function(dat, params) {
  
  p <- params
  n <- attr(dat, "n")
  sum(unlist(lapply(c(1:n), function(i) {
    
    dat_i <- filter(dat, id==i)
    J <- nrow(dat_i)
    
    # Calculate the set X_i to sum over
    X_set_i <- list(999)
    case_i <- case(dat$x,dat$v)
    
    # Calculate vectors for patient i
    x <- dat_i$x
    x_prev <- c(0, x[c(1:(J-1))])
    w <- subset(dat, select=c(w_sex,w_age)) # !!!!! TEMP
    names(w) <- c("w1","w2") # !!!!! TEMP
    y <- dat_i$y
    v <- dat_i$v
    u_prev <- c(0, dat_i$u[c(1:(J-1))])
    
    # Compute the likelihood for individual i
    f2 <- sum(unlist(lapply(X_set_i, function(x) {
      prod(unlist(lapply(c(1:J), function(j) {
        w_ij <- as.numeric(w[j,])
        f_x(x=x[j], x_prev=x_prev[j], w=w_ij) *
          f_y(y=y[j], x=x[j], w=w_ij) *
          f_v(v=v[j], u_prev=u_prev[j], w=w_ij)
      })))
    })))
    return(log(f2))
    
  })))
  
}



#' Calculate likelihood component f_x
#'
#' @param x Seroconversion indicator (time j)
#' @param x_prev Seroconversion indicator (time j-1)
#' @param w Vector of covariates (time j)
#' @return Numeric likelihood
f_x <- function(x, x_prev, w) {
  if (x==1) {
    if (x_prev==1) {
      return(1)
    } else {
      return(exp(a_x + sum(p$g_x*w)))
    }
  } else {
    if (x_prev==1) {
      return(0)
    } else {
      return(1 - exp(a_x + sum(p$g_x*w)))
    }
  }
}



#' Calculate likelihood component f_x
#'
#' @param y Event indicator (time j)
#' @param x Seroconversion indicator (time j)
#' @param w Vector of covariates (time j)
#' @return Numeric likelihood
f_y <- function(y, x, w) {
  explin <- exp(a_y + beta*x + sum(p$g_y*w))
  if (y==1) { return(explin) } else { return(1-explin) }
}



#' Calculate likelihood component f_v
#'
#' @param v Testing indicator (time j)
#' @param u_prev Variable U indicator (time j-1)
#' @param w Vector of covariates (time j)
#' @return Numeric likelihood
f_v <- function(v, u_prev, w) {
  if (v==1) {
    if (u_prev==1) {
      return(0)
    } else {
      return(exp(a_v + sum(p$g_v*w)))
    }
  } else {
    if (u_prev==1) {
      return(1)
    } else {
      return(1-exp(a_v + sum(p$g_v*w)))
    }
  }
}


