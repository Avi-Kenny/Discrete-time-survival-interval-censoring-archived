#' Alias for as.integer
#' 
In <- as.integer

#' Expit function
#' 
#' @param x Numeric input
#' @return Numeric output
expit <- function(x) {1/(1+exp(-x))}



#' Modified exp function (see scratch for derivation)
#' 
#' @param x Numeric input
#' @return Numeric output
exp2 <- (function() {
  
  expit <- function(x) {1/(1+exp(-x))}
  logit <- function(x) { log(x/(1-x)) }
  e <- -0.1
  ell <- logit(exp(e))
  x_0 <- e - (ell*exp(ell))/(exp(e)*(1+exp(ell))^2)
  k_0 <- exp(e-ell)*(1+exp(ell))^2
  exp2 <- function(x) {
    In(x<=e) * exp(x) +
      In(x>e) * expit(k_0*(x-x_0))
  }
  return(exp2)
})()
# exp2 <- exp



#' Compute "case" indicator variables
#' 
#' @param x Vector of seroconversion indicators for one individual
#' @param v Vector of testing indicators for one individual
#' @return An integer (1-4) representing the case, as follows:
#'     - Case 1: No testing data
#'     - Case 2: Only NEG tests
#'     - Case 3: NEG test then POS test
#'     - Case 4: first test POS
case <- function(x,v, warn=F) {
  
  case_i <- rep(NA, 4)
  sum_v <- sum(v)
  sum_vx <- sum(v*x)
  sum_v1x <- sum(v*(1-x))
  case_i[1] <- In(sum_v==0)
  case_i[2] <- In(sum_v>0 & sum_vx==0)
  case_i[3] <- In(sum_v>0 & sum_v1x>0 & sum_vx==1)
  case_i[4] <- In(sum_v==1 & sum_vx==1)
  if (sum(case_i)==0) {
    if (warn) { warning("No applicable case") }
    return(9)
  }
  if (sum(case_i)>1) {
    if (warn) { warning("Multiple cases matched") }
    return(9)
  }
  return(which(case_i==1))
  
}



#' Calculate time(s) of most recent negative test and/or positive test
#' 
#' @param case Case 
#' @param J Maximum time
#' @param x Vector of seroconversion indicators for one individual
#' @param v Vector of testing indicators for one individual
#' @return A list containing T_minus and T_plus
T_plusminus <- function(case, s_i, t_i, x, v) {
  
  len <- t_i-s_i+1
  if (!length(v)==len || !length(x)==len) { stop("Incorrect vector lengths.") }

  if (case %in% c(2,3)) {
    # Time of most recent negative test
    T_minus_index <- which.max(c(s_i:t_i)*v*(1-x))
    T_minus <- c(s_i:t_i)[T_minus_index]
  } else {
    T_minus <- 0
  }
  if (case %in% c(3,4)) {
    # Time of positive test
    T_plus_index <- which.max(v*x)
    T_plus <- c(s_i:t_i)[T_plus_index]
  } else {
    T_plus <- 0
  }
  
  return(list(T_minus=T_minus, T_plus=T_plus))
  
}



#' Compute Delta variable
#' 
#' @param case Case 
#' @param J Maximum time
#' @param x Vector of seroconversion indicators for one individual
#' @param v Vector of testing indicators for one individual
#' @return The vector Delta
g_delta <- function(case, s_i, t_i, T_minus, T_plus) {
  
  if (case==1) {
    d <- rep(0, t_i-s_i+1)
  } else if (case==2) {
    d <- c(rep(1,T_minus-s_i+1), rep(0,t_i-T_minus))
  } else if (case==3) {
    d <- c(rep(1,T_minus-s_i+1), rep(0,T_plus-T_minus-1), rep(1,t_i-T_plus+1))
  } else if (case==4) {
    d <- c(rep(0,T_plus-s_i), rep(1,t_i-T_plus+1))
  } else {
    stop("Error in computing g_delta.")
  }
  if (length(d)!=(t_i-s_i+1)) { stop("Delta is of the wrong length.") }
  
  return(d)
  
}



#' Helper function for debugging; prints timestamps
#' 
#' @param num Number
#' @param msg Message
chk <- function(num, msg="") {
  if (msg=="") {
    str <- paste0("Check ", num, ": ", Sys.time())
  } else {
    str <- paste0("Check ", num, " (", msg, "): ", Sys.time())
  }
  print(str)
}



#' Create a natural cubic spline basis
#' 
#' @param x Numeric input
#' @param i Index of spline column
#' @param m Model version number; corresponds to cfg$model_version
#' @param s Spline number; for cases in which models use multiple splines
construct_basis <- function(m, s) {
  
  if ((m %in% c(9,10)) && s==1) {
    
    grid <- seq(0,1,0.001)
    b <- Vectorize(function(x, i) {
      splines::ns(x=x, knots=c(0.25,0.5,0.75), Boundary.knots=c(0,1))[i]
    }, vectorize.args="x")
    y <- matrix(NA, nrow=length(grid), ncol=4)
    for (i in c(1:4)) { y[,i] <- sapply(grid, function(x) { b(x, i=i) }) }
    rm(b)

    return(function(x, i) {
      rows <- unlist(lapply(x, function(x) { which.min(abs(x-grid)) } ))
      return(y[rows,i])
    })
    
  } else {
    
    stop("Invalid choices for m and/or s.")
    
  }
}
