#' Expit function
#' 
#' @param x Numeric input
#' @return Numeric output
expit <- function(x) {1/(1+exp(-x))}



#' Compute "case" indicator variables
#' 
#' @param x Vector of seroconversion indicators for one individual
#' @param v Vector of testing indicators for one individual
#' @return An integer (1-4) representing the case, as follows:
#'     - Case 1: No testing data
#'     - Case 2: Only NEG tests
#'     - Case 3: NEG test then POS test
#'     - Case 4: first test POS
case <- function(x,v) {
  
  sum_v <- sum(v)
  sum_vx <- sum(v*x)
  sum_v1x <- sum(v*(1-x))
  case_i[1] <- as.integer(sum_v==0)
  case_i[2] <- as.integer(sum_v>0 & sum_vx==0)
  case_i[3] <- as.integer(sum_v>0 & sum_v1x>0 & sum_vx==1)
  case_i[4] <- as.integer(sum_v==1 & sum_vx==1)
  if (sum(case_i)==0) {
    stop("No applicable case")
  }
  if (sum(case_i)>1) {
    stop("Multiple cases matched")
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
T_plusminus <- function(case, J, x, v) {
  
  if (case %in% c(2,3)) {
    # Time of most recent negative test
    T_minus <- which.max(c(1:J)*v*(1-x))
  } else {
    T_minus <- 0
  }
  if (case %in% c(3,4)) {
    # Time of positive test
    T_plus <- which.max(v*x)
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
g_delta <- function(case, J, T_minus, T_plus) {
  
  if (case==1) {
    d <- rep(0, J)
  } else if (case==2) {
    d <- c(rep(1,T_minus), rep(0,J-T_minus))
  } else if (case==3) {
    d <- c(rep(0,T_plus-1), rep(1,J-T_plus+1))
  } else if (case==4) {
    d <- c(rep(1,T_minus), rep(1,T_plus-T_minus-1), rep(1,J-T_plus+1))
  }
  if (length(d)!=J) { stop("Delta is of the wrong length") }
  
  return(d)
  
}


