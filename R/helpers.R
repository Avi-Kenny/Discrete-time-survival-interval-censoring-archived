#' Alias for as.integer
#' 
In <- as.integer

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
  if (!length(v)==len || !length(x)==len) { stop("Incorrect vectors lengths.") }
  
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
