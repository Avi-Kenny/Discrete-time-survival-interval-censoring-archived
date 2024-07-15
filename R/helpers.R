#' Alias for as.integer
#' 
In <- as.integer

#' Expit function
#' 
#' @param x Numeric input
#' @return Numeric output
expit <- function(x) {1/(1+exp(-x))}



#' #' Modified exp function (see scratch for derivation)
#' #' 
#' #' @param x Numeric input
#' #' @return Numeric output
#' exp2 <- (function() {
#'   
#'   expit <- function(x) {1/(1+exp(-x))}
#'   logit <- function(x) { log(x/(1-x)) }
#'   e <- -0.1 # This value configurable but hard-coded
#'   ell <- logit(exp(e))
#'   x_0 <- e - (ell*exp(ell))/(exp(e)*(1+exp(ell))^2)
#'   k_0 <- exp(e-ell)*(1+exp(ell))^2
#'   exp2 <- function(x) {
#'     if (x<=e) {
#'       return(exp(x))
#'     } else {
#'       return(1/(1+exp(k_0*(x_0-x))))
#'     }
#'   }
#'   return(exp2)
#' })()



#' Inverse of complementary log-log link function
#' 
#' @param x Numeric input
#' @return Numeric output
icll <- function(x) { 1 - exp(-exp(x)) }



#' Compute "case" indicator variables
#' 
#' @param T_minus Calendar time of first negative test (0 = no NEG tests)
#' @param T_plus Calendar time of first positive test (0 = no POS tests)
#' @return An integer (1-4) representing the case, as follows:
#'     - Case 1: No testing data
#'     - Case 2: Only NEG tests
#'     - Case 3: NEG test then POS test
#'     - Case 4: first test POS
case <- function(T_minus, T_plus) {
  
  if (T_minus==0) {
    if (T_plus==0) { return(1) } else { return(4) }
  } else {
    if (T_plus==0) { return(2) } else { return(3) }
  }
  
  # return(dplyr::case_when(
  #   T_minus==0 & T_plus==0 ~ 1,
  #   T_minus!=0 & T_plus==0 ~ 2,
  #   T_minus!=0 & T_plus!=0 ~ 3,
  #   T_minus==0 & T_plus!=0 ~ 4,
  #   TRUE ~ 999
  # ))
  
}



#' Calculate time(s) of most recent negative test and/or positive test
#' 
#' @param J Maximum time
#' @param x Vector of seroconversion indicators for one individual
#' @param v Vector of testing indicators for one individual
#' @return A list containing T_minus and T_plus (0 = no NEG/POS tests)
T_plusminus <- function(s_i, t_i, x, v) {
  
  len <- t_i-s_i+1
  if (!length(v)==len || !length(x)==len) { stop("Incorrect vector lengths.") }
  
  some_tests_neg <- In(sum(v*(1-x))>0)
  some_tests_pos <- In(sum(v*x)>0)
  
  if (some_tests_neg) {
    # Time of most recent negative test
    T_minus_index <- which.max(c(s_i:t_i)*v*(1-x))
    T_minus <- c(s_i:t_i)[T_minus_index]
  } else {
    T_minus <- 0
  }
  if (some_tests_pos) {
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



#' Create a vector of ones and zeros
#' 
#' @param len Length of vector to output
#' @param num_ones Number of "1"s in the tail
uncompress <- function(len, num_ones) {
  c(rep(0,len-num_ones), rep(1,num_ones))
}



#' Create a natural cubic spline basis
#' 
#' @param which Which basis to construct
#' @param window_start Start year
construct_basis <- function(which, window_start=NA) {
  
  scale_age <- function(x) { x / 100 }
  scale_year <- function(x) { (x-window_start+1)/10 }
  
  if (which=="age (0-100), 4DF") {
    grid <- scale_age(seq(0,100, length.out=500))
    k <- scale_age(c(0,25,50,75,100))
  } else if (which=="age (13,20,30,60,90)") {
    grid <- scale_age(seq(13,90, length.out=500))
    k <- scale_age(c(13, 20, 30, 60, 90))
  } else if (which=="age (13,30,60,75,90)") {
    grid <- scale_age(seq(13,90, length.out=500))
    k <- scale_age(c(13, 30, 60, 75, 90))
  } else if (which=="age (13,32,52,71,90)") {
    grid <- scale_age(seq(13,90, length.out=500))
    k <- scale_age(seq(13,90, length.out=5))
  } else if (which %in% c("age (13,28,44,60,75)", "age (13,28,44,60,75) +i")) {
    grid <- scale_age(seq(13,75, length.out=500))
    k <- scale_age(round(seq(13,75, length.out=5)))
  } else if (which=="year (00,05,10,15,20)") {
    grid <- scale_year(seq(2000,2022, length.out=500))
    k <- scale_year(seq(2000,2020, length.out=5))
  } else if (which %in% c("year (10,13,17,20,23)", "year (10,13,17,20,23) +i")) {
    grid <- scale_year(seq(2010,2023, length.out=500))
    k <- scale_year(seq(2010,2023, length.out=5))
  } else if (which=="age (13,20,30,40,60)") {
    grid <- scale_age(seq(13,60, length.out=500))
    k <- scale_age(c(13,20,30,40,60))
  }
  
  if (substr(which, nchar(which)-1, nchar(which))=="+i") {
    num_df <- 5
    int <- TRUE
  } else {
    num_df <- 4
    int <- FALSE
  }
  
  b <- Vectorize(function(x, i) {
    splines::ns(x=x, knots=k[2:4], intercept=int, Boundary.knots=k[c(1,5)])[i]
  }, vectorize.args="x")
  y <- matrix(NA, nrow=length(grid), ncol=num_df)
  for (i in c(1:num_df)) { y[,i] <- sapply(grid, function(x) { b(x, i=i) }) }
  rm(b)
  
  return(function(x=NA, i=NA) {
    row <- unlist(lapply(x, function(x) { which.min(abs(x-grid)) } ))
    if (is.na(i)) { return(y[row,]) } else { return(y[row,i]) }
  })
  
}
