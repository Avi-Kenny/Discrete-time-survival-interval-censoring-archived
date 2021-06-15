#' Generate cohort events (seroconversion, ART initiation, testing, death)
#'
#' @param id Patient ID
#' @param b_age Age of patient at start_year
#' @param sex Sex of patient (0=female,1=male)
#' @param u Latent health behavior variable
#' @param start_year Start of cohort (Jan 1st, start_year)
#' @param end_year End of cohort (Jan 1st, end_year)
#' @param baseline_status Baseline cascade status; one of c("HIV-","HIV+ART-",
#'     "HIV+ART+"); currently unused !!!!!
#' @param params A list of Markov model parameters (alpha, beta, ...)
#' @return A list containing the following:
#'     - v: vector of testing indicators
#'     - x: vector of serostatus indicators
#'     - y: vector of outcome indicators
#'     - z: vector of ART status indicators
#'     - J: vector of outcome indicators
#' @notes Much of this code mirrors code in fit_stan.R; ensure the two are in
#'     sync with one another

generate_data_events <- function(
  id, b_age, sex, u, start_year, end_year, baseline_status, params
) {
  
  p <- params
  
  # Set baseline variables
  x <- v <- z <- y <- c()
  x_last <- z_last <- 0
  
  # Sample events
  # Note: this code mirrors the MCMC code
  for (t in 1:(12*(end_year-start_year))) {
    
    if (length(y)==0 || max(y, na.rm=TRUE)==0) {
      
      # Seroconversion
      p_sero <- ifelse(x_last==1, 1, expit(
        p$alpha0 + p$alpha1*sex + p$alpha2*(b_age+(t-1)/12) + p$alpha3*u
      ))
      x <- c(x, rbinom(n=1, size=1, prob=p_sero))
      
      # Testing
      # !!!!! Add a condition s.t. patient doesn't get tested after POS test
      p_test <- expit(
        p$beta0 + p$beta1*sex + p$beta2*(b_age+(t-1)/12) + p$beta3*u
      )
      v <- c(v, rbinom(n=1, size=1, prob=p_test))
      
      # ART
      p_art <- ifelse(z_last==1, 1,
        ifelse(x[length(x)]==0 || v[length(v)]==0, 0, expit(
          p$eta0 + p$eta1*sex + p$eta2*(b_age+(t-1)/12) + p$eta3*u
        ))
      )
      z <- c(z, rbinom(n=1, size=1, prob=p_art))
      
      # Outcome
      p_y <- min(0.99999, expit(
        p$gamma0 + p$gamma1*sex + p$gamma2*(b_age+(t-1)/12) + p$gamma3*u
      ) * exp(
        log(p$psi1)*x[length(x)]*(1-z[length(z)]) +
        log(p$psi2)*x[length(x)]*z[length(z)]
      ))
      y <- c(y, rbinom(n=1, size=1, prob=p_y))
      
      x_last <- x[length(x)]
      z_last <- z[length(z)]
      
    } else {
      
      # NA values coded as 9 (for Stan)
      v <- c(v, 9)
      x <- c(x, 9)
      y <- c(y, 9)
      z <- c(z, 9)
      
    }
    
  }
  
  # Add "testing case" to dataset
  #     Case 1: no testing data
  #     Case 2: most recent test was negative
  #     Case 3: negative test followed by a positive test
  #     Case 4: first test was positive
  T_i <- sum(v!=9)
  s <- as.integer(sum(v[1:T_i])>0)
  test_first <- s * min( (1:T_i) + T_i*(1-v[1:T_i]) )
  test_last <- s * max( (1:T_i)*v[1:T_i] )
  case <- ifelse(s==0, 1, ifelse(
    x[test_last]==0, 2, ifelse(
      x[test_first]==1, 4, 3
    )
  ))
  
  # Add last_neg_test and first_pos_test
  last_neg_test <- 0
  first_pos_test <- 0
  if (case==2) {
    last_neg_test <- test_last
  }
  if (case==3) {
    last_neg_test <- max( (1:T_i)*(v[1:T_i])*(1-x[1:T_i]) )
    first_pos_test <- min( (1:T_i) + T_i*(1-(x[1:T_i])*(v[1:T_i])) )
  }
  if (case==4) {
    first_pos_test <- test_first
  }
  
  # Calculate delta and delta*x
  miss <- rep(9,sum(x==9))
  if (case==1) {
    delta <- c(rep(0,T_i), miss)
  }
  if (case==2) {
    delta <- c(rep(1,last_neg_test), rep(0,T_i-last_neg_test), miss)
  }
  if (case==3) {
    delta <- c(rep(1,last_neg_test),
               rep(0,first_pos_test-last_neg_test-1),
               rep(1,T_i-first_pos_test+1),
               miss)
  }
  if (case==4) {
    delta <- c(rep(0,first_pos_test-1), rep(1,T_i-first_pos_test+1), miss)
  }
  deltax <- delta*x
  deltax <- ifelse(deltax==81,9,deltax)
  
  # !!!!! Condense code when porting to Stan
  
  return(list(id=id, v=v, x=x, y=y, z=z, T_i=T_i, last_neg_test=last_neg_test,
              first_pos_test=first_pos_test, test_first=test_first,
              test_last=test_last, case=case, delta=delta, deltax=deltax))
  
}
