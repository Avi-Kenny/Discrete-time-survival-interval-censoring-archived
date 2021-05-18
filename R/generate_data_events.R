#' Generate cohort events (seroconversion, ART initiation, testing, death)
#'
#' @param b_age Age of patient at start_year
#' @param sex Sex of patient (0=female,1=male)
#' @param start_year Start of cohort (Jan 1st, start_year)
#' @param end_year End of cohort (Jan 1st, end_year)
#' @param baseline_status Baseline cascade status; one of c("HIV-","HIV+ART-",
#'     "HIV+ART+")
#' @param params A list of the form list(alpha0=1,alpha1=1,alpha2=1,alpha3=1,
#'     beta0=1,beta1=1,beta2=1,beta3=1,eta0=1,eta1=1,eta2=1,eta3=1,gamma0=1,
#'     gamma1=1,gamma2=1,gamma3=1,psi1=1,psi2=1)
#' @return A list containing the following:
#'     - v: vector of testing indicators
#'     - x: vector of serostatus indicators
#'     - y: vector of outcome indicators
#'     - z: vector of ART status indicators
#'     - J: vector of outcome indicators
#' @notes
#'     - Parameters correspond to the Markov model

generate_data_events <- function(
  b_age, sex, u, start_year, end_year, baseline_status, params
) {
  
  p <- params
  
  # Set baseline variables
  # !!!!! For now, all patients are HIV- at baseline
  x <- v <- z <- y <- c()
  x_last <- 0
  z_last <- 0
  
  # Sample events
  # Note: this code mirrors the JAGS code
  for (j in 1:(12*(end_year-start_year))) {
    
    if (length(y)==0 || max(y, na.rm=TRUE)==0) {
      
      # Seroconversion
      p_sero <- ifelse(x_last==1, 1, expit(
        p$alpha0 + p$alpha1*sex + p$alpha2*(b_age+(j-1)/12) + p$alpha3*u
      ))
      x <- c(x, rbinom(n=1, size=1, prob=p_sero))
      
      # Testing
      p_test <- expit(
        p$beta0 + p$beta1*sex + p$beta2*(b_age+(j-1)/12) + p$beta3*u
      )
      v <- c(v, rbinom(n=1, size=1, prob=p_test))
      
      # ART
      p_art <- ifelse(z_last==1, 1,
        ifelse(x[length(x)]==0 || v[length(v)]==0, 0, expit(
          p$eta0 + p$eta1*sex + p$eta2*(b_age+(j-1)/12) + p$eta3*u
        ))
      )
      z <- c(z, rbinom(n=1, size=1, prob=p_art))
      
      # Outcome
      p_y <- min(0.99999, expit(
        p$gamma0 + p$gamma1*sex + p$gamma2*(b_age+(j-1)/12) + p$gamma3*u
      ) * exp(
        log(p$psi1)*x[length(x)]*(1-z[length(z)]) +
        log(p$psi2)*x[length(x)]*z[length(z)]
      ))
      y <- c(y, rbinom(n=1, size=1, prob=p_y))
      
      x_last <- x[length(x)]
      z_last <- z[length(z)]
      
    } else {
      
      v <- c(v, NA)
      x <- c(x, NA)
      y <- c(y, NA)
      z <- c(z, NA)
      
    }
    
  }
  
  dat <- list(v=v, x=x, y=y, z=z, J=sum(!is.na(y)))
  
  return (dat)
  
  # # Add "testing case" to dataset
  # #     Case 1: no testing data
  # #     Case 2: most recent test was negative
  # #     Case 3: negative test followed by a positive test
  # #     Case 4: first test was positive
  # dataset %<>% mutate(
  #   case = case_when(
  #     is.na(last_neg_test) & is.na(first_pos_test) ~ 1,
  #     !is.na(last_neg_test) & is.na(first_pos_test) ~ 2,
  #     !is.na(last_neg_test) & !is.na(first_pos_test) ~ 3,
  #     is.na(last_neg_test) & !is.na(first_pos_test) ~ 4,
  #   )
  # )
  
}
