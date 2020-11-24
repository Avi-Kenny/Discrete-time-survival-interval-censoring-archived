#' Perform imputation on a dataset with missingness
#'
#' @param p_sero_year A list of List of monthly seroconversion probabilities
#'     (see simulation constants)
#' @return A list of multinomial probabilities. The index of the list is the age
#'     of an individual in the dataset. For an individual of (exactly) age 33,
#'     the list value at index 33 will be a vector of length 34. This is a
#'     vector of multinomial probabilities, where first entry is the prob that
#'     the individual seroconverted by age 1 (MTCT), the second entry is the
#'     prob that the individual seroconverted by age 2, etc. The last entry is
#'     the prob that the individual never seroconverted.

construct_probs <- function(p_sero_year) {
  
  pm0 <- p_sero_year$male
  pf0 <- p_sero_year$female
  pm <- c(
    pm0[["1"]],
    rep(pm0[["2-10"]],9), rep(pm0[["11-15"]],5), rep(pm0[["16-20"]],5),
    rep(pm0[["21-25"]],5), rep(pm0[["26-30"]],5), rep(pm0[["31-35"]],5),
    rep(pm0[["36-40"]],5), rep(pm0[["41-45"]],5), rep(pm0[["46-50"]],5)
  )
  pf <- c(
    pf0[["1"]],
    rep(pf0[["2-10"]],9), rep(pf0[["11-15"]],5), rep(pf0[["16-20"]],5),
    rep(pf0[["21-25"]],5), rep(pf0[["26-30"]],5), rep(pf0[["31-35"]],5),
    rep(pf0[["36-40"]],5), rep(pf0[["41-45"]],5), rep(pf0[["46-50"]],5)
  )
  
  p_sero_m <- list(c(pm[1],1-pm[1]))
  p_sero_f <- list(c(pf[1],1-pf[1]))
  
  for (age in 2:50) {
    
    # Set next item to previous item (without last prob)
    p_sero_m[[age]] <- p_sero_m[[age-1]][1:(age-1)]
    p_sero_f[[age]] <- p_sero_f[[age-1]][1:(age-1)]
    
    # Calculate and set next prob
    next_prob_m <- prod(1-pm[1:(age-1)]) * pm[age]
    next_prob_f <- prod(1-pf[1:(age-1)]) * pf[age]
    p_sero_m[[age]] <- c(p_sero_m[[age]], next_prob_m)
    p_sero_f[[age]] <- c(p_sero_f[[age]], next_prob_f)
    
    # Set prob of not seroconverting
    p_sero_m[[age]] <- c(p_sero_m[[age]], 1-sum(p_sero_m[[age]]))
    p_sero_f[[age]] <- c(p_sero_f[[age]], 1-sum(p_sero_f[[age]]))
    
  }
  
  return(list(
    male = p_sero_m,
    female = p_sero_f
  ))
  
}
