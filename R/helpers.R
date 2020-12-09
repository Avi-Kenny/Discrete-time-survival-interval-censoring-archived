# #' Expit function
# #' @param x number
# #' @return Expit of x
# 
# expit <- function(x) { exp(x) / (1+exp(x)) }



#' Format converter for p_sero_year
#' @param p_sero_year number
#' @return p_sero_year, but with individual years instead of buckets
#' 
convert_p_sero <- function(p_sero_year) {
  
  new_list <- list()
  for (sex in c("male", "female")) {
    p <- p_sero_year[[sex]]
    new_probs <- c(
      rep(p[["1"]],1), rep(p[["2-10"]],9), rep(p[["11-15"]],5),
      rep(p[["16-20"]],5), rep(p[["21-25"]],5), rep(p[["26-30"]],5),
      rep(p[["31-35"]],5), rep(p[["36-40"]],5), rep(p[["41-45"]],5),
      rep(p[["46-50"]],5), rep(0, 50)
    )
    new_list[[sex]] <- new_probs
  }
  
  return (new_list)
  
}



#' Return discrete hazards of death, by age
#' @return A vector of discrete hazards, indexed by age
#' 
p_death_year <- function() {
  
  # !!!!! Need to update these numbers
  # !!!!! Make this a function of age as well?
  
  return (c(
    0.02, rep(0.01, 9), # 1-10
    rep(0.004, 10), # 11-20
    rep(0.004, 10), # 21-30
    rep(0.01, 10), # 31-40
    rep(0.02, 10), # 41-50
    rep(0.04, 10), # 51-60
    rep(0.05, 10), # 61-70
    rep(0.1, 10), # 71-80
    rep(0.2, 10), # 81-90
    rep(0.3, 9), # 91-99
    1
  ))
  
}



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
  
  pm <- p_sero_year$male
  pf <- p_sero_year$female
  
  p_sero_m <- list(c(pm[1],1-pm[1]))
  p_sero_f <- list(c(pf[1],1-pf[1]))
  
  for (age in 2:100) {
    
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
