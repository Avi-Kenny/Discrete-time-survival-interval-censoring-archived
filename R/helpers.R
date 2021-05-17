#' Expit function
#' 
#' @param x Numeric input
#' @return Numeric output
expit <- function(x) {1/(1+exp(-x))}



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
p_death_year <- function(mult) {
  
  # !!!!! Need to update these numbers
  probs <- c(
    rep(0.01, 9), # 1-9
    rep(0.002, 10), # 10-19
    rep(0.002, 10), # 20-29
    rep(0.004, 10), # 30-39
    rep(0.004, 10), # 40-49
    rep(0.01, 10), # 50-59
    rep(0.01, 10), # 60-69
    rep(0.03, 10), # 70-79
    rep(0.03, 10), # 80-89
    rep(0.1, 10), # 90-99
    rep(0.5, 10) # 100-109
  )
  
  return (mult*probs)
  
}



#' Construct multinomial probabilities
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

construct_m_probs <- function(p_sero_year) {
  
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
