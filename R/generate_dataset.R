#' Create a dataset resembling AHRI cohort
#'
#' @param num_patients Number of patients in cohort
#' @param start_year Year of cohort enrollment
#' @param end_year Last date at which we have cohort patient data
#' @param hazard_ratios List of hazard ratios for HIV+/ART- status and
#'     HIV+/ART+ status (relative to HIV-)
#' @param p_sero_year List of yearly discrete hazards of seroconversion, by sex
#' @param p_death_year List of yearly discrete hazards of death
#' @param u_mult List of hazard multipliers
#' @return A dataframe containing the following:
#'     !!!!! Update this
#'     patient_id: patient ID number
#'     dob: date of birth
#'     sex: sex (0=female, 1=male)
#'     alive: patient is alive (0=no, 1=yes)
#'     dod: date of death
#'     last_test_neg: date of most recent negative HIV test
#'     first_test_pos: date of first positive HIV test
#'     art_init: ART initiation date
#'     s: seroconversion date
#'     case: case code (see concept note appendix)
#'     hiv_status: HIV status code (see concept note appendix)
#'     stroke_date: date of stroke event (NA if no stroke occurred)
#' @notes
#'     - All dates are in CMC (Century Month Code) format, which is the number
#'       of months since January 1st, 1900

generate_dataset <- function(
  num_patients, start_year, end_year, hazard_ratios, p_sero_year, p_death_year,
  u_mult
) {
  
  # !!!!! TESTING: START !!!!!
  source("helpers.R")
  library(dplyr)
  library(tidyr)
  library(magrittr)
  library(data.table)
  num_patients = 10000
  start_year = 2000
  end_year = 2020
  hazard_ratios = list("hiv"=1.7, "art"=1.4)
  p_sero_year <- convert_p_sero(list(male = list("1"=0.03, "2-10"=0, "11-15"=0, "16-20"=0.01, "21-25"=0.02,"26-30"=0.03, "31-35"=0.02, "36-40"=0.01, "41-45"=0.005, "46-50"=0.005),female = list("1"=0.03, "2-10"=0, "11-15"=0.005, "16-20"=0.02, "21-25"=0.03,"26-30"=0.02, "31-35"=0.015, "36-40"=0.01, "41-45"=0.005, "46-50"=0)))
  p_death_year <- p_death_year() # !!!!!
  u_mult = list("sero"=1.5, "death"=1.3)
  # !!!!! TESTING: END !!!!!
  
  # !!!!! Remove unused variables
  
  # Generate probabilities for baseline serostatus
  m_probs = construct_probs(p_sero_year)
  
  # Generate baseline variables: patient_id, age, birthyear, sex, u
  d_patient_id <- c(1:num_patients)
  d_age <- sample(1:80, size=num_patients, replace=TRUE) # !!!!!
  d_birthyear <- start_year - d_age
  d_sex <- sample(c(0,1), size=num_patients, replace=TRUE)
  d_u <- rbinom(n=num_patients, size=1, prob=0.3) # !!!!! Set prob=0 in some replicates to remove correlations induced by U
  
  # Generate baseline variables: testing_prob
  # These will be the indices of the testing prob vector c(0,0.1,0.2,0.5,1)
  multinom_probs_u0 <- c(0.2,0.2,0.2,0.2,0.2) # !!!!! Pass this in
  multinom_probs_u1 <- c(0.2,0.2,0.2,0.2,0.2) # !!!!! Pass this in
  d_testing_prob_u0 <- rmultinom(n=num_patients, size=1, prob=multinom_probs_u0)
  d_testing_prob_u1 <- rmultinom(n=num_patients, size=1, prob=multinom_probs_u1)
  d_testing_prob_u0 <- apply(d_testing_prob_u0, 2, function(x){which(x==1)})
  d_testing_prob_u1 <- apply(d_testing_prob_u1, 2, function(x){which(x==1)})
  d_testing_prob_index <- (1-d_u)*d_testing_prob_u0 + d_u*d_testing_prob_u1
  testing_probs <- c(0,0.1,0.2,0.5,1)
  
  # !!!!! Baseline testing status ?
  
  # Generate baseline variables: cascade status and seroconversion date
  d_sero_age <- sapply(d_patient_id, function(i) {
    sex_mf <- ifelse(d_sex[i], "male", "female")
    probs <- m_probs[[sex_mf]][[min(d_age[i],50)]]
    pos <- which(as.numeric(rmultinom(n=1, size=1, prob=probs))==1) # !!!!! Incorporate u ?????
    if (pos==length(probs)) {
      return (NA)
    } else {
      return (pos)
    }
  })
  
  d_sero_year <- start_year - d_sero_age
  
  # For now, ART initiation is assumed to happen exactly 12 months after
  #     seroconversion
  d_art_init <- d_sero_year + 1
  
  # Sample cohort events
  df <- lapply(d_patient_id, function(i) {
    
    # Set baseline cascade status (one of "HIV-", "HIV+ART-", "HIV+ART+")
    d2_status <- ifelse(is.na(d_sero_year[i]), "HIV-", ifelse(
      d_sero_year[i]==start_year, "HIV+ART-", "HIV+ART+" # !!!!! check the condition
    ))
    
    # Set other baseline variables
    sex_mf <- ifelse(d_sex[i], "male", "female")
    d2_died <- 0
    d2_death_year <- NA
    d2_art_init <- d_art_init[i]
    d2_sero_age <- d_sero_age[i]
    d2_sero_year <- d_sero_year[i]
    d2_tests <- c()
    
    for (year in c((start_year+1):end_year)) {
      
      current_age <- year - d_birthyear[i]
      
      if (is.na(d2_death_year)) {
        
        # !!!!! ART starts one year after seroconversion
        # !!!!! Change this to "patient starts ART after first positive test"
        if (d2_status=="HIV+ART-") {
          d2_status <- "HIV+ART+"
          d2_art_init <- year
        }
        
        # Determine whether patient was tested
        p_test <- testing_probs[d_testing_prob_index[d_patient_id]]
        if (runif(1)<p_test) {
          d2_tests <- c(d2_tests, year)
        }
        
        # Determine whether patient seroconverted
        if (d2_status=="HIV-") {
          p_sero <- p_sero_year[[sex_mf]][current_age] * 
            ifelse(d_u[i],u_mult[["sero"]],1)
          if (p_sero>1) {
            p_sero <- 1
            warning("Probability of death exceeded one.")
          }
          if (runif(1)<p_sero) {
            d2_sero_age <- current_age
            d2_sero_year <- year
            d2_status <- "HIV+ART-"
          }
        }
        
        # Determine whether patient died
        hr <- ifelse(
          d2_status=="HIV+ART-", hazard_ratios$hiv,
            ifelse(d2_status=="HIV+ART+", hazard_ratios$art, 1)
        )
        p_death <- p_death_year[current_age] * hr *
          ifelse(d_u[i],u_mult[["death"]],1)
        if (p_death>1) {
          p_death <- 1
          warning("Probability of death exceeded one.")
        }
        
        if (runif(1)<p_death) {
          d2_died <- 1
          d2_death_year <- year
        }
        
      }
      
    }
    
    return (list(
      patient_id = d_patient_id[i],
      birthyear = d_birthyear[i],
      baseline_age = d_age[i],
      sex = d_sex[i],
      u = d_u[i],
      baseline_status = d2_status,
      died = d2_died,
      death_year = d2_death_year,
      art_init = d2_art_init,
      sero_age = d2_sero_age,
      sero_year = d2_sero_year,
      tests = d2_tests
    ))
    
  })
  
  df2 <- lapply(df, function(l) {
    l$tests <- NULL
    return(l)
  })
  
  df3 <- lapply(df, function(l) {
    return(l$tests)
  })
  
  df_patients <- rbindlist(df2)
  df_testing <- rbindlist(df3)
  
  attr(df, "start_year") <- start_year
  attr(df, "end_year") <- end_year
  
  
  
  # !!!!! Generate testing data
  
  # testing_case <- which(rmultinom(n=1, size=1, prob=c(0.1,0.3,0.6))==1)
  
  df %<>% mutate(
    
  )
  
  return (df)
  
  # testing case (case 1, case 2, etc.)
  
  # `testing_probs` is 
  testing_probs <- 
  # testing dates (last_test_neg, first_test_pos)
  
  

  return (data.frame(
    patient_id = d_patient_id
  ))

}
