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
  {
    source("helpers.R")
    library(dplyr)
    library(tidyr)
    library(magrittr)
    library(data.table)
    library(survival)
    num_patients <- 100
    start_year <- 2000
    end_year <- 2020
    hazard_ratios <- list("hiv"=1.7, "art"=1.4)
    p_sero_year <- convert_p_sero(list(male = list("1"=0.03, "2-10"=0, "11-15"=0, "16-20"=0.01, "21-25"=0.02,"26-30"=0.03, "31-35"=0.02, "36-40"=0.01, "41-45"=0.005, "46-50"=0.005),female = list("1"=0.03, "2-10"=0, "11-15"=0.005, "16-20"=0.02, "21-25"=0.03,"26-30"=0.02, "31-35"=0.015, "36-40"=0.01, "41-45"=0.005, "46-50"=0)))
    p_death_year <- p_death_year() # !!!!!
    u_mult <- list("sero"=1, "death"=1) # !!!!! U does not induce correlation
    # u_mult = list("sero"=1.5, "death"=1.3)
  }
  # !!!!! TESTING: END !!!!!
  
  # Generate probabilities for baseline serostatus
  m_probs <- construct_m_probs(p_sero_year)
  
  # Generate baseline variables: patient_id, age, birth_year, sex, u
  {
    d_patient_id <- c(1:num_patients)
    d_age <- sample(1:80, size=num_patients, replace=TRUE) # !!!!!
    d_birth_year <- start_year - d_age
    d_sex <- sample(c(0,1), size=num_patients, replace=TRUE)
    d_u <- rbinom(n=num_patients, size=1, prob=0.3) # !!!!! Set prob=0 in some replicates to remove correlations induced by U
  }
  
  # Generate baseline variables: testing_prob
  # These will be the indices of the testing prob vector c(0,0.1,0.2,0.5,1)
  {
    multinom_probs_u0 <- c(0.2,0.2,0.2,0.2,0.2) # !!!!! Pass this in
    multinom_probs_u1 <- c(0.2,0.2,0.2,0.2,0.2) # !!!!! Pass this in
    d_testing_prob_u0 <- rmultinom(n=num_patients, size=1, prob=multinom_probs_u0)
    d_testing_prob_u1 <- rmultinom(n=num_patients, size=1, prob=multinom_probs_u1)
    d_testing_prob_u0 <- apply(d_testing_prob_u0, 2, function(x){which(x==1)})
    d_testing_prob_u1 <- apply(d_testing_prob_u1, 2, function(x){which(x==1)})
    d_testing_prob_index <- (1-d_u)*d_testing_prob_u0 + d_u*d_testing_prob_u1
    testing_probs <- c(0,0.1,0.2,0.5,1)
  }
  
  # !!!!! Baseline testing status ?
  
  # Generate baseline variables: cascade status and seroconversion date
  # !!!!! Incorporate u ?
  d_sero_age <- sapply(d_patient_id, function(i) {
    sex_mf <- ifelse(d_sex[i], "male", "female")
    probs <- m_probs[[sex_mf]][[d_age[i]]] # !!!!! Changed
    sero_age <- which(as.numeric(rmultinom(n=1, size=1, prob=probs))==1)
    if (sero_age==length(probs)) {
      return (NA)
    } else {
      return (sero_age)
    }
  })
  
  d_sero_year <- start_year - d_sero_age
  
  # For baseline population, ART initiation is assumed to happen exactly 12
  #     months after seroconversion
  d_art_init <- d_sero_year + 1
  
  # Sample cohort events
  # !!!!! Can probably speed this up by putting all baseline info into a list of lists and lapplying over that
  dataset <- lapply(d_patient_id, function(i) {
    
    # Set baseline cascade status (one of "HIV-", "HIV+ART-", "HIV+ART+")
    d2_status <- ifelse(is.na(d_sero_year[i]), "HIV-", ifelse(
      d_sero_year[i]==start_year, "HIV+ART-", "HIV+ART+" # !!!!! check the condition
    ))
    d2_baseline_status <- d2_status
    
    # Set other baseline variables
    sex_mf <- ifelse(d_sex[i], "male", "female")
    d2_died <- 0
    d2_death_year <- NA
    d2_art_init <- d_art_init[i]
    d2_sero_year <- d_sero_year[i]
    d2_tests <- c()
    
    for (year in c((start_year+1):end_year)) {
      
      current_age <- year - d_birth_year[i]
      
      # Determine whether patient seroconverted
      if (d2_died==0 && d2_status=="HIV-") {
        p_sero <- p_sero_year[[sex_mf]][current_age] * 
          ifelse(d_u[i],u_mult[["sero"]],1)
        if (p_sero>1) {
          p_sero <- 1
          warning("Probability of seroconversion exceeded one.")
        }
        if (runif(1)<p_sero) {
          # Patient seroconverted
          d2_sero_year <- year
          d2_status <- "HIV+ART-"
        }
      }
      
      # Determine whether patient was tested
      if (d2_died==0 && d2_status != "HIV+ART+") {
        p_test <- testing_probs[d_testing_prob_index[i]]
        if (runif(1)<p_test) {
          # Patient was tested
          if (d2_status == "HIV+ART-") {
            d2_tests <- c(d2_tests, "pos")
            d2_art_init <- year
            d2_status <- "HIV+ART+"
          } else {
            d2_tests <- c(d2_tests, "neg")
          }
        } else {
          # Patient was not tested
          d2_tests <- c(d2_tests, "no")
        }
        
      } else {
        d2_tests <- c(d2_tests, NA)
      }
      
      # Determine whether patient died
      if (d2_died==0) {
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
          # Patient died
          d2_died <- 1
          d2_death_year <- year
        }
      }

    }
    
    neg_tests <- c((start_year+1):end_year)[which(d2_tests=="neg")]
    pos_test <- c((start_year+1):end_year)[which(d2_tests=="pos")]
    tests <- c(neg_tests,pos_test)
    
    d2_first_test <- ifelse(identical(tests,integer(0)), NA, min(tests))
    d2_last_neg_test <- ifelse(identical(neg_tests,integer(0)), NA,
                               max(neg_tests))
    d2_first_pos_test <- ifelse(identical(pos_test,integer(0)), NA,
                                pos_test)
    
    # Assume that baseline HIV+ART+ patients were tested at start_year-1
    if (d2_baseline_status=="HIV+ART+") {
      d2_first_pos_test <- start_year-1
      d2_first_test <- start_year-1
    }
    
    return (list(
      patient_id = d_patient_id[i],
      birth_year = d_birth_year[i],
      baseline_age = d_age[i],
      sex = d_sex[i],
      u = d_u[i],
      baseline_status = d2_baseline_status,
      endline_status = d2_status,
      died = d2_died,
      death_year = d2_death_year,
      art_init = d2_art_init,
      sero_year = d2_sero_year,
      first_test = d2_first_test,
      last_neg_test = d2_last_neg_test,
      first_pos_test = d2_first_pos_test
    ))
    
  })
  
  dataset <- rbindlist(dataset)
  attr(dataset, "start_year") <- start_year
  attr(dataset, "end_year") <- end_year
  
  # QA checks (run manually)
  if (FALSE) {
    
    # This should be zero
    qa_1 <- dataset_cp %>%
      filter(baseline_status!="HIV-" & !is.na(last_neg_test)) %>% nrow()
    if (qa_1!=0) { error("Dataset issue") }
    
    # !!!!! Add more QA checks
    
  }
  
  return (dataset)
  
}
