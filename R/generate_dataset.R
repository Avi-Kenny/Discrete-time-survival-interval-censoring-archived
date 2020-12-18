#' Create a dataset resembling AHRI cohort
#'
#' @param num_patients Number of patients in cohort
#' @param start_year Start of cohort (Jan 1st, start_year)
#' @param end_year End of cohort (Jan 1st, end_year)
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
  
  # Generate probabilities for baseline serostatus
  m_probs <- construct_m_probs(p_sero_year)
  
  # Generate baseline variables: patient_id, age, birth_year, sex, u
  # All patients assumed to be born on Jan 1st
  # Age represents number of completed years
  # All baseline variables represent values at Jan 1st of start_year
  {
    d_patient_id <- c(1:num_patients)
    d_age <- sample(1:80, size=num_patients, replace=TRUE) # !!!!! Change this distribution
    d_birth_year <- start_year - d_age
    d_sex <- sample(c(0,1), size=num_patients, replace=TRUE)
    d_u <- rbinom(n=num_patients, size=1, prob=0.3) # !!!!! Set prob=0 in some replicates to remove correlations induced by U
  }
  
  # Generate baseline variables: testing_prob
  # These will be the indices of the testing prob vector c(0,0.1,0.2,0.5,1)
  {
    multinom_probs_u0 <- c(0.2,0.2,0.2,0.2,0.2) # !!!!! Pass this in
    multinom_probs_u1 <- c(0.2,0.2,0.2,0.2,0.2) # !!!!! Pass this in
    d_testing_prob_u0 <- rmultinom(n=num_patients,size=1,prob=multinom_probs_u0)
    d_testing_prob_u1 <- rmultinom(n=num_patients,size=1,prob=multinom_probs_u1)
    d_testing_prob_u0 <- apply(d_testing_prob_u0, 2, function(x){which(x==1)})
    d_testing_prob_u1 <- apply(d_testing_prob_u1, 2, function(x){which(x==1)})
    d_testing_prob_index <- (1-d_u)*d_testing_prob_u0 + d_u*d_testing_prob_u1
    testing_probs <- c(0.1,0.1,0.1,0.1,0.1)
    # testing_probs <- c(0,0.1,0.2,0.5,1)
  }
  
  # Generate baseline seroconversion years
  # !!!!! Incorporate u ?
  d_sero_year <- sapply(d_patient_id, function(i) {
    sex_mf <- ifelse(d_sex[i], "male", "female")
    probs <- m_probs[[sex_mf]][[d_age[i]]]
    pos <- which(as.numeric(rmultinom(n=1, size=1, prob=probs))==1)
    if (pos==length(probs)) {
      return (NA)
    } else {
      return (d_birth_year[i]+pos-1)
    }
  })
  
  # !!!!! Add baseline testing data ?
  
  # Sample cohort events
  dataset <- lapply(d_patient_id, function(i) {
    
    # Set baseline cascade status (one of "HIV-", "HIV+ART-", "HIV+ART+")
    # !!!!! For now, all HIV+ patients are ART- at baseline
    d2_baseline_status <- ifelse(is.na(d_sero_year[i]), "HIV-", "HIV+ART-")
    d2_status <- d2_baseline_status
    
    # Set other baseline variables
    # For baseline population, all HIV+ individuals are assumed to be ART-
    sex_mf <- ifelse(d_sex[i], "male", "female")
    d2_died <- 0
    d2_death_year <- NA
    d2_art_init <- NA
    d2_sero_year <- d_sero_year[i]
    d2_tests <- c()
    
    for (year in c(start_year:(end_year-1))) {
      
      current_age <- year - d_birth_year[i]
      
      # Determine whether patient seroconverted
      if (d2_died==0 && d2_status=="HIV-") {
        p_sero <- p_sero_year[[sex_mf]][current_age+1] * 
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
      
      # !!!!! Just switched order of death and testing
      
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
      
    }
    
    neg_tests <- c(start_year:(end_year-1))[which(d2_tests=="neg")]
    pos_test <- c(start_year:(end_year-1))[which(d2_tests=="pos")]
    tests <- c(neg_tests,pos_test)
    
    d2_first_test <- ifelse(identical(tests,integer(0)), NA, min(tests))
    d2_last_neg_test <- ifelse(identical(neg_tests,integer(0)), NA,
                               max(neg_tests))
    d2_first_pos_test <- ifelse(identical(pos_test,integer(0)), NA,
                                pos_test)
    
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
  
  # Add "testing case" to dataset
  #     Case 1: no testing data
  #     Case 2: most recent test was negative
  #     Case 3: negative test followed by a positive test
  #     Case 4: first test was positive
  dataset %<>% mutate(
    case = case_when(
      is.na(last_neg_test) & is.na(first_pos_test) ~ 1,
      !is.na(last_neg_test) & is.na(first_pos_test) ~ 2,
      !is.na(last_neg_test) & !is.na(first_pos_test) ~ 3,
      is.na(last_neg_test) & !is.na(first_pos_test) ~ 4,
    )
  )
  
  attr(dataset, "start_year") <- start_year
  attr(dataset, "end_year") <- end_year
  
  # QA checks (run manually)
  if (FALSE) {
    
    # This should be zero
    dataset %>% filter(baseline_status!="HIV-" & !is.na(last_neg_test)) %>%
      nrow()

    # !!!!! Add more QA checks
    
  }
  
  return (dataset)
  
}
