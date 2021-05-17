
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
    # testing_probs <- c(0.1,0.1,0.1,0.1,0.1)
    testing_probs <- c(0,0.1,0.2,0.5,1)
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
    
    baseline_status <- ifelse(is.na(d_sero_year[i]), "HIV-", "HIV+ART-")
    
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
