#' Create a dataset resembling AHRI cohort
#'
#' @param num_patients Number of patients in cohort
#' @param start_year Year of cohort enrollment
#' @param end_year Last date at which we have cohort patient data
#' @param hr_hiv Hazard ratio of HIV+ ART- status (relative to HIV-)
#' @param hr_art Hazard ratio of HIV+ ART+ status (relative to HIV-)
#' @param m_probs List of multinomial probabilities returned by construct_probs()
#' @return A dataframe containing the following:
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
  num_patients, start_year, end_year, hr_hiv, hr_art, m_probs
) {
  
  # !!!!! Testing
  {
    p_sero_year <- list(male = list("1"=0.03, "2-10"=0, "11-15"=0, "16-20"=0.01, "21-25"=0.02,"26-30"=0.03, "31-35"=0.02, "36-40"=0.01, "41-45"=0.005, "46-50"=0.005),female = list("1"=0.03, "2-10"=0, "11-15"=0.005, "16-20"=0.02, "21-25"=0.03,"26-30"=0.02, "31-35"=0.015, "36-40"=0.01, "41-45"=0.005, "46-50"=0))
    num_patients <- 5
    start_year <- 2000
    end_year <- 2020
    hr_hiv <- 1.7
    hr_art <- 1.2
    m_probs <- construct_probs(p_sero_year)
  }
  
  # Generate patient_id, dob, sex
  # Assumes that sex is independent of DOB
  # !!!!! d_dob is placeholder; replace with actual distribution of patient ages
  d_patient_id <- c(1:num_patients)
  d_age <- sample(1:100, size=num_patients, replace=TRUE)
  d_birthyear <- end_year - d_age
  d_sex <- sample(c(0,1), size=num_patients, replace=TRUE)
  
  # Sample seroconversion age from multinomial
  d_sero_age <- sapply(d_patient_id, function(i) {
    age <- d_age[i]
    sex <- ifelse(d_sex[i], "male", "female")
    probs <- m_probs[[sex]][[min(age,50)]]
    mult <- as.numeric(rmultinom(n=1, size=1, prob=probs))
    pos <- which(mult==1)
    if (pos==length(probs)) {
      return (NA)
    } else {
      return (pos)
    }
  })
  
  d_sero_year <- end_year - d_sero_age
  
  # For now, ART initiation is assumed to happen exactly 12 months after
  #     seroconversion
  d_art_init <- d_sero_year + 1
  
  # alive + dod // this is our outcome for now; depends on cascade status
  d_dods <- rep(NA, num_patients)
  for (year in c((start_year+1):end_year)) {
    
    # !!!!! CONTINUE HERE
    
    p_death <- 0.05 # !!!!! should be a function of age/sex
    p_death_hiv <- hr_hiv * p_death
    p_death_art <- hr_art * p_death
    
    # probs <- 
    
    unifs <- runif(num_patients)
    # deaths <- 
    
  }
  
  
  # testing case (case 1, case 2, etc.)
  # testing dates (last_test_neg, first_test_pos)
  
  
  
  
  
  
  # !!!!! Recycle code below this point
  
  # Populate data frame
  # !!!!! Vectorize this if possible
  for (i in 1:num_patients) {
    
    # Calculate `date of death`
    if (!alive) {
      dod <- ifelse(
        dob != end_date,
        sample(dob:end_date, size=1),
        dob
      )
    } else {
      dod <- NA
    }
    
    # Calculate `case type`
    case <- sample(c(1,2,3,5), size=1)
    case <- ifelse(case==3 && !is.na(dod),4,case)
    
    # Calculate `date of last HIV- test`, `date of first HIV+ test`
    last_date <- ifelse(is.na(dod),end_date,dod)
    switch(case,
           
           # Case 1: First HIV test was pos
           "1" = {
             last_test_neg <- NA
             first_test_pos <- ifelse(
               dob != last_date,
               sample(dob:last_date, size=1),
               dob
             )
             art_init <- ifelse(runif(1)<0.5,
                                ifelse(
                                  first_test_pos != last_date,
                                  sample(first_test_pos:last_date, size=1),
                                  first_test_pos
                                ),
                                NA
             )
           },
           
           # Case 2: 1+ neg tests followed by a pos test
           "2" = {
             last_test_neg <- ifelse(
               dob != last_date,
               sample(dob:last_date, size=1),
               dob
             )
             first_test_pos <- ifelse(
               last_test_neg != last_date,
               sample(last_test_neg:last_date, size=1),
               last_test_neg
             )
             art_init <- ifelse(runif(1)<0.5,
                                ifelse(
                                  first_test_pos != last_date,
                                  sample(first_test_pos:last_date, size=1),
                                  first_test_pos
                                ),
                                NA
             )
           },
           
           # Case 3: 1+ neg tests and no pos test (and is alive)
           "3" = {
             last_test_neg <- ifelse(
               dob != last_date,
               sample(dob:last_date, size=1),
               dob
             )
             first_test_pos <- NA
             art_init <- NA
           },
           
           # Case 4: 1+ neg tests and no pos test (and is dead)
           "4" = {
             last_test_neg <- ifelse(
               dob != last_date,
               sample(dob:last_date, size=1),
               dob
             )
             first_test_pos <- NA
             art_init <- NA
           },
           
           # Case 5: no testing data
           "5" = {
             last_test_neg <- NA
             first_test_pos <- NA
             art_init <- NA
           }
           
    )
    
    # Add row to data frame
    df[i,] <- list(i, dob, sex, alive, dod, last_test_neg,
                   first_test_pos, art_init, NA, case)
    
  }
  
  # Generate "true" seroconversion dates
  df$s <- impute_sero_dates(
    df = df,
    end_date = (2019-1900)*12
  )
  
  return (data.frame(
    patient_id = d_patient_id,
    age = d_age,
    birthyear = d_birthyear,
    sex = d_sex,
  ))

}
