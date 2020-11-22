#' Create a dataset resembling AHRI cohort
#'
#' @param num_patients Number of patients in cohort
#' @param end_date !!!!! TO DO
#' @param rel_risk !!!!! incorporate
#' @param missingness !!!!! incorporate
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
#' @todo
#'     - !!!!! Incorporate rel_risk and missingness parameters (both currently unused)
#'     - Add other predictor variables (e.g. SES)
#'     - Ensure marginal or joint distributions of all variables matches the AHRI cohort
#'     - alive and dod should incorporate rel_risk
#'     - case should match proportion of case types in cohort
#'     - last_test_neg, first_test_pos should incorporate missingness and reflect cohort testing rates

generate_dataset <- function(num_patients, end_date, rel_risk, missingness) {
  
  # !!!!! Testing
  end_date <- 1440
  num_patients <- 5
  
  # Generate patient_id, dob, sex
  # Assumes that sex is independent of DOB
  df <- data.frame(
    patient_id = c(1:num_patients),
    dob = sample((end_date-(12*100)):end_date, size=num_patients, replace=TRUE),
    sex = sample(c(0,1), size=num_patients, replace=TRUE)
  )
  
  # !!!!! Add to df:
    # "alive" = integer(),
    # "dod" = integer(),
    # "last_test_neg" = integer(),
    # "first_test_pos" = integer(),
    # "art_init" = integer(),
    # "s" = integer(),
    # "case" = integer(),
    # "hiv_status" = integer(),
    # "stroke_date" = integer()
    # !!!!! add predictors
  
  
  
  
  
  
  # !!!!! Recycle code below this point
  
  # Populate data frame
  # !!!!! Vectorize this if possible
  for (i in 1:num_patients) {
    
    # Calculate `patient is alive`
    alive <- rbinom(n=1, size=1, prob=0.9)
    
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
  
  return (df)

}
