#' Generate baseline data
#'
#' @param num_patients Number of patients in cohort
#' @param start_year Start of cohort (Jan 1st, start_year)
#' @return A dataframe, one row per patient, containing the following fields:
#'     - id: patient ID variable
#'     - b_age: baseline age (in completed years at start_year)
#'     - birth_year: birth year (Jan 1st of year)
#'     - sex: sex (0=female, 1=male)
#'     - u: unmeasured "health behavior" variable
#' @notes
#'     - All dates are in CMC (Century Month Code) format, which is the number
#'       of months since January 1st, 1900

generate_data_baseline <- function(
  num_patients, start_year
) {
  
  # Generate baseline variables: patient_id, age, birth_year, sex, u
  # All patients assumed to be born on Jan 1st
  # Age represents number of completed years
  # All baseline variables represent values at Jan 1st of start_year
  # !!!!! Change the age distribution
  # !!!!! Add baseline status/testing data
  {
    id <- c(1:num_patients)
    b_age <- sample(1:80, size=num_patients, replace=TRUE)
    birth_year <- start_year - b_age
    sex <- sample(c(0,1), size=num_patients, replace=TRUE)
    u <- rnorm(n=num_patients)
  }
  
  dat <- data.frame(id=id, b_age=b_age, birth_year=birth_year, sex=sex, u=u)
  attr(dat, "start_year") <- start_year
  attr(dat, "num_patients") <- num_patients
  
  return (dat)
  
}
