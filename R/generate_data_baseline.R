#' Generate baseline data
#'
#' @param num_patients Number of patients in cohort
#' @param start_year Start of cohort (Jan 1st, start_year)
#' @return !!!!!
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
  {
    id <- c(1:num_patients)
    age <- sample(1:80, size=num_patients, replace=TRUE) # !!!!! Change this distribution
    birth_year <- start_year - d_age
    sex <- sample(c(0,1), size=num_patients, replace=TRUE)
    u <- rnorm(n=num_patients)
  }
  
  return (data.frame(id=id, age=age, birth_year=birth_year, sex=sex, u=u))
  
}
