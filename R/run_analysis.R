#' Run Cox PH analysis
#'
#' @param dataset A dataset returned by either generate_dataset() or
#'     perform_imputation()
#' @return A list containing the following:
#'     est: point estimate of exposure coefficient
#'     se: standard error of exposure coefficient estimate

run_analysis <- function(dataset, start_year, end_year) {
  
  # Transform dataset
  {
    # dataset <- df # !!!!!
    
    # Data cleaning
    dataset$start_year <- attributes(dataset)$start_year
    dataset$death_year <- replace_na(dataset$death_year,
                                     attributes(dataset)$end_year)
    dataset %<>% rename("end_year" = death_year)
    
    # Put dataset into "counting process" format
    dataset_cp <- survSplit(
      formula = Surv(start_year, end_year, died) ~.,
      data = dataset,
      cut = c(start_year:end_year)
    )
    
    # Create time-varying exposure variable
    dataset_cp %<>% mutate(
      art_status = replace_na(ifelse(start_year>=art_init, 1, 0),0),
      age = start_year - birthyear
    )
    
    # Create time-varying age bucket variable
    dataset_cp %<>% mutate(
      age_bin = case_when(
        age %in% c(1:9) ~ 0,
        age %in% c(10:19) ~ 1,
        age %in% c(20:29) ~ 2,
        age %in% c(30:39) ~ 3,
        age %in% c(40:49) ~ 4,
        age %in% c(50:59) ~ 5,
        age %in% c(60:69) ~ 6,
        age %in% c(70:79) ~ 7,
        age %in% c(80:89) ~ 8,
        age %in% c(90:99) ~ 9,
        age %in% c(100:109) ~ 10
      )
    )
    
  }
  
  # Fit time-varying Cox model
  fit <- coxph(
    Surv(start_year, end_year, died) ~ art_status + factor(age_bin) +
                                       cluster(patient_id),
    data = dataset_cp
  )
  
  # fit <- coxph(
  #   Surv(start_time, end_time, had_event)~factor(hiv_status),
  #   data = dataset
  # )
  summ <- summary(fit)$coefficients
  results <- list(
    exp_est = exp(summ["art_status","coef"]),
    est = summ["art_status","coef"],
    se = summ["art_status","se(coef)"]
  )
  print(results) # !!!!!
  
  return(results)
  
}
