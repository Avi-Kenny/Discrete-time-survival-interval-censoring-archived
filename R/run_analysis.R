#' Run Cox PH analysis
#'
#' @param dataset_cp A dataset returned by transform_dataset()
#' @param method One of c("ideal", "censor", "mi")
#' @param include_no_testers Binary; whether to include individuals in the
#'     analysis who have no testing data; this is always TRUE for
#'     method=="censor"
#' @return A list containing the following:
#'     est_hiv: point estimate of HIV+ART- exposure coefficient
#'     se_hiv: standard error of HIV+ART- exposure coefficient
#'     est_art: point estimate of HIV+ART+ exposure coefficient
#'     se_art: standard error of HIV+ART+ exposure coefficient

run_analysis <- function(dataset_cp, method, include_no_testers) {
  
  # !!!!! `method` argument currently ignored for "ideal" and "mi"
  
  # Create "censor" dataset
  if (method=="censor") {
    
    # Add an ID row
    dataset_cp <- cbind("obs_id"=c(1:nrow(dataset_cp)),dataset_cp)
    
    # Exclude patients with no testing data
    dataset_cp %<>% filter(case!=1)
    
    # Exclude observation time prior to the first test
    dataset_cp %<>% filter(start_year>=first_test)
    
    # Exclude observation time after the last negative test for case 2
    dataset_cp %<>% filter(
      !(replace_na(case==2 & start_year>last_neg_test,FALSE))
    )
    
    # !!!!! Check censoring manually
    
  }
  
  if (method %in% c("ideal", "mi")) {
    
    # !!!!! Archive this functionality; it creates severe selection bias
    
    if (include_no_testers==FALSE) {
      dataset_cp %<>% filter(case!=1)
    }
    
  }
  
  # Fit time-varying Cox model ("ideal")
  fit <- coxph(
    Surv(start_year, end_year, died) ~ factor(casc_status) + age_bin +
      cluster(patient_id),
    data = dataset_cp
  )
  summ <- summary(fit)$coefficients
  
  results <- list(
    est_hiv = summ["factor(casc_status)HIV+ART-","coef"],
    se_hiv = summ["factor(casc_status)HIV+ART-","se(coef)"],
    est_art = summ["factor(casc_status)HIV+ART+","coef"],
    se_art = summ["factor(casc_status)HIV+ART+","se(coef)"]
  )
  
  return(results)
  
}
