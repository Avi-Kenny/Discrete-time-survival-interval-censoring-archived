#' Run Cox PH analysis
#'
#' @param dat_cp A dataset returned by transform_dataset()
#' @param options Placeholder; currently unused
#' @return A list containing the following:
#'     est_hiv: point estimate of HIV+ART- exposure coefficient
#'     se_hiv: standard error of HIV+ART- exposure coefficient
#'     est_art: point estimate of HIV+ART+ exposure coefficient
#'     se_art: standard error of HIV+ART+ exposure coefficient

run_analysis <- function(dat_cp, options=list()) {
  
  # Create "censor" dataset
  if (!is.null(options$method) && options$method=="censor") {
    
    # Add an ID row
    dat_cp <- cbind("obs_id"=c(1:nrow(dat_cp)),dat_cp)
    
    # Exclude patients with no testing data
    dat_cp %<>% filter(case!=1)
    
    # Exclude observation time prior to the first test
    dat_cp %<>% filter(start_year>=first_test)
    
    # Exclude observation time after the last negative test for case 2
    dat_cp %<>% filter(
      !(replace_na(case==2 & start_year>last_neg_test,FALSE))
    )
    
    # !!!!! Check censoring manually
    
  }
  
  # Fit time-varying Cox model ("ideal")
  fit <- coxph(
    Surv(start_time, end_time, y) ~ factor(casc_status) + age + cluster(id),
    data = dat_cp
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
