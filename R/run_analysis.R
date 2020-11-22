#' Run Cox PH analysis
#'
#' @param dataset A dataset returned by either generate_dataset() or
#'     perform_imputation()
#' @return A list containing the following:
#'     est: point estimate of exposure coefficient
#'     se: standard error of exposure coefficient estimate

run_analysis <- function(dataset) {
  
  fit <- coxph(
    Surv(start_time, end_time, had_event)~factor(hiv_status),
    data = dataset
  )
  summ <- summary(fit)$coefficients
  results <- list(
    est = summ["exposure","coef"],
    se = summ["exposure","se(coef)"]
  )
  
  return(results)
  
}
