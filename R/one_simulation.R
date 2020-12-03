#' Run a single simulation replicate
#'
#' @return A list containing the following:
#'     est_complete: HR estimate from complete data
#'     est_missing: HR estimate from MI procedure
#'     se_complete: SE of HR estimate from complete data
#'     se_missing: SE of HR estimate from MI procedure

one_simulation <- function() {
  
  # Create dataset with known cascade status dates
  df_complete <- generate_dataset(
    num_patients = C$num_patients,
    start_year = C$start_year,
    end_year = C$end_year,
    hazard_ratios = list("hiv"=L$hr_hiv, "art"=L$hr_art),
    p_sero_year = C$p_sero_year,
    u_mult = list("sero"=L$u_mult_sero, "death"=L$u_mult_death)
  )
  
  # Create second dataset by imposing missingness structure
  df_missing <- impose_missingness(df_complete)
  
  # Perform MI on second dataset
  dfs_mi <- list()
  for (i in 1:(C$m)) {
    dfs_mi[[i]] <- perform_imputation(df_missing)
  }
  
  # Perform Cox PH regression on complete dataset
  results_complete <- run_analysis(df_complete)
  
  # Perform Cox PH regression on imputed datasets
  results_mi <- list()
  for (i in 1:m) {
    results_mi[[i]] <- run_analysis(df_complete)
  }
  
  # Combine MI estimates (using "Rubin's Rules")
  # !!!!! Check this
  # !!!!! May need to multiply var(estimates) by m/(m-1)
  estimates <- sapply(results_mi, function(x) {x$est})
  vars <- sapply(results_mi, function(x) {x$se^2})
  results_mi_combined <- list(
    est = mean(estimates),
    se = sqrt( var(estimates) + mean(vars) )
  )
  
  return (list(
    est_complete = results_complete$est,
    est_missing = results_mi_combined$est,
    se_complete = results_complete$se,
    se_missing = results_mi_combined$se
  ))
  
}
