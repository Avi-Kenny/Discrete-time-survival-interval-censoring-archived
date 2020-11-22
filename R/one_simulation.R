#' Run a single simulation replicate
#'
#' @param x TO DO
#' @return A list containing the following:
#'     x: descr
#'     y: descr

one_simulation <- function(rel_risk) {
  
  # Set local variables
  # num_patients <- 1000
  # end_date <- (2019-1900)*12
  m <- 5 # Number of MI replicates
  
  # Create dataset with known cascade status dates
  df_complete <- generate_dataset(
    rel_risk = rel_risk,
    num_patients, end_date,
    rel_risk, missingness
  )
  
  # Create second dataset by imposing missingness structure
  df_missing <- impose_missingness(df_complete)
  
  # Perform MI on second dataset
  dfs_mi <- list()
  for (i in 1:m) {
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
