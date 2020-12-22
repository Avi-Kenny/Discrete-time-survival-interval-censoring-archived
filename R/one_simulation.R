#' Run a single simulation replicate
#'
#' @return A list containing the following:
#'     est_complete: HR estimate from complete data
#'     est_missing: HR estimate from MI procedure
#'     se_complete: SE of HR estimate from complete data
#'     se_missing: SE of HR estimate from MI procedure

one_simulation <- function() {
  
  # Create dataset with known cascade status dates
  dataset <- generate_dataset(
    num_patients = C$num_patients,
    start_year = C$start_year,
    end_year = C$end_year,
    hazard_ratios = list("hiv"=L$hr_hiv, "art"=L$hr_art),
    p_sero_year = C$p_sero_year,
    p_death_year = p_death_year(L$p_death_mult),
    u_mult = list("sero"=L$u_mult_sero, "death"=L$u_mult_death)
  )
  
  if (L$method=="ideal") {
    # Transform data and run Cox PH analysis
    dataset_cp <- transform_dataset(dataset)
    results <- run_analysis(
      dataset_cp = dataset_cp,
      method = "ideal"
    )
  }
  
  if (L$method=="censor") {
    # Remove unknown seroconversion info
    # dataset %<>% mutate(sero_year=NA)
    
    # !!!!! Need to remove sero info and impute; ask Mark what he does currently
    
    # Transform data and run Cox PH analysis
    dataset_cp <- transform_dataset(dataset)
    results <- run_analysis(
      dataset_cp = dataset_cp,
      method = "censor"
    )
  }
  
  if (L$method=="mi") {
    
    # Perform MI on second dataset
    datasets_mi <- list()
    for (i in 1:C$m) {
      datasets_mi[[i]] <- perform_imputation(dataset, C$p_sero_year)
    }
    
    # Transform data and run Cox PH analysis on imputed datasets
    results_mi <- list()
    for (i in 1:C$m) {
      results_mi[[i]] <- run_analysis(
        dataset_cp = transform_dataset(datasets_mi[[i]]),
        method = "mi"
      )
    }
    
    # Combine MI estimates
    est_hiv <- sapply(results_mi, function(res) {res$est_hiv})
    var_hiv <- sapply(results_mi, function(res) {res$se_hiv^2})
    est_art <- sapply(results_mi, function(res) {res$est_art})
    var_art <- sapply(results_mi, function(res) {res$se_art^2})
    results <- list(
      est_hiv = mean(est_hiv),
      se_hiv = sqrt( mean(var_hiv) + (1+(1/C$m))*var(est_hiv) ),
      est_art = mean(est_art),
      se_art = sqrt( mean(var_art) + (1+(1/C$m))*var(est_art) )
    )
    
  }
  
  return (results)
  
}
