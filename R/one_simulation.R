#' Run a single simulation replicate
#'
#' @return A list containing the following:
#'     est_complete: HR estimate from complete data
#'     est_missing: HR estimate from MI procedure
#'     se_complete: SE of HR estimate from complete data
#'     se_missing: SE of HR estimate from MI procedure

one_simulation <- function() {
  
  # Generate baseline data
  dat_baseline <- generate_data_baseline(
    num_patients = C$num_patients,
    start_year = C$start_year
  )
  
  # # !!!!! Testing
  # dat_baseline <- data.frame(
  #   id = 1:20,
  #   b_age = sample(1:80, size=20, replace=TRUE),
  #   sex = sample(c(0,1), size=20, replace=TRUE),
  #   u = rnorm(n=20)
  # )
  
  # Generate event data
  params <- list(
    alpha0=-5,  alpha1=0.1,  alpha2=0.05,  alpha3=0.2,
    beta0=-5,   beta1=0.1,   beta2=0.05,   beta3=0.2,
    eta0=-5,    eta1=0.1,    eta2=0.05,    eta3=0.2,
    gamma0=-5,  gamma1=0.1,  gamma2=0.05,  gamma3=0.2,
    psi1=L$hr_hiv,
    psi2=L$hr_art
  )
  dat_events <- apply(
    X = dat_baseline,
    MARGIN = 1,
    FUN = function(r) {
      generate_data_events(
        b_age = r[["b_age"]],
        sex = r[["sex"]],
        u = r[["u"]],
        start_year = C$start_year,
        end_year = C$end_year,
        # baseline_status, # !!!!! Right now everyone starts as HIV-
        params = params
      )
    }
  )
  
  # !!!!! Continue
  
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
