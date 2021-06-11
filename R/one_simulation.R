#' Run a single simulation replicate
#'
#' @return A list containing the following:
#'     est_complete: HR estimate from complete data
#'     est_missing: HR estimate from MI procedure
#'     se_complete: SE of HR estimate from complete data
#'     se_missing: SE of HR estimate from MI procedure

one_simulation <- function() {
  
  # !!!!! Testing
  # C <- list(num_patients=10, start_year=2000, end_year=2001, m=5)
  C <- list(num_patients=100, start_year=2000, end_year=2002, m=5)
  L <- list(hr_hiv=1.3, hr_art=0.6)
  
  # Set parameters
  params <- list(
    alpha0=-4,  alpha1=0.1,  alpha2=0.05,  alpha3=0.2,
    beta0=-4,   beta1=0.1,   beta2=0.05,   beta3=0.2,
    eta0=-4,    eta1=0.1,    eta2=0.05,    eta3=0.2,
    gamma0=-4,  gamma1=0.1,  gamma2=0.05,  gamma3=0.2,
    psi1=L$hr_hiv,
    psi2=L$hr_art
  )
  
  # Generate baseline data
  dat_baseline <- generate_data_baseline(
    num_patients = C$num_patients,
    start_year = C$start_year
  )
  
  # Generate event data
  # !!!!! For now, all patients are HIV- at baseline
  dat_events <- apply(X=dat_baseline, MARGIN=1, FUN=function(r) {
    generate_data_events(
      b_age = r[["b_age"]],
      sex = r[["sex"]],
      u = r[["u"]],
      start_year = C$start_year,
      end_year = C$end_year,
      params = params
    )
  })
  attr(dat_events, "end_year") <- C$end_year
  
  # Transform data to JAGS format
  dat_jags <- transform_jags(
    dat_baseline = dat_baseline,
    dat_events = dat_events
  )
  
  # Set MCMC params
  mcmc <- list(n.adapt=1000, n.burn=1000, n.iter=1000, thin=1, n.chains=2)
  
  # Fit the model in JAGS
  fit <- fit_jags(
    dat = dat_jags,
    mcmc = mcmc
  )
  
  # Take m samples from the posterior
  # theta_m <- posterior_param_sample(fit=fit, size=C$m)
  theta_m <- list(                             # !!!!!
    # alpha0 = alpha0,                         # !!!!!
    psi1 = rnorm(C$m, mean=L$hr_hiv, sd=0.1),  # !!!!!
    psi2 = rnorm(C$m, mean=L$hr_art, sd=0.1)   # !!!!!
  )                                            # !!!!!
  
  # !!!!! Continue
  
  if (L$method=="ideal") {
    # Transform data and run Cox PH analysis
    dat_cp <- transform_dataset(
      dat_baseline = dat_baseline,
      dat_events = dat_events
    )
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
