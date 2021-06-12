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
  L <- list(method="mi", hr_hiv=1.3, hr_art=0.6)
  
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
  # !!!!! This is generated as a "true cohort" rather than an "open cohort"
  dat_baseline <- generate_data_baseline(
    num_patients = C$num_patients,
    start_year = C$start_year
  )
  
  # Generate event data
  # !!!!! For now, all patients are HIV- at baseline
  dat_events <- apply(dat_baseline, MARGIN=1, function(r) {
    generate_data_events(
      b_age = r[["b_age"]],
      sex = r[["sex"]],
      u = r[["u"]],
      start_year = C$start_year,
      end_year = C$end_year,
      baseline_status = NA,
      params = params
    )
  })
  attr(dat_events, "end_year") <- C$end_year
  
  # Transform data to JAGS format
  dat_mcmc <- transform_mcmc(
    dat_baseline = dat_baseline,
    dat_events = dat_events
  )
  
  # Set MCMC params
  mcmc <- list(n.adapt=1000, n.burn=1000, n.iter=1000, thin=1, n.chains=2)
  
  # # Fit the model in JAGS
  # # !!!!! Migrating to Stan
  # fit <- fit_jags(
  #   dat = dat_mcmc,
  #   mcmc = mcmc
  # )
  
  # Take m samples from the posterior
  # theta_m <- posterior_param_sample(fit=fit, size=C$m)
  # !!!!! Temp: START
  {
    psi1_psample <- rnorm(C$m, mean=L$hr_hiv, sd=0.1)
    psi2_psample <- rnorm(C$m, mean=L$hr_art, sd=0.1)
    theta_m <- list()
    for (i in 1:C$m) {
      theta_m[[i]] <- list(
        psi1 = psi1_psample[i],
        psi2 = psi2_psample[i]
      )
    }
  }
  # !!!!! Temp: END
  
  if (L$method=="ideal") {
    # Transform data and run Cox PH analysis
    dat_cp <- transform_dataset(
      dat_baseline = dat_baseline,
      dat_events = dat_events
    )
    results <- run_analysis(
      dat_cp = dat_cp,
      method = "ideal"
    )
  }
  
  if (L$method=="censor") {
    # !!!!! Ask Mark what he does currently
  }
  
  if (L$method=="mi") {
    
    # Perform MI on second dataset
    datasets_mi <- list()
    for (i in 1:C$m) {
      datasets_mi[[i]] <- perform_imputation(
        # !!!!! Make sure we are memoising within perform_imputation
        dataset,
        C$p_sero_year,
        theta_m = theta_m[[i]]
      )
    }
    
    # Transform data and run Cox PH analysis on imputed datasets
    results_mi <- list()
    for (i in 1:C$m) {
      dat_cp <- transform_dataset(
        datasets_mi[[i]] # !!!!! not in correct format
      )
      results_mi[[i]] <- run_analysis(
        dat_cp = dat_cp,
        method = "mi"
      )
    }
    
    # Combine MI estimates using "Rubin's rules"
    v <- function(results_mi, attr) {
      sapply(results_mi, function(r) { r[[attr]] })
    }
    est_hiv <- v(results_mi, "est_hiv")
    est_art <- v(results_mi, "est_art")
    var_hiv <- (v(results_mi, "se_hiv"))^2
    var_art <- (v(results_mi, "se_art"))^2
    results <- list(
      est_hiv = mean(est_hiv),
      se_hiv = sqrt( mean(var_hiv) + (1+(1/C$m))*var(est_hiv) ),
      est_art = mean(est_art),
      se_art = sqrt( mean(var_art) + (1+(1/C$m))*var(est_art) )
    )
    
  }
  
  return (results)
  
}
