#' Run a single simulation replicate
#'
#' @return A list containing the following:
#'     est_complete: HR estimate from complete data
#'     est_missing: HR estimate from MI procedure
#'     se_complete: SE of HR estimate from complete data
#'     se_missing: SE of HR estimate from MI procedure

one_simulation <- function() {
  
  # Generate data
  dat <- generate_data(
    n = 50,
    # n = 500,
    max_time = 10,
    # max_time = 100,
    params = list(a = c(log(10),log(1.05)),
                  b = c(log(1.2),log(1.001),log(1.5)))
    # params = list(a = c(log(1.3),log(1.002)),
    #               b = c(log(1.2),log(1.001),log(1.5)))
  )
  
  # Run Cox PH model
  model <- coxph(
    Surv(t_start, t_end, d) ~ z_sex + z_age + x + cluster(id),
    data = dat
  )
  summary(model)
  summ <- summary(model)$coefficients
  
  dat2 <- impose_missingness(dat) # !!!!!
  
  # Generate an incomplete dataset from a complete dataset
  # Probability of being assigned to a "case" is independent of all other vars
  impose_missingness <- function(dat) {
    
    # Generate the "testing case" variable
    #     Case 1: no testing data
    #     Case 2: most recent test was negative
    #     Case 3: negative test followed by a positive test
    #     Case 4: first test was positive
    n <- attr(dat, "n")
    # case <- sample(c(1:4), size=n, replace=T)
    dat$case <- rep(NA, nrow(dat))
    dat$test_reg <- rep(NA, nrow(dat))
    
    # Loop through patients
    for (i in c(1:n)) {
      
      # Extract variables from dataset
      rows <- which(dat$id==i)
      x_i <- dat[rows,"x"]
      n_obs_i <- length(rows)
      max_time <- attr(dat, "max_time")
      
      # Assign to one of three testing regimens:
      #   Reg 1: 0 tests (prob 0.1)
      #   Reg 2: ceil(max_time/8) tests (prob 0.45)
      #   Reg 3: ceil(max_time/4) tests (prob 0.45)
      test_reg <- sample(c(1:3), size=1, prob=c(0.1,0.45,0.45))
      if (test_reg==2) {
        num_tests <- ceiling(max_time/8)
      } else if (test_reg==3) {
        num_tests <- ceiling(max_time/4)
      } else {
        num_tests <- 0
      }
      tests <- sort(sample(c(1:max_time), size=num_tests))
      tests <- tests[tests<=n_obs_i]
      if (length(tests)==0) {
        case <- 1 # No testing data
        x_i <- rep(NA, n_obs_i)
      } else {
        if (x_i[tests[1]]==1) {
          case <- 4 # First test POS
          x_i[c(1:round(tests[1]-1))] <- NA
        } else if (x_i[tests[length(tests)]]==0) {
          case <- 2 # most recent test NEG
          if (n_obs_i>tests[length(tests)]) {
            x_i[c(1:n_obs_i)>tests[length(tests)]] <- NA
          }
        } else {
          case <- 3 # NEG test then POS test
          ind_neg <- max(tests[x_i[tests]==0])
          ind_pos <- min(tests[x_i[tests]==1])
          if (ind_pos-ind_neg>1) {
            x_i[c(round(ind_neg+1):round(ind_pos-1))] <- NA
          }
        }
      }
      
      dat[rows,"x"] <- x_i
      dat[rows,"test_reg"] <- test_reg
      dat[rows,"case"] <- case
      
      # # Set new X variable value based on case
      # if (case[i]==1) {
      #   x_new <- rep(NA, n_obs_i)
      # } else if (case[i]==2) {
      #   
      #   x_new <- 999 # !!!!!
      # } else if (case[i]==3) {
      #   x_new <- 999 # !!!!!
      # } else if (case[i]==4) {
      #   x_new <- 999 # !!!!!
      # }
      
      # # Add case variable to dataframe
      # dat[rows,"case"] <- case[i]
      
    }
    
    return (dat)
    
  }
  
  # !!!!! Testing
  # C <- list(num_patients=10, start_year=2000, end_year=2001, m=5)
  # C <- list(num_patients=5000, start_year=2000, end_year=2002, m=5)
  # L <- list(method="mi", hr_hiv=1.4, hr_art=0.7)
  
  # Generate baseline data
  # !!!!! This is generated as a "true cohort" rather than an "open cohort"
  dat_baseline <- generate_data_baseline(
    num_patients = C$num_patients,
    start_year = C$start_year
  )
  
  # Set parameters
  params <- list(
    alpha0=-4,  alpha1=0.1,  alpha2=0.05,  alpha3=0, # alpha3=0.2
    beta0=-4,   beta1=0.1,   beta2=0.05,   beta3=0,  # beta3=0.2
    eta0=-4,    eta1=0.1,    eta2=0.05,    eta3=0,   # eta3=0.2
    gamma0=-4,  gamma1=0.1,  gamma2=0.05,  gamma3=0, # gamma3=0.2
    psi1=L$hr_hiv,
    psi2=L$hr_art
  )
  
  # Generate event data
  # !!!!! For now, all patients are HIV- at baseline
  dat_events <- apply(dat_baseline, MARGIN=1, function(r) {
    generate_data_events(
      id = r[["id"]],
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
  
  # # Transform data to JAGS format
  # dat_mcmc <- transform_mcmc(
  #   dat_baseline = dat_baseline,
  #   dat_events = dat_events
  # )
  # 
  # # Set MCMC params
  # mcmc <- list(n.adapt=1000, n.burn=1000, n.iter=1000, thin=1, n.chains=2)
  # 
  # # Fit the model in Stan
  # fit <- fit_stan(
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
      theta_m[[i]] <- params
      theta_m[[i]]$psi1 <- psi1_psample[i]
      theta_m[[i]]$psi2 <- psi2_psample[i]
    }
  }
  # !!!!! Temp: END
  
  if (L$method=="ideal") {
    # Transform data and run Cox PH analysis
    dat_cp <- transform_dataset(
      dat_baseline = dat_baseline,
      dat_events = dat_events
    )
    results <- run_analysis(dat_cp=dat_cp)
    # print(exp(results$est_hiv)); print(exp(results$est_art)); # !!!!!
  }
  
  if (L$method=="mi") {
    
    # Perform MI on second dataset
    dat_events_mi <- list()
    for (i in 1:C$m) {
      dat_events_mi[[i]] <- perform_imputation(
        dat_baseline = dat_baseline,
        dat_events = dat_events,
        theta_m = theta_m[[i]]
      )
    }
    
    results_mi <- list()
    for (i in 1:C$m) {
      # Transform data and run Cox PH analysis
      dat_cp <- transform_dataset(
        dat_baseline = dat_baseline,
        dat_events = dat_events_mi[[i]]
      )
      results_mi[[i]] <- run_analysis(dat_cp=dat_cp)
      # print(paste("Replicate:",i)) # !!!!!
      # print(paste("HIV:", exp(results_mi[[i]]$est_hiv))) # !!!!!
      # print(paste("ART:", exp(results_mi[[i]]$est_art))) # !!!!!
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
  
  if (L$method=="censor") {
    # !!!!! Ask Mark what he does currently
  }
  
  # Add CI bounds
  results$hiv_ci_l <- results$est_hiv - 1.96*results$se_hiv
  results$hiv_ci_u <- results$est_hiv + 1.96*results$se_hiv
  results$art_ci_l <- results$est_art - 1.96*results$se_art
  results$art_ci_u <- results$est_art + 1.96*results$se_art
  
  return (results)
  
}
