#' Run a single simulation replicate
#'
#' @return A list containing the following:
#'     est_complete: HR estimate from complete data
#'     est_missing: HR estimate from MI procedure
#'     se_complete: SE of HR estimate from complete data
#'     se_missing: SE of HR estimate from MI procedure

one_simulation <- function() {
  
  dat <- generate_data(
    n = L$n,
    max_time = L$max_time,
    params = L$params
  )
  
  # Log-likelihood for a single individual
  #' @param x_j X_j Serostatus at time j
  #' @param x_j1 X_{j-1} serostatus at time j-1
  #' @param y_j Y_j outcome indicator at time j
  #' @param z_j Z_j covariate vector at time j
  #' @param prm Parameter vector: (alpha_j, gamma, beta_j, xi, beta_x)
  loglik_ij <- function(x_j, x_j1, y_j, z_j, prm) {
    
    if (x_j1==1) {
      piece_x <- 1
    } else {
      # !!!!! Add a "t" argument for linear baseline hazard
      exp_lin_x <- exp(prm[1]+prm[2]*z_j[1]+prm[3]*z_j[2])
      piece_x <- ifelse(x_j==1, exp_lin_x, 1-exp_lin_x)
    }
    
    exp_lin_y <- exp(prm[4]+prm[5]*z_j[1]+prm[6]*z_j[2]+prm[7]*x_j)
    piece_y <- ifelse(y_j==1, exp_lin_y, 1-exp_lin_y)
    
    return(log(max(piece_x,1e-8))+log(max(piece_y,1e-8)))
    
  }
  
  # Log-likelihood across individuals and time
  negloglik <- function(prm) {
    
    -1 * sum(apply(
      X = dat,
      MARGIN = 1,
      FUN = function(r) {
        loglik_ij(
          x_j = r[["x"]],
          x_j1 = r[["x_prev"]],
          y_j = r[["y"]],
          z_j = c(r[["z_sex"]], r[["z_age"]]),
          prm = prm
        )
      }
    ))
    
  }
  
  # Run optimizer
  opt <- solnp(
    pars = log(c(0.002, 1, 1, 0.002, 1, 1, 1)),
    fun = negloglik
  )
  if (opt$convergence!=0) { warning("solnp() did not converge") }
  
  # Extract estimates and compute SEs
  lik_ests <- exp(opt$pars)
  lik_se <- sqrt(diag(solve(opt$hessian)))
  lik_ci_lo <- exp(opt$pars-1.96*lik_se)
  lik_ci_hi <- exp(opt$pars+1.96*lik_se)
  
  # Get estimates and SEs from Cox model
  model <- coxph(
    Surv(t_start, t_end, y) ~ z_sex + z_age + x + cluster(id),
    data = dat
  )
  cox_ests <- as.numeric(summary(model)$conf.int[,1])
  cox_ci_lo <- as.numeric(summary(model)$conf.int[,3])
  cox_ci_hi <- as.numeric(summary(model)$conf.int[,4])
  
  return(list(
    lik_beta_x_est = lik_ests[7],
    lik_beta_x_ci_lo = lik_ci_lo[7],
    lik_beta_x_ci_hi = lik_ci_hi[7],
    lik_gamma_1_est = lik_ests[5],
    lik_gamma_1_lo = lik_ci_lo[5],
    lik_gamma_1_hi = lik_ci_hi[5],
    lik_gamma_2_est = lik_ests[6],
    lik_gamma_2_lo = lik_ci_lo[6],
    lik_gamma_2_hi = lik_ci_hi[6],
    cox_beta_x_est = cox_ests[3],
    cox_beta_x_ci_lo = cox_ci_lo[3],
    cox_beta_x_ci_hi = cox_ci_hi[3],
    cox_gamma_1_est = cox_ests[1],
    cox_gamma_1_lo = cox_ci_lo[1],
    cox_gamma_1_hi = cox_ci_hi[1],
    cox_gamma_2_est = cox_ests[2],
    cox_gamma_2_lo = cox_ci_lo[2],
    cox_gamma_2_hi = cox_ci_hi[2]
  ))
  
  
  
  # !!!!! CONTINUE !!!!!
  
  # dat2 <- impose_missingness(dat) # !!!!!
  # 
  # 
  # # !!!!! Testing
  # # C <- list(num_patients=10, start_year=2000, end_year=2001, m=5)
  # # C <- list(num_patients=5000, start_year=2000, end_year=2002, m=5)
  # # L <- list(method="mi", hr_hiv=1.4, hr_art=0.7)
  # 
  # # Generate baseline data
  # # !!!!! This is generated as a "true cohort" rather than an "open cohort"
  # dat_baseline <- generate_data_baseline(
  #   num_patients = C$num_patients,
  #   start_year = C$start_year
  # )
  # 
  # # Set parameters
  # params <- list(
  #   alpha0=-4,  alpha1=0.1,  alpha2=0.05,  alpha3=0, # alpha3=0.2
  #   beta0=-4,   beta1=0.1,   beta2=0.05,   beta3=0,  # beta3=0.2
  #   eta0=-4,    eta1=0.1,    eta2=0.05,    eta3=0,   # eta3=0.2
  #   gamma0=-4,  gamma1=0.1,  gamma2=0.05,  gamma3=0, # gamma3=0.2
  #   psi1=L$hr_hiv,
  #   psi2=L$hr_art
  # )
  # 
  # # Generate event data
  # # !!!!! For now, all patients are HIV- at baseline
  # dat_events <- apply(dat_baseline, MARGIN=1, function(r) {
  #   generate_data_events(
  #     id = r[["id"]],
  #     b_age = r[["b_age"]],
  #     sex = r[["sex"]],
  #     u = r[["u"]],
  #     start_year = C$start_year,
  #     end_year = C$end_year,
  #     baseline_status = NA,
  #     params = params
  #   )
  # })
  # attr(dat_events, "end_year") <- C$end_year
  
  # # Take m samples from the posterior
  # # theta_m <- posterior_param_sample(fit=fit, size=C$m)
  # # !!!!! Temp: START
  # {
  #   psi1_psample <- rnorm(C$m, mean=L$hr_hiv, sd=0.1)
  #   psi2_psample <- rnorm(C$m, mean=L$hr_art, sd=0.1)
  #   theta_m <- list()
  #   for (i in 1:C$m) {
  #     theta_m[[i]] <- params
  #     theta_m[[i]]$psi1 <- psi1_psample[i]
  #     theta_m[[i]]$psi2 <- psi2_psample[i]
  #   }
  # }
  # # !!!!! Temp: END
  
  # if (L$method=="ideal") {
  #   # Transform data and run Cox PH analysis
  #   dat_cp <- transform_dataset(
  #     dat_baseline = dat_baseline,
  #     dat_events = dat_events
  #   )
  #   results <- run_analysis(dat_cp=dat_cp)
  #   # print(exp(results$est_hiv)); print(exp(results$est_art)); # !!!!!
  # }
  
  # if (L$method=="mi") {
  #   
  #   # Perform MI on second dataset
  #   dat_events_mi <- list()
  #   for (i in 1:C$m) {
  #     dat_events_mi[[i]] <- perform_imputation(
  #       dat_baseline = dat_baseline,
  #       dat_events = dat_events,
  #       theta_m = theta_m[[i]]
  #     )
  #   }
  #   
  #   results_mi <- list()
  #   for (i in 1:C$m) {
  #     # Transform data and run Cox PH analysis
  #     dat_cp <- transform_dataset(
  #       dat_baseline = dat_baseline,
  #       dat_events = dat_events_mi[[i]]
  #     )
  #     results_mi[[i]] <- run_analysis(dat_cp=dat_cp)
  #     # print(paste("Replicate:",i)) # !!!!!
  #     # print(paste("HIV:", exp(results_mi[[i]]$est_hiv))) # !!!!!
  #     # print(paste("ART:", exp(results_mi[[i]]$est_art))) # !!!!!
  #   }
  #   
  #   # Combine MI estimates using "Rubin's rules"
  #   v <- function(results_mi, attr) {
  #     sapply(results_mi, function(r) { r[[attr]] })
  #   }
  #   est_hiv <- v(results_mi, "est_hiv")
  #   est_art <- v(results_mi, "est_art")
  #   var_hiv <- (v(results_mi, "se_hiv"))^2
  #   var_art <- (v(results_mi, "se_art"))^2
  #   results <- list(
  #     est_hiv = mean(est_hiv),
  #     se_hiv = sqrt( mean(var_hiv) + (1+(1/C$m))*var(est_hiv) ),
  #     est_art = mean(est_art),
  #     se_art = sqrt( mean(var_art) + (1+(1/C$m))*var(est_art) )
  #   )
  #   
  # }
  
  # return (res)
  
}
