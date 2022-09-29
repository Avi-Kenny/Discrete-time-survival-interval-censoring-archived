#' Run a single simulation replicate
#'
#' @return A list containing the following:
#'     est_complete: HR estimate from complete data
#'     est_missing: HR estimate from MI procedure
#'     se_complete: SE of HR estimate from complete data
#'     se_missing: SE of HR estimate from MI procedure

one_simulation <- function() {
  
  chk(0, "Start")
  
  # Generate dataset
  dat <- generate_data(
    n = L$n,
    max_time = L$max_time,
    params = L$params
  )
  chk(1, "Data generated")
  
  # Add x_prev column
  dat %<>% arrange(id, t_end)
  dat$x_prev <- ifelse(
    dat$id==c(0,dat$id[c(1:(length(dat$id)-1))]),
    c(0,dat$x[c(1:(length(dat$x)-1))]),
    0
  )
  
  if (F) {
    
    # Run optimizer
    opt <- solnp(
      pars = log(c(0.002, 1, 1, 0.002, 1, 1, 1)),
      fun = negloglik_full
    )
    if (opt$convergence!=0) { warning("solnp() did not converge") }
    
  } # Rsolnp optimizer
  
  if (F) {
    
    par <- log(c(0.003,1.2,1.1,0.002,1.1,1,1.3,0.4,1.1,1.2)) # Starting values
    # par <- log(c(0.005,1.3,1.2,0.003,1.2,1.1,1.5,0.4,1.2,1.1)) # True values (except baseline testing rate)
    
    # Run optimizer (full data structure)
    chk(2, "negloglik_full: START")
    fnc_full <- function(par) { negloglik_full(dat, par) }
    opt_full <- optim(par=par[1:7], fn=fnc_full)
    chk(2, "negloglik_full: optimizer done")
    # hessian_full <- optimHess(par=opt_full$par, fn=fnc_full) # !!!!!
    lik_full <- list(ests=opt_full$par)
    chk(2, "negloglik_full: hessian done")
    # lik_full <- list(ests=exp(opt_full$par))
    # lik_full$se <- sqrt(diag(solve(hessian_full))) # !!!!!
    # lik_full$ci_lo <- exp(lik_full$ests-1.96*lik_full$se) # !!!!!
    # lik_full$ci_hi <- exp(lik_full$ests+1.96*lik_full$se) # !!!!!
    
  } # DEBUG: full data structure
  
  if (F) {
    
    par <- log(c(0.003,0.002,1.3,0.4)) # Starting values (a_x, a_y, beta, a_v)
    fnc_miss <- function(par) { negloglik_miss_nocovariates(dat, par) }
    opt_miss <- optim(par=par, fn=fnc_miss) # method="BFGS"
    lik_miss <- list(ests=opt_miss$par)
    # lik_miss <- list(ests=exp(opt_miss$par))
    res <- list(
      lik_M_a_x_est = lik_miss$ests[1],
      lik_M_a_y_est = lik_miss$ests[2],
      lik_M_beta_est = lik_miss$ests[3],
      lik_M_a_v_est = lik_miss$ests[4]
    )
    
  } # DEBUG: missing data structure, no covariates
  
  if (F) {
    
    # Get estimates and SEs from Cox model
    chk(10, "Cox: START")
    model <- coxph(
      Surv(t_start, t_end, y) ~ w_sex + w_age + x + cluster(id),
      data = dat
    )
    cox_full <- list(
      ests = as.numeric(summary(model)$conf.int[,1]),
      ci_lo = as.numeric(summary(model)$conf.int[,3]),
      ci_hi = as.numeric(summary(model)$conf.int[,4])
    )
    chk(11, "Cox: END")
    
  } # DEBUG: Cox model (ideal data structure)
  
  # Set initial parameters
  # True values: log(c(0.005,1.3,1.2,0.003,1.2,1.1,1.5))
  # par <- log(c(0.003,1.2,1.1,0.002,1.1,1,1.3)) # Starting values
  par <- log(c(0.003,1.2,1.1,0.002,1.3,1,1.3)) # Starting values # !!!!! 2022-09-28
  names(par) <- c("a_x", "g_x1", "g_x2", "a_y", "g_y1", "g_y2", "beta")
  
  # pars_to_fix <- c(1,2,3,4,6) # !!!!! 2022-09-28
  # par <- par[-pars_to_fix] # !!!!! 2022-09-28
  
  # system.time({
  chk(2, "negloglik_miss: START")
  fnc_miss <- function(par) { negloglik_miss(dat, par) }
  opt_miss <- optim(par=par, fn=fnc_miss)
  chk(2, "negloglik_miss: optimizer done")
  # hessian_miss <- optimHess(par=opt_miss$par, fn=fnc_miss)
  hessian_miss <- hessian(func=fnc_miss, x=opt_miss$par)
  chk(2, "negloglik_miss: hessian done")
  lik_miss <- list(ests=opt_miss$par)
  lik_miss$se <- sqrt(diag(solve(hessian_miss)))
  chk(2, "negloglik_miss: SE extraction done")
  lik_miss$ci_lo <- lik_miss$ests-1.96*lik_miss$se
  lik_miss$ci_hi <- lik_miss$ests+1.96*lik_miss$se
  chk(2, "negloglik_miss: END")
  # })
  
  res <- list() # !!!!! 2022-09-28
  for (i in c(1:length(par))) { # !!!!! 2022-09-28
    res[[paste0("lik_M_",names(par)[i],"_est")]] <- as.numeric(opt_miss$par[i]) # !!!!! 2022-09-28
    res[[paste0("lik_M_",names(par)[i],"_se")]] <- sqrt(diag(solve(hessian_miss)))[i] # !!!!! 2022-09-28
  } # !!!!! 2022-09-28
  
  # res <- list(
  #   lik_M_beta_est = lik_miss$ests[7],
  #   lik_M_beta_ci_lo = lik_miss$ci_lo[7],
  #   lik_M_beta_ci_hi = lik_miss$ci_hi[7],
  #   lik_M_g_y1_est = lik_miss$ests[5],
  #   lik_M_g_y1_lo = lik_miss$ci_lo[5],
  #   lik_M_g_y1_hi = lik_miss$ci_hi[5],
  #   lik_M_g_y2_est = lik_miss$ests[6],
  #   lik_M_g_y2_lo = lik_miss$ci_lo[6],
  #   lik_M_g_y2_hi = lik_miss$ci_hi[6]
  # )
  
  return(res)
  
}
