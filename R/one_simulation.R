#' Run a single simulation replicate
#'
#' @return A list containing the following:
#'     est_complete: HR estimate from complete data
#'     est_missing: HR estimate from MI procedure
#'     se_complete: SE of HR estimate from complete data
#'     se_missing: SE of HR estimate from MI procedure

one_simulation <- function() {
  
  chk(0, "START")
  
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
  # True values: log(c(0.005,1.3,1.2,0.003,1.2,1.1,1.5,0.6,1,1))
  par <- log(c(0.003,1.2,1.1,0.002,1.3,1,1.3,0.8,0.999,1.001)) # Starting values
  names(par) <- c("a_x", "g_x1", "g_x2", "a_y", "g_y1", "g_y2", "beta_x",
                  "beta_z", "t_x", "t_y")
  
  chk(2, "construct_negloglik_miss: START")
  negloglik_miss <- construct_negloglik_miss(dat)
  chk(2, "construct_negloglik_miss: END")
  chk(2, "optim: START")
  opt_miss <- optim(par=par, fn=negloglik_miss)
  browser() # !!!!!
  chk(2, "optim: END")
  chk(2, "hessian: START")
  hessian_miss <- hessian(func=negloglik_miss, x=opt_miss$par)
  chk(2, "hessian: END")
  lik_miss <- list(ests=opt_miss$par)
  lik_miss$se <- sqrt(diag(solve(hessian_miss)))
  lik_miss$ci_lo <- lik_miss$ests-1.96*lik_miss$se
  lik_miss$ci_hi <- lik_miss$ests+1.96*lik_miss$se

  res <- list() # !!!!! 2022-09-28
  for (i in c(1:length(par))) { # !!!!! 2022-09-28
    res[[paste0("lik_M_",names(par)[i],"_est")]] <- as.numeric(opt_miss$par[i]) # !!!!! 2022-09-28
    res[[paste0("lik_M_",names(par)[i],"_se")]] <- sqrt(diag(solve(hessian_miss)))[i] # !!!!! 2022-09-28
  } # !!!!! 2022-09-28
  
  # Old code; archive???
  if (F) {
    res <- list(
      lik_M_beta_x_est = lik_miss$ests[7],
      lik_M_beta_x_ci_lo = lik_miss$ci_lo[7],
      lik_M_beta_x_ci_hi = lik_miss$ci_hi[7],
      lik_M_g_y1_est = lik_miss$ests[5],
      lik_M_g_y1_lo = lik_miss$ci_lo[5],
      lik_M_g_y1_hi = lik_miss$ci_hi[5],
      lik_M_g_y2_est = lik_miss$ests[6],
      lik_M_g_y2_lo = lik_miss$ci_lo[6],
      lik_M_g_y2_hi = lik_miss$ci_hi[6]
    )
  }
  
  chk(3, "END")
  
  return(res)
  
}
