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
  
  # # Run optimizer
  # opt <- solnp(
  #   pars = log(c(0.002, 1, 1, 0.002, 1, 1, 1)),
  #   fun = negloglik
  # )
  # if (opt$convergence!=0) { warning("solnp() did not converge") }
  
  func <- function(par) { negloglik(par) }
  opt <- optim(par=log(c(0.002, 1, 1, 0.002, 1, 1, 1)), fn=func)
  h <- optimHess(par=opt$par, fn=func)
  lik_full <- list(ests=exp(opt$par))
  lik_full$se <- sqrt(diag(solve(h)))
  lik_full$ci_lo <- exp(lik_full$ests-1.96*lik_full$se)
  lik_full$ci_hi <- exp(lik_full$ests+1.96*lik_full$se)
  
  # Get estimates and SEs from Cox model
  model <- coxph(
    Surv(t_start, t_end, y) ~ w_sex + w_age + x + cluster(id),
    data = dat
  )
  cox_full <- list(
    ests = as.numeric(summary(model)$conf.int[,1]),
    ci_lo = as.numeric(summary(model)$conf.int[,3]),
    ci_hi = as.numeric(summary(model)$conf.int[,4])
  )
  
  return(list(
    lik_F_beta_est = lik_full$ests[7],
    lik_F_beta_ci_lo = lik_full$ci_lo[7],
    lik_F_beta_ci_hi = lik_full$ci_hi[7],
    lik_F_g_x_1_est = lik_full$ests[5],
    lik_F_g_x_1_lo = lik_full$ci_lo[5],
    lik_F_g_x_1_hi = lik_full$ci_hi[5],
    lik_F_g_x_2_est = lik_full$ests[6],
    lik_F_g_x_2_lo = lik_full$ci_lo[6],
    lik_F_g_x_2_hi = lik_full$ci_hi[6],
    cox_beta_est = cox_full$ests[3],
    cox_beta_ci_lo = cox_full$ci_lo[3],
    cox_beta_ci_hi = cox_full$ci_hi[3],
    cox_g_x_1_est = cox_full$ests[1],
    cox_g_x_1_lo = cox_full$ci_lo[1],
    cox_g_x_1_hi = cox_full$ci_hi[1],
    cox_g_x_2_est = cox_full$ests[2],
    cox_g_x_2_lo = cox_full$ci_lo[2],
    cox_g_x_2_hi = cox_full$ci_hi[2]
  ))
  
}
