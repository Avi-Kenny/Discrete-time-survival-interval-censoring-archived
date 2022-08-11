#' Run a single simulation replicate
#'
#' @return A list containing the following:
#'     est_complete: HR estimate from complete data
#'     est_missing: HR estimate from MI procedure
#'     se_complete: SE of HR estimate from complete data
#'     se_missing: SE of HR estimate from MI procedure

one_simulation <- function() {
  
  # Generate dataset
  dat <- generate_data(
    n = L$n,
    max_time = L$max_time,
    params = L$params
  )
  
  # Add x_prev column
  dat %<>% arrange(id, t_end)
  dat$x_prev <- ifelse(
    dat$id==c(0,dat$id[c(1:(length(dat$id)-1))]),
    c(0,dat$x[c(1:(length(dat$x)-1))]),
    0
  )
  
  # Rsolnp optimizer
  if (F) {
    # # Run optimizer
    # opt <- solnp(
    #   pars = log(c(0.002, 1, 1, 0.002, 1, 1, 1)),
    #   fun = negloglik_full
    # )
    # if (opt$convergence!=0) { warning("solnp() did not converge") }
  }
  
  # !!!!! DEBUGGING !!!!!
  if (F) {
    
    par <- log(c(0.003,0.002,1.3,0.4)) # Starting values (a_x, a_y, beta, a_v)
    fnc_miss <- function(par) { negloglik_miss2(dat, par) }
    opt_miss <- optim(par=par, fn=fnc_miss) # method="BFGS"
    lik_miss <- list(ests=opt_miss$par)
    # lik_miss <- list(ests=exp(opt_miss$par))
    res <- list(
      lik_M_a_x_est = lik_miss$ests[1],
      lik_M_a_y_est = lik_miss$ests[2],
      lik_M_beta_est = lik_miss$ests[3],
      lik_M_a_v_est = lik_miss$ests[4]
    )
    
  } else {
    
    # Set initial parameters
    par <- log(c(0.003,1.2,1.001,0.002,1.1,1.0005,1.3,0.4,1.1,1.0005)) # Starting values
    # par <- log(c(0.005,1.3,1.002,0.003,1.2,1.001,1.5,0.4,1.2,1.001)) # True values (except baseline testing rate)
    
    # Run optimizer (full data structure)
    fnc_full <- function(par) { negloglik_full(dat, par) }
    opt_full <- optim(par=par[1:7], fn=fnc_full)
    # hessian_full <- optimHess(par=opt_full$par, fn=fnc_full) # !!!!!
    lik_full <- list(ests=opt_full$par)
    # lik_full <- list(ests=exp(opt_full$par))
    # lik_full$se <- sqrt(diag(solve(hessian_full))) # !!!!!
    # lik_full$ci_lo <- exp(lik_full$ests-1.96*lik_full$se) # !!!!!
    # lik_full$ci_hi <- exp(lik_full$ests+1.96*lik_full$se) # !!!!!
    
    # Run optimizer (missing data structure)
    fnc_miss <- function(par) { negloglik_miss(dat, par) }
    opt_miss <- optim(par=par, fn=fnc_miss)
    # hessian_miss <- optimHess(par=opt_miss$par, fn=fnc_miss) # !!!!!
    lik_miss <- list(ests=opt_miss$par)
    # lik_miss <- list(ests=exp(opt_miss$par))
    # lik_miss$se <- sqrt(diag(solve(hessian_miss))) # !!!!!
    # lik_miss$ci_lo <- exp(lik_miss$ests-1.96*lik_miss$se) # !!!!!
    # lik_miss$ci_hi <- exp(lik_miss$ests+1.96*lik_miss$se) # !!!!!
    
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
    
    res <- list(
      lik_M_beta_est = lik_miss$ests[7],
      # lik_M_beta_ci_lo = lik_miss$ci_lo[7],
      # lik_M_beta_ci_hi = lik_miss$ci_hi[7],
      lik_M_g_y1_est = lik_miss$ests[5],
      # lik_M_g_y1_lo = lik_miss$ci_lo[5],
      # lik_M_g_y1_hi = lik_miss$ci_hi[5],
      lik_M_g_y2_est = lik_miss$ests[6],
      # lik_M_g_y2_lo = lik_miss$ci_lo[6],
      # lik_M_g_y2_hi = lik_miss$ci_hi[6],
      lik_F_beta_est = lik_full$ests[7],
      # lik_F_beta_ci_lo = lik_full$ci_lo[7],
      # lik_F_beta_ci_hi = lik_full$ci_hi[7],
      lik_F_g_y1_est = lik_full$ests[5],
      # lik_F_g_y1_lo = lik_full$ci_lo[5],
      # lik_F_g_y1_hi = lik_full$ci_hi[5],
      lik_F_g_y2_est = lik_full$ests[6],
      # lik_F_g_y2_lo = lik_full$ci_lo[6],
      # lik_F_g_y2_hi = lik_full$ci_hi[6],
      cox_beta_est = cox_full$ests[3],
      # cox_beta_ci_lo = cox_full$ci_lo[3],
      # cox_beta_ci_hi = cox_full$ci_hi[3],
      cox_g_y1_est = cox_full$ests[1],
      # cox_g_y1_lo = cox_full$ci_lo[1],
      # cox_g_y1_hi = cox_full$ci_hi[1],
      cox_g_y2_est = cox_full$ests[2]
      # cox_g_y2_lo = cox_full$ci_lo[2],
      # cox_g_y2_hi = cox_full$ci_hi[2]
    )
    
    # !!!!! DEBUGGING !!!!!
    if (T) {
      
      res$lik_M_a_x_est <- lik_miss$ests[1]
      res$lik_F_a_x_est <- lik_full$ests[1]
      
      res$lik_M_g_x1_est <- lik_miss$ests[2]
      res$lik_F_g_x1_est <- lik_full$ests[2]
      res$lik_M_g_x2_est <- lik_miss$ests[3]
      res$lik_F_g_x2_est <- lik_full$ests[3]
      
      res$lik_M_a_y_est <- lik_miss$ests[4]
      res$lik_F_a_y_est <- lik_full$ests[4]
      
      res$lik_M_a_v_est <- lik_miss$ests[8]
      res$lik_M_g_v1_est <- lik_miss$ests[9]
      res$lik_M_g_v2_est <- lik_miss$ests[10]
      
    }
    
  }
  
  return(res)
  
}
