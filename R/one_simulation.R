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
  
  # Set initial parameter estimate - should roughly (but not exactly) equal the
  # true parameters
  par <- log(c(
    a_x=0.003, g_x1=1.2, g_x2=1.1, a_y=0.002, g_y1=1.3, g_y2=1, beta_x=1.3,
    beta_z=0.8, t_x=0.999, t_y=1.001, a_s=0.03, t_s=1.001, g_s1=1.8, g_s2=1.7
  ))
  
  chk(2, "construct_negloglik_miss: START")
  negloglik_miss <- construct_negloglik_miss(dat)
  chk(2, "construct_negloglik_miss: END")
  chk(2, "optim: START")
  opt_miss <- optim(par=par, fn=negloglik_miss)
  chk(2, "optim: END")
  chk(2, "hessian: START")
  hessian_miss <- hessian(func=negloglik_miss, x=opt_miss$par)
  hessian_inv <- solve(hessian_miss)
  chk(2, "hessian: END")

  res <- list()
  for (i in c(1:length(par))) {
    res[[paste0("lik_M_",names(par)[i],"_est")]] <- as.numeric(opt_miss$par[i])
    res[[paste0("lik_M_",names(par)[i],"_se")]] <- sqrt(diag(hessian_inv))[i]
  }
  
  chk(3, "END")
  
  return(res)
  
}
