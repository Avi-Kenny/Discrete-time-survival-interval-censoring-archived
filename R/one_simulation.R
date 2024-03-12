#' Run a single simulation replicate
#'
#' @return A list containing the following:
#'     est_complete: HR estimate from complete data
#'     est_missing: HR estimate from MI procedure
#'     se_complete: SE of HR estimate from complete data
#'     se_missing: SE of HR estimate from MI procedure

one_simulation <- function() {
  
  # Set this flag to TRUE to speed up code (but with worse optim performance)
  speedy <- T
  if (speedy) {
    maxit <- 200
    hess_r <- 2
  } else {
    maxit <- 500
    hess_r <- 4
  }
  
  chk(0, "START")
  
  # Generate dataset
  dat <- generate_data(
    n = L$n,
    max_time = L$max_time,
    params = L$params
  )
  chk(1, "Data generated")
  
  # # Add x_prev column
  # # !!!!! This might not be needed
  # dat %<>% arrange(id, t_end)
  # dat$x_prev <- ifelse(
  #   dat$id==c(0,dat$id[c(1:(length(dat$id)-1))]),
  #   c(0,dat$x[c(1:(length(dat$x)-1))]),
  #   0
  # )
  
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
      Surv(t_start, t_end, y) ~ w_1 + w_2 + x + cluster(id),
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
  par_init <- log(c(
    # a_x=0.003, g_x1=1.2, g_x2=1.1, a_y=0.002, g_y1=1.3, g_y2=1, beta_x=1.3, # Monthly
    # beta_z=0.8, t_x=0.999, t_y=1.001, a_s=0.03, t_s=1.001, g_s1=1.8, g_s2=1.7 # Monthly
    a_x=0.03, g_x1=1.2, g_x2=1.1, a_y=0.02, g_y1=1.3, g_y2=1, beta_x=1.3, # Yearly
    beta_z=0.8, t_x=0.999, t_y=1.001, a_s=0.03, t_s=1.001, g_s1=1.8, g_s2=1.7 # Yearly
  ))
  
  chk(2, "construct_negloglik_miss: START")
  negloglik_miss <- construct_negloglik_miss(dat)
  chk(2, "construct_negloglik_miss: END")
  chk(2, "optim: START")
  opt_miss <- stats::optim(
    par = par_init,
    fn = negloglik_miss,
    method = "Nelder-Mead",
    control = list(maxit=maxit)
  )
  chk(2, "optim: END")
  chk(2, "hessian: START")
  hessian_miss <- numDeriv::hessian(
    func = negloglik_miss,
    x = opt_miss$par,
    method = "Richardson",
    method.args = list(r=hess_r)
  )
  hessian_inv <- solve(hessian_miss)
  chk(2, "hessian: END")
  
  res <- list()
  pn <- names(par_init)
  for (i in c(1:length(par_init))) {
    res[[paste0("lik_M_",pn[i],"_est")]] <- as.numeric(opt_miss$par[i])
    res[[paste0("lik_M_",pn[i],"_se")]] <- sqrt(diag(hessian_inv))[i]
  }
  
  chk(3, "END")
  
  return(res)
  
}
