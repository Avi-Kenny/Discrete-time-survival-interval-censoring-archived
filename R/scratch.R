

# Checking source of TEMP ERROR A1
if (F) {
  
  L <- list(
    n = 500,
    max_time = 70,
    params = list(
      a_x=log(0.005), a_y=log(0.003), a_v=log(0.7), g_x=c(0,0),
      g_y=c(0,0), g_v=c(0,0), beta=log(1.5)
    )
    
    # n = 400,
    # max_time = 50,
    # params = list(
    #   a_x=log(0.005), a_y=log(0.003), a_v=log(0.1),
    #   g_x=c(log(1.3),log(1.002)), g_y=c(log(1.2),log(1.001)),
    #   g_v=c(log(1.2),log(1.001)), beta=log(1.5)
    # )
  )
  
  for (k in c(1:100)) {
    
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
    
    # Convert parameter vector to a named list
    p <- log(c(0.002, 1, 1, 0.002, 1, 1, 1, 0.001, 1, 1))
    params <- list(a_x=p[1], g_x=c(p[2],p[3]), a_y=p[4], g_y=c(p[5],p[6]),
                   beta=p[7], a_v=p[8], g_v=c(p[9],p[10]))
    
    # Compute the negative likelihood across individuals
    n <- attr(dat, "n")
    TMP <- lapply(c(1:n), function(i) {
      dat_i <- filter(dat, id==i)
      J <- nrow(dat_i)
      x <- dat_i$x
      x_prev <- c(0, x[c(1:(J-1))])
      if (!all(x_prev==dat_i$x_prev)) {
        print("TEMP Error A1")
        browser()
      }
      return(0)
    })
    
  }
  
  
}


# Misc
if (F) {
  
  L <- list(
    n = 50,
    max_time = 100,
    params = list(g_x = c(log(1.3),log(1.002)),
                  g_y = c(log(1.2),log(1.001)),
                  g_v = c(log(1.2),log(1.001)),
                  beta = log(1.5))
  )
  
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
  
  # Testing the likelihood functions
  par <- log(c(0.002, 1, 1, 0.002, 1, 1, 1, 0.001, 1, 1))
  negloglik_miss(dat, par)
  
  # par <- log(c(0.002, 1, 1, 0.002, 1, 1, 1, 0.001, 1, 1))
  # fnc_miss <- function(par) { negloglik_miss(dat, par) }
  # opt_miss <- optim(par=par, fn=fnc_miss)
  # hessian_miss <- optimHess(par=opt_miss$par, fn=fnc_miss)
  # lik_miss <- list(ests=exp(opt_miss$par))
  # lik_miss$se <- sqrt(diag(solve(hessian_miss)))
  # lik_miss$ci_lo <- exp(lik_miss$ests-1.96*lik_miss$se)
  # lik_miss$ci_hi <- exp(lik_miss$ests+1.96*lik_miss$se)
  
}
