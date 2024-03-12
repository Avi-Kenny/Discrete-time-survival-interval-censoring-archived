
# New DGM parameters
if (F) {
  
  # a_s = 0.05, # Yearly
  # a_x = 0.005, # Yearly
  # a_y = 0.003, # Yearly
  # a_z = 0.01, # Yearly
  
  new_alpha <- function(alpha) { log(1-(1-exp(alpha))^12) }
  
  old_a <- log(0.01)
  new_a <- new_alpha(old_a)
  
  print(paste("Prob (mth):", round(exp(old_a),4)))
  print(paste("Prob (yr):", round(exp(new_a),4)))
  
}

# Tweaking Hessian tuning params
if (F) {
  
  negloglik_miss <- function(x) {
    Sys.sleep(0.1)
    (x[1]-1)^2 + (x[2]-2)^2 + (x[3]-3)^2
  }
  
  # No Hessian
  system.time({
    opt_miss <- stats::optim(
      par = c(0,0,0),
      fn = negloglik_miss,
      method = "Nelder-Mead",
      control = list(maxit=500,
                     reltol=1e-8)
    )
    print(opt_miss)
  })
  
  # Hessian (single call)
  system.time({
    opt_miss <- stats::optim(
      par = c(0,0,0),
      fn = negloglik_miss,
      method = "Nelder-Mead",
      control = list(maxit=500,
                     reltol=1e-8),
      hessian = T
    )
    print(opt_miss)
  })
  
  # Hessian (separate calls)
  system.time({
    opt_miss <- stats::optim(
      par = c(0,0,0),
      fn = negloglik_miss,
      method = "Nelder-Mead",
      control = list(maxit=20,
                     reltol=1e-8)
    )
    print(opt_miss)
    
  })
  system.time({
    opt_hess <- stats::optimHess(
      par = opt_miss$par,
      fn = negloglik_miss,
      control = list(maxit=20,
                     reltol=0.01)
    )
    print(opt_hess)
    
  })
  
  system.time({
    hessian_miss <- numDeriv::hessian(
      func = negloglik_miss,
      x = opt_miss$par,
      method = "Richardson",
      method.args = list(
        eps = 0.0001,
        d = 0.0001, # d gives the fraction of x to use for the initial numerical approximation
        zero.tol = sqrt(.Machine$double.eps/7e-7),
        r = 4, # r gives the number of Richardson improvement iterations (repetitions with successly smaller d
        v = 2, # v gives the reduction factor
        show.details = F
      )
    )
    print(round(hessian_miss, 4))
  })
  
}

# TEMP debugging (new link function)
if (F) {
  
  grid <- seq(-3,0.5,0.001)
  expit <- function(x) { 1/(1+exp(-x)) }
  logit <- function(x) { log(x/(1-x)) }
  e <- -0.1
  ell <- logit(exp(e))
  x_0 <- e - (ell*exp(ell))/(exp(e)*(1+exp(ell))^2)
  k_0 <- exp(e-ell)*(1+exp(ell))^2
  exp2 <- Vectorize(function(x) {
    In(x<=e) * exp(x) +
    In(x>e) * expit(k_0*(x-x_0))
  })
  df <- data.frame(
    x = rep(grid,2),
    y = c(sapply(grid,exp), sapply(grid,exp2)),
    which = rep(c("exp", "exp2"), each=length(grid))
  )
  ggplot(df, aes(x=x, y=y, color=which)) + geom_line()
  ggplot(df, aes(x=x, y=y, color=which)) + geom_line() +
    lims(x=c(-0.3,0.1), y=c(0.7,1.1))
  
  cloglog <- function(x) { log(-log(1-x)) }
  cloglog_i <- function(x) { 1 - exp(-exp(x)) }
  df <- data.frame(
    x = rep(grid,3),
    y = c(sapply(grid,exp), sapply(grid,exp2), sapply(grid,cloglog_i)),
    which = rep(c("exp", "exp2", "cloglog_i"), each=length(grid))
  )
  ggplot(df, aes(x=x, y=y, color=which)) + geom_line()
  ggplot(df, aes(x=x, y=y, color=which)) + geom_line() +
    lims(x=c(-0.3,0.1), y=c(0.7,1.1))
  
  
}

# Plotting parameter profiles
if (F) {
  
  # Try 2D plots of parameter surface
  
  fnc_miss(opt_miss$par) # 1007.201
  
  # P1 49.88202
  # P2 0.004576096
  # P3 9.425506
  # P4 25.86428
  # P5 0.003669396
  # P6 7.77546
  # P7 6.810297
  
  # Param 1
  {
    grid <- seq(-6.33825644, -4.33825644, length.out=11)
    evals <- sapply(grid, function(p1) {
      fnc_miss(c(p1, opt_miss$par[2:7]))
    })
    df <- data.frame(x=grid, y=evals)
    lm_fit <- lm(y~x+I(x^2), data=df)
    cf <- lm_fit$coefficients
    fn_evals <- sapply(grid, function(x) {
      cf[1] + cf[2]*x + cf[3]*x^2
    })
    ggplot(df, aes(x=x, y=y)) +
      geom_point() +
      geom_line(data=data.frame(x=grid, y=fn_evals), color="green")
    print(cf[3]) # 49.88202
  }
  
  # Param 7
  {
    grid <- seq(-0.5, 0.5, length.out=11)
    evals <- sapply(grid, function(p7) {
      fnc_miss(c(opt_miss$par[1:6], p7))
    })
    df <- data.frame(x=grid, y=evals)
    lm_fit <- lm(y~x+I(x^2), data=df)
    cf <- lm_fit$coefficients
    fn_evals <- sapply(grid, function(x) {
      cf[1] + cf[2]*x + cf[3]*x^2
    })
    ggplot(df, aes(x=x, y=y)) +
      geom_point() +
      geom_line(data=data.frame(x=grid, y=fn_evals), color="green")
    print(cf[3]) # 6.810297
  }
  
  # Param k
  {
    grid_lims <- list(c(0.002,0.008),c(1,1.5),c(1.2,2),c(0.001,0.005),c(0.7,1.4),c(1,2),c(1,1.7))
    grid_lims <- lapply(grid_lims,log)
    k <- 6
    grid <- seq(grid_lims[[k]][1], grid_lims[[k]][2], length.out=9)
    evals <- sapply(grid, function(p) {
      prm <- replace(opt_miss$par, k, p)
      fnc_miss(prm)
    })
    df <- data.frame(x=grid, y=evals)
    lm_fit <- lm(y~x+I(x^2), data=df)
    cf <- lm_fit$coefficients
    fn_evals <- sapply(grid, function(x) {
      cf[1] + cf[2]*x + cf[3]*x^2
    })
    ggplot(df, aes(x=x, y=y)) +
      geom_point() +
      geom_line(data=data.frame(x=grid, y=fn_evals), color="green") +
      labs(title=paste("Param",k))
    print(cf[3])
  }
  
}

# Try reading in AHRI dataset
if (F) {
  
  library(readstata13)
  dat <- read.dta13("../AHRI Surveillance Datasets/SurveillanceEpisodes.dta")
  
  dat2 <- dat[,1:9]
  
  # Testing whether years are "skipped"
  for (i in c(2:500)) {
    year_prev <- dat2[round(i-1),"CalendarYear"]
    year_curr <- dat2[i,"CalendarYear"]
    ind_prev <- dat2[round(i-1),"IndividualId"]
    ind_curr <- dat2[i,"IndividualId"]
    if (ind_prev==ind_curr && year_curr-year_prev>1) {
      print(paste(year_prev,year_curr,ind_prev,ind_curr))
    }
  }
  
  # Number of unique ind-year combinations
  # About 2 million
  dat2 %<>% mutate(
    key = paste0(IndividualId,"-",CalendarYear)
  )
  length(unique(dat2$key))
  length(unique(dat2$IndividualId))
  
}

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
