#' Negative log-likelihood across individuals and time
#' @param dat A dataset returned by generate_dataset()
#' @param par Vector of parameters governing the distribution.
#' @return Numeric likelihood
#' @notes This corresponds to the missing data structure
construct_negloglik_miss <- function(dat, parallelize=FALSE, cl=NULL) {
  
  # Construct a data object specific to each individual
  n <- attr(dat, "n")
  dat_objs <- lapply(c(1:n), function(i) {
    
    d <- list()
    dat_i <- dat[dat$id==i,] # Does dat_i need to be returned with r?
    
    # Start and end times
    d$s_i <- attr(dat, "s_i")[i]
    d$t_i <- attr(dat, "t_i")[i]
    
    # Data vectors
    d$w <- subset(dat_i, select=names(dat)[substr(names(dat), 1, 2)=="w_"])
    d$y <- dat_i$y
    d$z <- dat_i$z
    d$v <- dat_i$v
    d$d <- dat_i$d
    d$u <- dat_i$u
    d$cal_time <- dat_i$t_end
    d$cal_time_sc <- d$cal_time / 100 # Rescaling is done to prevent optimization issues
    
    # Calculate the set X_i to sum over
    d$X_i_set <- list()
    for (j in c(1:(d$t_i-d$s_i+2))) {
      x_ <- c(rep(0,d$t_i-d$s_i-j+2), rep(1,j-1))
      case_i_ <- case(x_,d$v)
      T_pm <- T_plusminus(case=case_i_, s_i=d$s_i, t_i=d$t_i, x=x_, v=d$v)
      if (case_i_!=9 &&
          all(d$u==x_*d$d) &&
          all(d$d==g_delta(case=case_i_, s_i=d$s_i, t_i=d$t_i,
                           T_minus=T_pm$T_minus, T_plus=T_pm$T_plus))
      ) {
        d$X_i_set <- c(d$X_i_set, list(x_))
      }
    }
    
    return(d)
    
  })
  
  negloglik_miss <- function(par) {
    
    print(paste("negloglik_miss() called:",Sys.time())) # !!!!!
    # fn_calls <<- fn_calls+1 # !!!!!

    # Convert parameter vector to a named list
    p <- as.numeric(par)
    params <- list(a_x=p[1], g_x=c(p[2],p[3]), a_y=p[4], g_y=c(p[5],p[6]),
                   beta_x=p[7], beta_z=p[8], t_x=p[9], t_y=p[10],
                   a_s=p[11], t_s=p[12], g_s=c(p[13],p[14]))
    
    lik_fn <- function(i) {
      
      # Extract data for individual i
      d <- dat_objs[[i]]
      
      # Compute the likelihood for individual i
      w <- t(d$w)
      f2 <- sum(unlist(lapply(d$X_i_set, function(x) {
        x_prev <- c(0,x[1:(length(x)-1)])
        prod(unlist(lapply(c(1:(d$t_i-d$s_i+1)), function(j) {
          # Note: here, j is the index of time within an individual rather than
          #       a calendar time index
          w_ij <- as.numeric(w[,j])
          return(
            f_x(x=x[j], x_prev=x_prev[j], w=w_ij, j=d$cal_time_sc[j],
                s=In(d$cal_time[j]==d$s_i), params=params) * # Maybe create this indicator variable (vector) earlier
              f_y(y=d$y[j], x=x[j], w=w_ij, z=d$z[j], j=d$cal_time_sc[j],
                  params=params)
          )
        })))
      })))
      # if (is.nan(f2) || is.nan(f2)) { browser() } # Debugging
      if (f2<=0) {
        f2 <- 1e-10
        # warning("Likelihood of zero")
      }
      
      return(f2)
      
    }
    
    # !!!!! Debugging
    if (F) {
      print("Debugging: START")
      print("cl")
      print(cl)
      
      n_batches <- 4
      folds <- cut(c(1:n), breaks=n_batches, labels=FALSE)
      batches <- lapply(c(1:n_batches), function(batch) {
        c(1:n)[folds==batch]
      })
      # cl <- parallel::makeCluster(10) # !!!!!
      parallel::clusterExport(cl, c("batches", "lik_fn"), envir=environment())
      
      fnc <- function(x) { Sys.sleep(4) }
      # fnc <- function(x) { Sys.sleep(0.0001) }
      t_1 <- system.time({
        r1 <- lapply(c(1:70), fnc)
      })
      t_2 <- system.time({
        r2 <- parallel::parLapply(cl, c(1:70), fnc)
      })
      print("t_1")
      print(t_1)
      print("t_2")
      print(t_2)
      stop("hey")
      
      
      t1 <- system.time({
        v1 <- -1 * sum(log(unlist(lapply(c(1:n), lik_fn))))
      })
      t2 <- system.time({
        v2 <- -1 * sum(log(unlist(parallel::parLapply(cl, c(1:n), lik_fn))))
      })
      t3 <- system.time({
        v3 <- -1 * sum(unlist(parallel::parLapply(
          cl,
          c(1:n_batches),
          function(batch) {
            sum(log(unlist(lapply(batches[[batch]], lik_fn))))
          })
        ))
      })
      
      print("Serial")
      print(v1)
      print(t1)
      print("Parallel")
      print(v2)
      print(t2)
      print(paste0("Parallel (", n_batches, " batches)"))
      print(v3)
      print(t3)
      stop("TEMP STOP")
      
      # n_cores <- parallel::detectCores() - 1
      # print(paste0("Using ", n_cores, " cores."))
      # cl <- parallel::makeCluster(n_cores)
      # parallel::clusterExport(cl, ls(.GlobalEnv))
      
    }

    # Compute the negative likelihood across individuals
    if (parallelize) {
      
      n_batches <- length(cl)
      folds <- cut(c(1:n), breaks=n_batches, labels=FALSE)
      batches <- lapply(c(1:n_batches), function(batch) {
        c(1:n)[folds==batch]
      })
      parallel::clusterExport(cl, c("batches"), envir=environment())
      return(-1 * sum(unlist(parallel::parLapply(
        cl,
        c(1:n_batches),
        function(batch) { sum(log(unlist(lapply(batches[[batch]], lik_fn)))) })
      )))
      
      # return(-1 * sum(log(unlist(parallel::parLapply(cl, c(1:n), lik_fn)))))
      
    } else {
      
      return(-1 * sum(log(unlist(lapply(c(1:n), lik_fn)))))
      
    }
    
  }
  
  return(negloglik_miss)
  
}



#' Calculate likelihood component f_x
#'
#' @param x Seroconversion indicator (time j)
#' @param x_prev Seroconversion indicator (time j-1)
#' @param w Vector of covariates (time j)
#' @param params Named list of parameters
#' @return Numeric likelihood
f_x <- function(x, x_prev, w, j, s, params) {
  if (s==0) {
    if (x==1) {
      if (x_prev==1) {
        return(1)
      } else {
        return(exp(params$a_x + params$t_x*j + sum(params$g_x*w)))
      }
    } else {
      if (x_prev==1) {
        return(0)
      } else {
        return(1 - exp(params$a_x + params$t_x*j + sum(params$g_x*w)))
      }
    }
  } else {
    expitlin <- expit(params$a_s + params$t_s*j + sum(params$g_s*w))
    if (x==1) {
      return(expitlin)
    } else {
      return(1-expitlin)
    }
  }
}



#' Calculate likelihood component f_y
#'
#' @param y Event indicator (time j)
#' @param x Seroconversion indicator (time j)
#' @param w Vector of covariates (time j)
#' @param z ART indicator (time j)
#' @param params Named list of parameters
#' @return Numeric likelihood
f_y <- function(y, x, w, z, j, params) {
  explin <- exp(params$a_y + params$t_y*j + sum(params$g_y*w) + params$beta_x*x +
                  params$beta_z*z)
  if (y==1) { return(explin) } else { return(1-explin) }
}
