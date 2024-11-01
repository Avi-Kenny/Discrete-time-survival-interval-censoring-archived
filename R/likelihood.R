#' Transform dataset
#' @param dat A dataset returned by generate_dataset()
#' @param model_version Model version, as specified in MAIN.R
transform_dataset <- function(dat, model_version=0, window_start, window_end) {
  
  # Construct spline bases
  b9 <- construct_basis("age (13,20,30,40,60)")
  b10 <- construct_basis("year (10,13,16,19,22)", window_start=window_start,
                         window_end=window_end)
  b12 <- construct_basis("year (17,...,22)", window_start=window_start,
                         window_end=window_end)
  
  # This procedure assumes the dataset is sorted by id
  row_map <- new.env()
  id_prev <- -1
  i <- 1
  row_start <- 1
  for (id in dat$id) {
    if (id!=id_prev) {
      row_map[[as.character(id_prev)]] <- c(row_start:(i-1))
      row_start <- i
    }
    id_prev <- id
    i <- i + 1
  }
  row_map[[as.character(id_prev)]] <- c(row_start:(i-1))
  
  # Extract and remove attributes
  n <- attr(dat, "n")
  s_i <- attr(dat, "s_i")
  t_i <- attr(dat, "t_i")
  attr(dat, "n") <- NULL
  attr(dat, "s_i") <- NULL
  attr(dat, "t_i") <- NULL
  
  # Construct a dataframe specific to each individual
  dat_objs <- lapply(c(1:n), function(i) {
    
    d <- list()
    d$dat_i <- dat[row_map[[as.character(i)]],]
    
    # Start and end times
    d$s_i <- s_i[i]
    d$t_i <- t_i[i]
    
    # Create new variables
    d$dat_i$init_visit <- 0
    d$dat_i[1,"init_visit"] <- 1
    d$dat_i$cal_time_sc <- d$dat_i$t_end / 10 # Rescaling originally done to prevent optimization issues; check if this is still needed
    
    # Apply spline bases to dataframe
    if (model_version %in% c(29:33)) {
      d$dat_i$b9_1 <- signif(sapply(d$dat_i$w_1, function(w_1) { b9(w_1,1) }),4)
      d$dat_i$b9_2 <- signif(sapply(d$dat_i$w_1, function(w_1) { b9(w_1,2) }),4)
      d$dat_i$b9_3 <- signif(sapply(d$dat_i$w_1, function(w_1) { b9(w_1,3) }),4)
      d$dat_i$b9_4 <- signif(sapply(d$dat_i$w_1, function(w_1) { b9(w_1,4) }),4)
      d$dat_i$b10_1 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b10(j,1) }),4)
      d$dat_i$b10_2 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b10(j,2) }),4)
      d$dat_i$b10_3 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b10(j,3) }),4)
      d$dat_i$b10_4 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b10(j,4) }),4)
    } else if (model_version==34) {
      d$dat_i$b9_1 <- signif(sapply(d$dat_i$w_1, function(w_1) { b9(w_1,1) }),4)
      d$dat_i$b9_2 <- signif(sapply(d$dat_i$w_1, function(w_1) { b9(w_1,2) }),4)
      d$dat_i$b9_3 <- signif(sapply(d$dat_i$w_1, function(w_1) { b9(w_1,3) }),4)
      d$dat_i$b9_4 <- signif(sapply(d$dat_i$w_1, function(w_1) { b9(w_1,4) }),4)
      d$dat_i$b12_1 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b12(j,1) }),4)
      d$dat_i$b12_2 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b12(j,2) }),4)
      d$dat_i$b12_3 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b12(j,3) }),4)
      d$dat_i$b12_4 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b12(j,4) }),4)
    }
    
    # Calculate the set X_i to sum over
    d$X_i_set <- list()
    for (j in c(1:(d$t_i-d$s_i+2))) {
      
      x_ <- c(rep(0,d$t_i-d$s_i-j+2), rep(1,j-1))
      T_pm <- T_plusminus(s_i=d$s_i, t_i=d$t_i, x=x_, v=d$dat_i$v)
      case_i_ <- case(T_pm$T_minus, T_pm$T_plus)
      
      if (case_i_!=9 &&
          all(d$dat_i$u==x_*d$dat_i$delta) &&
          all(d$dat_i$delta==g_delta(case=case_i_, s_i=d$s_i, t_i=d$t_i,
                                     T_minus=T_pm$T_minus, T_plus=T_pm$T_plus))
      ) {

        if (length(x_)==1) {
          x_prev <- 0
        } else {
          x_prev <- c(0,x_[1:(length(x_)-1)])
        }
        
        d$X_i_set <- c(d$X_i_set, list(list(
          x = sum(x_),
          x_prev = max(sum(x_)-1,0)
        )))
        
      }
      
    }
    
    # # Convert to matrix
    # d$dat_i <- as.matrix(d$dat_i)
    
    return(d)
    
  })
  
  rm(dat)
  
  return(dat_objs)
  
}

#' Negative log-likelihood across individuals and time
#' @param parallelize Whether or not to parallelize computation
#' @param model_version Model version, as specified in MAIN.R
#' @return Numeric likelihood
#' @notes This corresponds to the missing data structure
construct_negloglik <- function(
  inds, parallelize=FALSE, model_version=0, use_counter=F, temp=FALSE
) {
  
  if (!identical(temp, FALSE)) { dat_objs_wrapper <- temp }
  
  # cl <- parallel::makeCluster(cfg$sim_n_cores)
  # objs_to_export <- c("f_x", "f_y", "icll", "lik_fn2", "inds", "batches",
  #                     "uncompress")
  # parallel::clusterExport(cl, objs_to_export, envir=.GlobalEnv)
  
  negloglik <- function(par) {
    
    if (use_counter) {
      counter <<- counter + 1
      if (counter%%10==0) {
        print(paste0(counter, " negloglik calls processed: ",Sys.time()))
      }
    }

    # Convert parameter vector to a named list
    p <- as.numeric(par)
    if (model_version==7) {
      params <- list(a_x=p[1], g_x=p[2:3], t_x=p[4], a_s=p[5], g_s=p[6:7], t_s=p[8], beta_x=p[9], a_y=p[10], g_y=p[11:12], t_y=p[13])
    } else if (model_version==29) {
      params <- list(a_x=p[1], g_x=p[2:5], t_x=p[6:9], a_s=p[10], g_s=p[11:14], beta_x=p[15:20], a_y=p[21], g_y=p[22:25], t_y=p[26:29])
    } else if (model_version==30) {
      params <- list(a_x=p[1], g_x=p[2:5], t_x=p[6:9], a_s=p[10], g_s=p[11:14], beta_x=p[15:18], a_y=p[19], g_y=p[20:23], t_y=p[24:27])
    } else if (model_version==31) {
      params <- list(a_x=p[1], g_x=p[2:5], t_x=p[6:9], a_s=p[10], g_s=p[11:14], beta_x=p[15:20], a_y=p[21], g_y=p[22:25], t_y=p[26:29])
    } else if (model_version==32) {
      params <- list(a_x=p[1], g_x=p[2:5], t_x=p[6:9], a_s=p[10], g_s=p[11:14], beta_x=p[15:23], a_y=p[24], g_y=p[25:28], t_y=p[29:32])
    } else if (model_version==33) {
      params <- list(a_x=p[1], g_x=p[2:5], t_x=p[6:9], a_s=p[10], g_s=p[11:14], beta_x=p[15:23], a_y=p[24], g_y=p[25:28], t_y=p[29:32])
    } else if (model_version==34) {
      params <- list(a_x=p[1], g_x=p[2:5], t_x=p[6:9], a_s=p[10], g_s=p[11:14], beta_x=p[15:18], a_y=p[19], g_y=p[20:23], t_y=p[24:27])
    }
    
    # Compute the negative likelihood across individuals
    if (parallelize) {
      
      # Original code
      nll <- -1 * sum(log(unlist(
        parallel::parLapply(cl, dat_objs_wrapper, function(d) {
          lapply(d, function(d2) { lik_fn2(d2, params, inds) })
        })
      )))
      return(nll)
      
      if (F) {
        # # Debugging 2
        # nll <- -1 * sum(log(unlist(
        #   parallel::parLapply(cl, dat_objs_wrapper, function(d) {
        #     lapply(d, function(d2) { Sys.sleep(0.0001); return(1); })
        #   })
        # )))
        # return(nll)
        
        # # !!!!! DEBUGGING 1
        # if (F) {
        #   nll <- -1 * sum(log(unlist(
        #     lapply(dat_objs_wrapper, function(d) {
        #       st <- system.time({
        #         x <- lapply(d, function(d2) { lik_fn2(d2, params, inds) })
        #       })
        #       print(paste0("Length(d): ", length(d)))
        #       print("Time:")
        #       print(st)
        #     })
        #   )))
        #   return(nll)
        # }
      }

    } else {
      
      nll <- -1 * sum(log(unlist(
        lapply(dat_objs_wrapper, function(d) {
          lapply(d, function(d2) { lik_fn2(d2, params, inds) })
        })
      )))
      return(nll)
      
    }
    
  }
  
  return(negloglik)
  
}



#' Calculate likelihood for individual i
#'
#' @param d One component of the list dat_objs
#' @param params TO DO
#' @param inds TO DO
#' @return Numeric likelihood
#' @note dat_objs is accessed globally
lik_fn2 <- function(d, params, inds) {
  
  {
    
    # Precompute values
    f_xy_vals <- apply(X=d$dat_i, MARGIN=1, FUN = function(r) {
      w_ij <- as.numeric(r[inds$w])
      j <- r[["cal_time_sc"]]
      spl_ij <- r[inds$spl]
      
      f_x_00 <- f_x(x=0, x_prev=0, w=w_ij, j=j, s=r[["init_visit"]], spl=spl_ij,
                    params=params)
      f_x_01 <- f_x(x=1, x_prev=0, w=w_ij, j=j, s=r[["init_visit"]], spl=spl_ij,
                    params=params)
      f_x_11 <- f_x(x=1, x_prev=1, w=w_ij, j=j, s=r[["init_visit"]], spl=spl_ij,
                    params=params)
      f_y_0 <- f_y(y=r[["y"]], x=0, w=w_ij, z=r[["z"]], j=j, spl=spl_ij,
                   params=params)
      f_y_1 <- f_y(y=r[["y"]], x=1, w=w_ij, z=r[["z"]], j=j, spl=spl_ij,
                   params=params)
      
      return(c(f_x_00*f_y_0, f_x_01*f_y_1, f_x_11*f_y_1))
    })
    
    len_obs <- round(d$t_i - d$s_i + 1)
    
    f2_fnc <- function(x_obj) {
      
      x <- uncompress(len_obs, x_obj$x)
      x_prev <- uncompress(len_obs, x_obj$x_prev)
      
      return(prod(sapply(X=c(1:length(x)), FUN = function(j) {
        x_prev_j <- x_prev[j]
        x_j <- x[j]
        if (x_prev_j==0 && x_j==0) {
          return(f_xy_vals[1,j])
        } else if (x_prev_j==0 && x_j==1) {
          return(f_xy_vals[2,j])
        } else if (x_prev_j==1 && x_j==1) {
          return(f_xy_vals[3,j])
        }
      })))
      
    }
    
    f2 <- sum(unlist(lapply(d$X_i_set, f2_fnc)))
      
  }
  
  # if (is.na(f2) || is.nan(f2)) {
  #   browser()
  # } # Debugging
  if (f2<=0) {
    f2 <- 1e-10
    # warning("Likelihood of zero")
  }
  
  return(f2)
  
}



#' Calculate likelihood component f_x
#'
#' @param x Seroconversion indicator (time j)
#' @param x_prev Seroconversion indicator (time j-1)
#' @param w Vector of covariates (time j)
#' @param params Named list of parameters
#' @return Numeric likelihood

if (cfg$model_version==7) {
  
  f_x <- function(x, x_prev, w, j, s, spl, params) {
    if (s==0) {
      if (x==1 && x_prev==1) {
        return(1)
      } else {
        prob <- icll(params$a_x + params$t_x*j + sum(params$g_x*w))
        if (x==1) { return(prob) } else { return(1-prob) }
      }
    } else {
      prob <- icll(params$a_s + params$t_s*j + sum(params$g_s*w))
      if (x==1) { return(prob) } else { return(1-prob) }
    }
  }
  
} else if (cfg$model_version %in% c(29:33)) {
  
  f_x <- function(x, x_prev, w, j, s, spl, params) {
    if (s==0) {
      if (x==1 && x_prev==1) {
        return(1)
      } else {
        prob <- icll(
          params$a_x + params$t_x[1]*spl[["b10_1"]] +
            params$t_x[2]*spl[["b10_2"]] + params$t_x[3]*spl[["b10_3"]] +
            params$t_x[4]*spl[["b10_4"]] +
            params$g_x[1]*spl[["b9_1"]] + params$g_x[2]*spl[["b9_2"]] +
            params$g_x[3]*spl[["b9_3"]] + params$g_x[4]*spl[["b9_4"]]
        )
        if (x==1) { return(prob) } else { return(1-prob) }
      }
    } else {
      prob <- icll(
        params$a_s + params$g_s[1]*spl[["b9_1"]] + params$g_s[2]*spl[["b9_2"]] +
          params$g_s[3]*spl[["b9_3"]] + params$g_s[4]*spl[["b9_4"]]
      )
      if (x==1) { return(prob) } else { return(1-prob) }
    }
  }
  
} else if (cfg$model_version==34) {
  
  f_x <- function(x, x_prev, w, j, s, spl, params) {
    if (s==0) {
      if (x==1 && x_prev==1) {
        return(1)
      } else {
        prob <- icll(
          params$a_x + params$t_x[1]*spl[["b12_1"]] +
            params$t_x[2]*spl[["b12_2"]] + params$t_x[3]*spl[["b12_3"]] +
            params$t_x[4]*spl[["b12_4"]] +
            params$g_x[1]*spl[["b9_1"]] + params$g_x[2]*spl[["b9_2"]] +
            params$g_x[3]*spl[["b9_3"]] + params$g_x[4]*spl[["b9_4"]]
        )
        if (x==1) { return(prob) } else { return(1-prob) }
      }
    } else {
      prob <- icll(
        params$a_s + params$g_s[1]*spl[["b9_1"]] + params$g_s[2]*spl[["b9_2"]] +
          params$g_s[3]*spl[["b9_3"]] + params$g_s[4]*spl[["b9_4"]]
      )
      if (x==1) { return(prob) } else { return(1-prob) }
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

if (cfg$model_version==7) {
  
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      params$a_y + params$t_y*j + sum(params$g_y*w) + params$beta_x*x
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
} else if (cfg$model_version==29) {
  
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      x*(
        (params$beta_x[1] + params$beta_x[2]*j) +
          (params$beta_x[3] + params$beta_x[4]*j) * max(w[1]-0.3,0) +
          (params$beta_x[5] + params$beta_x[6]*j) * max(w[1]-0.45,0)
      ) +
        params$a_y + params$g_y[1]*spl[["b9_1"]] +
        params$g_y[2]*spl[["b9_2"]] + params$g_y[3]*spl[["b9_3"]] +
        params$g_y[4]*spl[["b9_4"]] + params$t_y[1]*spl[["b10_1"]] +
        params$t_y[2]*spl[["b10_2"]] + params$t_y[3]*spl[["b10_3"]] +
        params$t_y[4]*spl[["b10_4"]]
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
} else if (cfg$model_version==30) {
  
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      x*(params$beta_x[1] + params$beta_x[2]*j +
           params$beta_x[3]*w[1] + params$beta_x[4]*j*w[1]) +
        params$a_y + params$g_y[1]*spl[["b9_1"]] +
        params$g_y[2]*spl[["b9_2"]] + params$g_y[3]*spl[["b9_3"]] +
        params$g_y[4]*spl[["b9_4"]] + params$t_y[1]*spl[["b10_1"]] +
        params$t_y[2]*spl[["b10_2"]] + params$t_y[3]*spl[["b10_3"]] +
        params$t_y[4]*spl[["b10_4"]]
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
} else if (cfg$model_version==31) {
  
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      x*(
        (params$beta_x[1] + params$beta_x[2]*j) +
          (params$beta_x[3] + params$beta_x[4]*j) * max(w[1]-0.2,0) +
          (params$beta_x[5] + params$beta_x[6]*j) * min(max(w[1]-0.4,0), 0.5)
      ) +
        params$a_y + params$g_y[1]*spl[["b9_1"]] +
        params$g_y[2]*spl[["b9_2"]] + params$g_y[3]*spl[["b9_3"]] +
        params$g_y[4]*spl[["b9_4"]] + params$t_y[1]*spl[["b10_1"]] +
        params$t_y[2]*spl[["b10_2"]] + params$t_y[3]*spl[["b10_3"]] +
        params$t_y[4]*spl[["b10_4"]]
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
} else if (cfg$model_version==32) {
  
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      x*(
        (params$beta_x[1] + params$beta_x[2]*j + params$beta_x[3]*max(j-0.6,0)) +
          w[1] * (
            params$beta_x[4] + params$beta_x[5]*j + params$beta_x[6]*max(j-0.6,0)
          ) +
          max(w[1]-0.4,0) * (
            params$beta_x[7] + params$beta_x[8]*j + params$beta_x[9]*max(j-0.6,0)
          )
      ) +
        params$a_y + params$g_y[1]*spl[["b9_1"]] +
        params$g_y[2]*spl[["b9_2"]] + params$g_y[3]*spl[["b9_3"]] +
        params$g_y[4]*spl[["b9_4"]] + params$t_y[1]*spl[["b10_1"]] +
        params$t_y[2]*spl[["b10_2"]] + params$t_y[3]*spl[["b10_3"]] +
        params$t_y[4]*spl[["b10_4"]]
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
} else if (cfg$model_version==33) {
  
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      x*(
        1 * (
          params$beta_x[1] + params$beta_x[2]*j + params$beta_x[3]*j^2
        ) +
          w[1] * (
            params$beta_x[4] + params$beta_x[5]*j + params$beta_x[6]*j^2
          ) +
          (w[1])^2 * (
            params$beta_x[7] + params$beta_x[8]*j + params$beta_x[9]*j^2
          )
      ) +
        params$a_y + params$g_y[1]*spl[["b9_1"]] +
        params$g_y[2]*spl[["b9_2"]] + params$g_y[3]*spl[["b9_3"]] +
        params$g_y[4]*spl[["b9_4"]] + params$t_y[1]*spl[["b10_1"]] +
        params$t_y[2]*spl[["b10_2"]] + params$t_y[3]*spl[["b10_3"]] +
        params$t_y[4]*spl[["b10_4"]]
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
} else if (cfg$model_version==34) {
  
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      x*(params$beta_x[1] + params$beta_x[2]*j +
           params$beta_x[3]*w[1] + params$beta_x[4]*j*w[1]) +
        params$a_y + params$g_y[1]*spl[["b9_1"]] +
        params$g_y[2]*spl[["b9_2"]] + params$g_y[3]*spl[["b9_3"]] +
        params$g_y[4]*spl[["b9_4"]] + params$t_y[1]*spl[["b12_1"]] +
        params$t_y[2]*spl[["b12_2"]] + params$t_y[3]*spl[["b12_3"]] +
        params$t_y[4]*spl[["b12_4"]]
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
}
