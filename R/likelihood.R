#' Transform dataset
#' @param dat A dataset returned by generate_dataset()
#' @param model_version Model version, as specified in MAIN.R
transform_dataset <- function(dat, model_version=0) {
  
  # Construct spline bases
  b2 <- construct_basis("age (13,20,30,60,90)")
  b3 <- construct_basis("age (13,30,60,75,90)")
  b4 <- construct_basis("year (00,05,10,15,20)")
  
  n <- attr(dat, "n")
  
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
  
  # Construct a dataframe specific to each individual
  dat_objs <- lapply(c(1:n), function(i) {
    
    d <- list()
    # d$dat_i <- dat[dat$id==i,]
    d$dat_i <- dat[row_map[[as.character(i)]],]
    
    # Start and end times
    d$s_i <- attr(dat, "s_i")[i]
    d$t_i <- attr(dat, "t_i")[i]
    
    # Create new variables
    d$dat_i$init_visit <- 0
    d$dat_i[1,"init_visit"] <- 1
    d$dat_i$cal_time_sc <- d$dat_i$t_end / 10 # Rescaling originally done to prevent optimization issues; check if this is still needed
    
    # Apply spline bases to dataframe
    # !!!!! This section can be sped up by modifying construct_basis()
    if (model_version %in% c(13:18)) {
      d$dat_i$b2_1 <- signif(sapply(d$dat_i$w_2, function(w_2) { b2(w_2,1) }),4)
      d$dat_i$b2_2 <- signif(sapply(d$dat_i$w_2, function(w_2) { b2(w_2,2) }),4)
      d$dat_i$b2_3 <- signif(sapply(d$dat_i$w_2, function(w_2) { b2(w_2,3) }),4)
      d$dat_i$b2_4 <- signif(sapply(d$dat_i$w_2, function(w_2) { b2(w_2,4) }),4)
      d$dat_i$b3_1 <- signif(sapply(d$dat_i$w_2, function(w_2) { b3(w_2,1) }),4)
      d$dat_i$b3_2 <- signif(sapply(d$dat_i$w_2, function(w_2) { b3(w_2,2) }),4)
      d$dat_i$b3_3 <- signif(sapply(d$dat_i$w_2, function(w_2) { b3(w_2,3) }),4)
      d$dat_i$b3_4 <- signif(sapply(d$dat_i$w_2, function(w_2) { b3(w_2,4) }),4)
      d$dat_i$b4_1 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b4(j,1) }),4)
      d$dat_i$b4_2 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b4(j,2) }),4)
      d$dat_i$b4_3 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b4(j,3) }),4)
      d$dat_i$b4_4 <- signif(sapply(d$dat_i$cal_time_sc, function(j) { b4(j,4) }),4)
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
        d$X_i_set <- c(d$X_i_set, list(x_))
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
construct_negloglik <- function(parallelize=FALSE, model_version=0) {
  
  negloglik <- function(par) {
    
    counter <<- counter + 1
    if (counter%%10==0) {
      print(paste0(counter, " negloglik calls processed: ",Sys.time())) # !!!!!
    }
    # print(sort(sapply(ls(parent.frame(2)),function(x) { # !!!!!
    #   format(object.size(get(x)), units="MB") # !!!!!
    # }))) # !!!!!
    # browser() # !!!!!

    # Convert parameter vector to a named list
    p <- as.numeric(par)
    if (model_version==0) {
      params <- list(a_x=p[1], g_x=p[2:3], a_y=p[4], g_y=p[5:6], beta_x=p[7], beta_z=p[8], t_x=p[9], t_y=p[10], a_s=p[11], t_s=p[12], g_s=p[13:14])
    } else if (model_version==1) {
      params <- list(a_x=p[1], a_y=p[2], beta_x=p[3], beta_z=p[4], a_s=p[5])
    } else if (model_version==2) {
      params <- list(a_x=p[1], a_y=p[2], beta_x=p[3], beta_z=p[4], a_s=p[5], g_x1=p[6], g_y1=p[7], g_s1=p[8])
    } else if (model_version %in% c(3,4)) {
      params <- list(a_x=p[1], g_x=p[2:3], a_y=p[4], g_y=p[5:6], beta_x=p[7], beta_z=p[8], a_s=p[9], g_s=p[10:11])
    } else if (model_version==5) {
      params <- list(a_x=p[1], g_x=p[2:3], a_y=p[4], g_y=p[5:6], beta_x=p[7], beta_z=p[8], t_y=p[9], a_s=p[10], g_s=p[11:12])
    } else if (model_version==6) {
      params <- list(a_x=p[1], g_x=p[2:3], a_y=p[4], g_y=p[5:6], beta_x=p[7], beta_z=p[8], t_x=p[9], t_y=p[10], a_s=p[11], g_s=p[12:13])
    } else if (model_version==7) {
      params <- list(a_x=p[1], g_x=p[2:3], a_y=p[4], g_y=p[5:6], beta_x=p[7], beta_z=p[8], t_x=p[9], t_y=p[10], a_s=p[11], t_s=p[12], g_s=p[13:14])
    } else if (model_version==8) {
      params <- list(a_x=p[1], g_x=p[2:3], a_y=p[4], g_y=p[5:8], beta_x=p[9], beta_z=p[10], t_x=p[11], t_y=p[12], a_s=p[13], t_s=p[14], g_s=p[15:16])
    } else if (model_version==9) {
      params <- list(a_x=p[1], g_x=p[2:3], a_y=p[4], g_y=p[5:9], beta_x=p[10], beta_z=p[11], t_x=p[12], t_y=p[13], a_s=p[14], t_s=p[15], g_s=p[16:17])
    } else if (model_version==10) {
      params <- list(beta_z=p[1], a_y=p[2], g_y=p[3:7], t_y=p[8])
    } else if (model_version %in% c(11,12)) {
      params <- list(a_x=p[1], g_x=p[2:6], t_x=p[7], a_s=p[8], g_s=p[9:10], t_s=p[11], beta_x=p[12], beta_z=p[13], a_y=p[14], g_y=p[15:19], t_y=p[20])
    } else if (model_version==13) {
      params <- list(a_x=p[1], g_x=p[2:6], t_x=p[7:10], a_s=p[11], g_s=p[12:13], t_s=p[14], beta_x=p[15], beta_z=p[16], a_y=p[17], g_y=p[18:22], t_y=p[23])
    } else if (model_version==14) {
      params <- list(a_x=p[1], g_x=p[2:6], t_x=p[7:10], a_s=p[11], g_s=p[12:13], t_s=p[14], beta_x=p[15], beta_z=p[16], a_y=p[17], g_y=p[18:22], t_y=p[23:26])
    } else if (model_version==15) {
      params <- list(a_x=p[1], g_x=p[2:9], t_x=p[10:13], a_s=p[14], g_s=p[15:16], t_s=p[17], beta_x=p[18], beta_z=p[19], a_y=p[20], g_y=p[21:25], t_y=p[26:29])
    } else if (model_version==16) {
      params <- list(a_x=p[1], g_x=p[2:9], t_x=p[10:13], a_s=p[14], g_s=p[15:16], t_s=p[17:20], beta_x=p[21], beta_z=p[22], a_y=p[23], g_y=p[24:28], t_y=p[29:32])
    } else if (model_version==17) {
      params <- list(a_x=p[1], g_x=p[2:9], t_x=p[10:13], a_s=p[14], g_s=p[15:19], beta_x=p[20], beta_z=p[21], a_y=p[22], g_y=p[23:27], t_y=p[28:31])
    } else if (model_version==18) {
      params <- list(a_x=p[1], g_x=p[2:9], t_x=p[10:13], a_s=p[14], g_s=p[15:19], beta_x=p[20:21], beta_z=p[22:23], a_y=p[24], g_y=p[25:29], t_y=p[30:33])
    }
    
    # Compute the negative likelihood across individuals
    if (parallelize) {
      
      # parallel::clusterEvalQ(cl, sink(paste0("C:/Users/ak811/Desktop/Avi/Research/HIVMI/output", Sys.getpid(), ".txt"))) # !!!!!
      # parallel::parLapply(cl, c(1:5), function(i) {
      #   print("Check 1")
      #   print(pryr::mem_used())
      # })
      # print(pryr::object_size(dat_objs))
      # parallel::clusterExport(cl, c("dat_objs"), envir=.GlobalEnv)

      # return(-1 * sum(unlist(parallel::parLapply(
      #   cl,
      #   c(1:cfg$sim_n_cores),
      #   function(batch) {
      #     print("START BATCH")
      #     print("mem_used() START") # !!!!!
      #     print(pryr::mem_used()) # !!!!!
      #     inner_sum <- 0
      #     for (i in batches[[batch]]) {
      #       lik_i <- lik_fn(i, params, inds)
      #       print("mem_used() BEFORE") # !!!!!
      #       print(pryr::mem_used()) # !!!!!
      #       print("mem_used() AFTER") # !!!!!
      #       # browser() # !!!!!
      #       # mem_used() # !!!!!
      #       # print("Mem change:")
      #       # print(mem_change({ lik_i <- lik_fn(i) }))
      #       inner_sum <- inner_sum + log(lik_i)
      #     }
      #     return(inner_sum)
      #     # sum(log(unlist(lapply(batches[[batch]], lik_fn))))
      #   })
      # )))
      
      # # !!!!! Testing v1
      # if (Sys.getenv("avi_test")=="type1") {
      #   print("Type 1")
      #   nll <- -1 * sum(log(unlist(parallel::parLapply(cl, dat_objs, function(d) {
      #     lik_fn2(d, params, inds)
      #   }))))
      #   print(paste0("NLL: ", nll))
      #   print(pryr::mem_used())
      #   return(nll)
      # }
      
      # # !!!!! TEMPORARY; for testing
      # nll <- -1 * sum(log(unlist(
      #   lapply(dat_objs_wrapper, function(d) {
      #     lapply(d, function(d2) { lik_fn2(d2, params, inds) })
      #   })
      # )))
      
      
      
      # !!!!! Testing v2
      # if (Sys.getenv("avi_test")=="type2") {
        # print("Type 2")
        nll <- -1 * sum(log(unlist(
          parallel::parLapply(cl, dat_objs_wrapper, function(d) {
            lapply(d, function(d2) { lik_fn2(d2, params, inds) })
          })
        )))
        # print(paste0("NLL: ", nll))
        # print(pryr::mem_used())
        return(nll)
      # }
      
    } else {
      
      stop("This section needs to be updated.")
      
      # return(-1 * sum(log(unlist(lapply(c(1:n), lik_fn)))))
      
      inner_sum <- 0
      for (i in c(1:n)) {
        lik_i <- lik_fn(i, params, inds)
        # browser() # !!!!!
        print(mem_used()) # !!!!!
        # print("Mem change:")
        # print(mem_change({ lik_i <- lik_fn(i) }))
        inner_sum <- inner_sum + log(lik_i)
      }
      return(-1*inner_sum)
    }
    
  }
  
  return(negloglik)
  
}



#' Calculate likelihood for individual i
#'
#' @param i Index for an individual
#' @param params TO DO
#' @param inds TO DO
#' @return Numeric likelihood
#' @note dat_objs is accessed globally
lik_fn <- function(i, params, inds) {
  
  print(paste0("i: ", i))
  # Extract data for individual i
  d <- dat_objs[[i]]
  
  # Compute the likelihood for individual i
  f2 <- sum(unlist(lapply(d$X_i_set, function(x) {
    d$dat_i$x <- x
    if (length(x)==1) {
      d$dat_i$x_prev <- 0
    } else {
      d$dat_i$x_prev <- c(0,x[1:(length(x)-1)]) # !!!!! This can be precomputed
    }
    return(prod(apply(X=d$dat_i, MARGIN=1, FUN = function(r) {
      x <- r[["x"]]
      w_ij <- as.numeric(r[inds$w])
      j <- r[["cal_time_sc"]]
      spl_ij <- r[inds$spl]
      return(
        f_x(x=x, x_prev=r[["x_prev"]], w=w_ij, j=j,
            s=r[["init_visit"]], spl=spl_ij, params=params) *
          f_y(y=r[["y"]], x=x, w=w_ij, z=r[["z"]], j=j,
              spl=spl_ij, params=params)
      )
    })))
  })))
  
  browser() # !!!!!

  if (f2<=0) {
    f2 <- 1e-10
    # warning("Likelihood of zero")
  }
  
  return(f2)
  
}



#' Calculate likelihood for individual i
#'
#' @param d One component of the list dat_objs
#' @param params TO DO
#' @param inds TO DO
#' @return Numeric likelihood
#' @note dat_objs is accessed globally
lik_fn2 <- function(d, params, inds) {
  
  # print(paste0("i: ", i))
  # Extract data for individual i
  # d <- dat_objs[[i]]
  
  # !!!!! New code
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
    
    f2 <- sum(unlist(lapply(d$X_i_set, function(x) {
      if (length(x)==1) { # !!!!! This can be precomputed
        x_prev <- 0
      } else {
        x_prev <- c(0,x[1:(length(x)-1)])
      }
      
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
      
    })))
    
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
if (cfg$model_version==0) {
  
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
  
} else if (cfg$model_version==1) {
  
  f_x <- function(x, x_prev, w, j, s, spl, params) {
    if (s==0) {
      if (x==1 && x_prev==1) {
        return(1)
      } else {
        prob <- icll(params$a_x)
        if (x==1) { return(prob) } else { return(1-prob) }
      }
    } else {
      prob <- icll(params$a_s)
      if (x==1) { return(prob) } else { return(1-prob) }
    }
  }
  
} else if (cfg$model_version==2) {
  
  f_x <- function(x, x_prev, w, j, s, spl, params) {
    if (s==0) {
      if (x==1 && x_prev==1) {
        return(1)
      } else {
        prob <- icll(params$a_x + params$g_x1*w[1])
        if (x==1) { return(prob) } else { return(1-prob) }
      }
    } else {
      prob <- icll(params$a_s + params$g_s1*w[1])
      if (x==1) { return(prob) } else { return(1-prob) }
    }
  }
  
} else if (cfg$model_version %in% c(3,4,5)) {
  
  f_x <- function(x, x_prev, w, j, s, spl, params) {
    if (s==0) {
      if (x==1 && x_prev==1) {
        return(1)
      } else {
        prob <- icll(params$a_x + sum(params$g_x*w))
        if (x==1) { return(prob) } else { return(1-prob) }
      }
    } else {
      prob <- icll(params$a_s + sum(params$g_s*w))
      if (x==1) { return(prob) } else { return(1-prob) }
    }
  }
  
} else if (cfg$model_version==6) {
  
  f_x <- function(x, x_prev, w, j, s, spl, params) {
    if (s==0) {
      if (x==1 && x_prev==1) {
        return(1)
      } else {
        prob <- icll(params$a_x + params$t_x*j + sum(params$g_x*w))
        if (x==1) { return(prob) } else { return(1-prob) }
      }
    } else {
      prob <- icll(params$a_s + sum(params$g_s*w))
      if (x==1) { return(prob) } else { return(1-prob) }
    }
  }
  
} else if (cfg$model_version %in% c(7,8,9)) {
  
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
  
} else if (cfg$model_version==10) {
  
  f_x <- function(x, x_prev, w, j, s, spl, params) { 1 }
  
} else if (cfg$model_version==11) {
  
  b1 <- construct_basis("age (0-100), 4DF")
  f_x <- function(x, x_prev, w, j, s, spl, params) {
    if (s==0) {
      if (x==1 && x_prev==1) {
        return(1)
      } else {
        prob <- icll(
          params$a_x + params$t_x*j + params$g_x[1]*w[1] +
            params$g_x[2]*b1(w[2],1) + params$g_x[3]*b1(w[2],2) +
            params$g_x[4]*b1(w[2],3) + params$g_x[5]*b1(w[2],4)
        )
        if (x==1) { return(prob) } else { return(1-prob) }
      }
    } else {
      prob <- icll(params$a_s + params$t_s*j + sum(params$g_s*w))
      if (x==1) { return(prob) } else { return(1-prob) }
    }
  }
  
} else if (cfg$model_version==12) {
  
  f_x <- (function() {
    b2 <- construct_basis("age (13,20,30,60,90)")
    f_x <- function(x, x_prev, w, j, s, spl, params) {
      if (s==0) {
        if (x==1 && x_prev==1) {
          return(1)
        } else {
          prob <- icll(
            params$a_x + params$t_x*j + params$g_x[1]*w[1] +
              params$g_x[2]*b2(w[2],1) + params$g_x[3]*b2(w[2],2) +
              params$g_x[4]*b2(w[2],3) + params$g_x[5]*b2(w[2],4)
          )
          if (x==1) { return(prob) } else { return(1-prob) }
        }
      } else {
        prob <- icll(params$a_s + params$t_s*j + sum(params$g_s*w))
        if (x==1) { return(prob) } else { return(1-prob) }
      }
    }
    return(f_x)
  })()
  
} else if (cfg$model_version %in% c(13,14)) {
  
  f_x <- function(x, x_prev, w, j, s, spl, params) {
    if (s==0) {
      if (x==1 && x_prev==1) {
        return(1)
      } else {
        prob <- icll(
          params$a_x + params$t_x[1]*spl[["b4_1"]] +
            params$t_x[2]*spl[["b4_2"]] +
            params$t_x[3]*spl[["b4_3"]] + params$t_x[4]*spl[["b4_4"]] +
            params$g_x[1]*w[1] + params$g_x[2]*spl[["b2_1"]] +
            params$g_x[3]*spl[["b2_2"]] + params$g_x[4]*spl[["b2_3"]] +
            params$g_x[5]*spl[["b2_4"]]
        )
        if (x==1) { return(prob) } else { return(1-prob) }
      }
    } else {
      prob <- icll(params$a_s + params$t_s*j + sum(params$g_s*w))
      if (x==1) { return(prob) } else { return(1-prob) }
    }
  }
  
} else if (cfg$model_version==15) {
  
  f_x <- function(x, x_prev, w, j, s, spl, params) {
    if (s==0) {
      if (x==1 && x_prev==1) {
        return(1)
      } else {
        prob <- icll(
          params$a_x + params$t_x[1]*spl[["b4_1"]] +
            params$t_x[2]*spl[["b4_2"]] + params$t_x[3]*spl[["b4_3"]] +
            params$t_x[4]*spl[["b4_4"]] +
            w[1]*(
              params$g_x[1]*spl[["b2_1"]] + params$g_x[2]*spl[["b2_2"]] +
                params$g_x[3]*spl[["b2_3"]] + params$g_x[4]*spl[["b2_4"]]
            ) +
            (1-w[1])*(
              params$g_x[5]*spl[["b2_1"]] + params$g_x[6]*spl[["b2_2"]] +
                params$g_x[7]*spl[["b2_3"]] + params$g_x[8]*spl[["b2_4"]]
            )
        )
        if (x==1) { return(prob) } else { return(1-prob) }
      }
    } else {
      prob <- icll(params$a_s + params$t_s*j + sum(params$g_s*w))
      if (x==1) { return(prob) } else { return(1-prob) }
    }
  }
  
} else if (cfg$model_version==16) {
  
  f_x <- function(x, x_prev, w, j, s, spl, params) {
    if (s==0) {
      if (x==1 && x_prev==1) {
        return(1)
      } else {
        prob <- icll(
          params$a_x + params$t_x[1]*spl[["b4_1"]] +
            params$t_x[2]*spl[["b4_2"]] + params$t_x[3]*spl[["b4_3"]] +
            params$t_x[4]*spl[["b4_4"]] +
            w[1]*(
              params$g_x[1]*spl[["b2_1"]] + params$g_x[2]*spl[["b2_2"]] +
                params$g_x[3]*spl[["b2_3"]] + params$g_x[4]*spl[["b2_4"]]
            ) +
            (1-w[1])*(
              params$g_x[5]*spl[["b2_1"]] + params$g_x[6]*spl[["b2_2"]] +
                params$g_x[7]*spl[["b2_3"]] + params$g_x[8]*spl[["b2_4"]]
            )
        )
        if (x==1) { return(prob) } else { return(1-prob) }
      }
    } else {
      prob <- icll(
        params$a_s + params$t_s[1]*spl[["b4_1"]] + params$t_s[2]*spl[["b4_2"]] +
          params$t_s[3]*spl[["b4_3"]] + params$t_s[4]*spl[["b4_4"]] +
          sum(params$g_s*w)
      )
      if (x==1) { return(prob) } else { return(1-prob) }
    }
  }
  
} else if (cfg$model_version %in% c(17,18)) {
  
  f_x <- function(x, x_prev, w, j, s, spl, params) {
    if (s==0) {
      if (x==1 && x_prev==1) {
        return(1)
      } else {
        prob <- icll(
          params$a_x + params$t_x[1]*spl[["b4_1"]] +
            params$t_x[2]*spl[["b4_2"]] + params$t_x[3]*spl[["b4_3"]] +
            params$t_x[4]*spl[["b4_4"]] +
            w[1]*(
              params$g_x[1]*spl[["b2_1"]] + params$g_x[2]*spl[["b2_2"]] +
                params$g_x[3]*spl[["b2_3"]] + params$g_x[4]*spl[["b2_4"]]
            ) +
            (1-w[1])*(
              params$g_x[5]*spl[["b2_1"]] + params$g_x[6]*spl[["b2_2"]] +
                params$g_x[7]*spl[["b2_3"]] + params$g_x[8]*spl[["b2_4"]]
            )
        )
        if (x==1) { return(prob) } else { return(1-prob) }
      }
    } else {
      prob <- icll(
        params$a_s + params$g_s[1]*w[1] + params$g_s[2]*spl[["b3_1"]] +
          params$g_s[3]*spl[["b3_2"]] + params$g_s[4]*spl[["b3_3"]] + 
          params$g_s[5]*spl[["b3_4"]]
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
if (cfg$model_version==0) {
  
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      params$a_y + params$t_y*j + sum(params$g_y*w) + params$beta_x*x*(1-z) +
        params$beta_z*x*z
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
} else if (cfg$model_version==1) {
  
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(params$a_y + params$beta_x*x*(1-z) + params$beta_z*x*z)
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
} else if (cfg$model_version==2) {
  
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      params$a_y + params$g_y1*w[1] + params$beta_x*x*(1-z) + params$beta_z*x*z
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
} else if (cfg$model_version==3) {
  
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      params$a_y + sum(params$g_y*w) + params$beta_x*x*(1-z) + params$beta_z*x*z
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
} else if (cfg$model_version==4) {
  
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      params$a_y + sum(params$g_y*w) + params$beta_x*x*(1-z) + params$beta_z*x*z
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
} else if (cfg$model_version %in% c(5,6,7)) {
  
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      params$a_y + params$t_y*j + sum(params$g_y*w) + params$beta_x*x*(1-z) +
        params$beta_z*x*z
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
} else if (cfg$model_version==8) {
  
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      params$a_y + params$t_y*j + params$g_y[1]*w[1] + params$g_y[2]*w[2] +
        params$g_y[3]*(w[2]^2) + params$g_y[4]*(w[2]^3) +
        params$beta_x*x*(1-z) + params$beta_z*x*z
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
} else if (cfg$model_version %in% c(9,11)) {
  
  b1 <- construct_basis("age (0-100), 4DF")
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      params$a_y + params$t_y*j + params$g_y[1]*w[1] +
        params$g_y[2]*b1(w[2],1) + params$g_y[3]*b1(w[2],2) +
        params$g_y[4]*b1(w[2],3) + params$g_y[5]*b1(w[2],4) +
        params$beta_x*x*(1-z) + params$beta_z*x*z
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
} else if (cfg$model_version==10) {
  
  # Note: there was a bug in this model when it was first run (beta_x excluded)
  b1 <- construct_basis("age (0-100), 4DF")
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      params$beta_x*x*(1-z) + params$beta_z*x*z + params$a_y +
        params$g_y[1]*w[1] + params$g_y[2]*b1(w[2],1) +
        params$g_y[3]*b1(w[2],2) + params$g_y[4]*b1(w[2],3) +
        params$g_y[5]*b1(w[2],4) + params$t_y*j
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }

} else if (cfg$model_version %in% c(12,13)) {
  
  # Note: there was a bug in this model when it was first run (beta_x excluded)
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      params$beta_x*x*(1-z) + params$beta_z*x*z + params$a_y +
        params$g_y[1]*w[1] + params$g_y[2]*spl[["b3_1"]] +
        params$g_y[3]*spl[["b3_2"]] + params$g_y[4]*spl[["b3_3"]] +
        params$g_y[5]*spl[["b3_4"]] + params$t_y*j
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
} else if (cfg$model_version %in% c(14:17)) {
  
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      params$beta_x*x*(1-z) + params$beta_z*x*z + params$a_y +
        params$g_y[1]*w[1] + params$g_y[2]*spl[["b3_1"]] +
        params$g_y[3]*spl[["b3_2"]] + params$g_y[4]*spl[["b3_3"]] +
        params$g_y[5]*spl[["b3_4"]] + params$t_y[1]*spl[["b4_1"]] +
        params$t_y[2]*spl[["b4_2"]] + params$t_y[3]*spl[["b4_3"]] +
        params$t_y[4]*spl[["b4_4"]]
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
} else if (cfg$model_version==18) {
  
  f_y <- function(y, x, w, z, j, spl, params) {
    prob <- icll(
      (params$beta_x[1]+params$beta_x[2]*j)*x*(1-z) +
        (params$beta_z[1]+params$beta_z[2]*j)*x*z + params$a_y +
        params$g_y[1]*w[1] + params$g_y[2]*spl[["b3_1"]] +
        params$g_y[3]*spl[["b3_2"]] + params$g_y[4]*spl[["b3_3"]] +
        params$g_y[5]*spl[["b3_4"]] + params$t_y[1]*spl[["b4_1"]] +
        params$t_y[2]*spl[["b4_2"]] + params$t_y[3]*spl[["b4_3"]] +
        params$t_y[4]*spl[["b4_4"]]
    )
    if (y==1) { return(prob) } else { return(1-prob) }
  }
  
}
