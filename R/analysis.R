##################.
##### Config #####
##################.

cfg2$opt_maxit <- 5000
cfg2$opt_r <- 2 # 2 for speed, 4 for accuracy
cfg2$opt_reltol <- 1e-5



###############################.
##### Data quality checks #####
###############################.

if (F) {
  
  chk(2, "DQA: START")
  
  # Setup
  # !!!!! Needs to change; first_hiv_pos_dt and last_hiv_neg_dt vars removed
  dqa <- function(test) { if (test==F) { stop("DQA Error") } }
  dat_grp2 <- dat_raw %>% dplyr::group_by(id) %>%
    dplyr::summarize(
      unique_T_plus = (function(vec){ length(unique(vec)) })(first_hiv_pos_dt), # !!!!!
      unique_T_minus = (function(vec){ length(unique(vec)) })(last_hiv_neg_dt) # !!!!!
    )
  
  # Tests on raw data
  dqa(identical(sort(unique(dat_raw$sex)), c("Female", "Male")))
  dqa(sum(dat_raw$died==0, na.rm=T) + sum(dat_raw$died==1, na.rm=T)
      == length(unique(dat_raw$id)))
  dqa(max(dat_grp2$unique_T_plus)==1)
  dqa(max(dat_grp2$unique_T_minus)==1)
  
  # Tests on processed data
  dqa(length(unique(dat_prc$id))==max(dat_prc$id))
  dqa(sum(dat_grp$T_plus-dat_grp$T_minus<0, na.rm=T)==0)
  dqa(sum(dat_grp$case==999)==0)
  
  # Tests comparing dat_raw vs. dat_prc
  dqa(sum(dat_prc$y==1, na.rm=T)==sum(dat_raw$died==1, na.rm=T))
  dqa(sum(dat_prc$y==0, na.rm=T)==
        sum(dat_raw$died==0, na.rm=T) + sum(is.na(dat_raw$died)))
  
  
  
  # Distribution of "year at entry"
  dat_prc %>%
    group_by(id) %>%
    mutate(min_year=min(year)) %>%
    filter(year==min_year) %>%
    xtabs(~year, data=.)
  
  # Distribution of "year of test"
  xtabs(~ResultDate, data=dat_prc)
  
  # Distribution of "year of first test"
  dat_prc %>%
    group_by(id) %>%
    filter(!is.na(ResultDate)) %>%
    mutate(min_year=min(year)) %>%
    filter(year==min_year) %>%
    xtabs(~year, data=.)
  
  # Distribution of "year of POS test"
  dat_prc %>%
    group_by(id) %>%
    filter(!is.na(ResultDate) & HIVResult=="P") %>%
    xtabs(~year, data=.)
  
  # Distribution of "year of first POS test"
  dat_prc %>%
    group_by(id) %>%
    filter(!is.na(ResultDate) & HIVResult=="P") %>%
    mutate(min_year=min(year)) %>%
    filter(year==min_year) %>%
    xtabs(~year, data=.)
  
  # Out of "year at entry" years, which years had HIV tests?
  dat_prc %>%
    group_by(id) %>%
    mutate(min_year=min(year)) %>%
    filter(year==min_year) %>%
    filter(!is.na(ResultDate)) %>%
    xtabs(~year, data=.)
  
  
  # !!!!! New DQA: checking "initial sero" model
  dat_2 <- dat %>%
    group_by(id) %>%
    mutate(min_time=min(t_end), max_time=max(t_end))
  dat_3 <- dplyr::filter(dat_2, t_end==min_time)
  dat_4 <- dplyr::filter(dat_3, v==1)
  dat_4 %<>% dplyr::mutate(
    age_below_50 = In(age<=50),
    year_00_07 = In(t_end<=7),
    year_08_16 = In(t_end>8 & t_end<=16),
    year_17_22 =  In(t_end>16)
  )
  
  # % Positive (out of testers)
  mean(dplyr::filter(dat_4, T)$u) # Overall 34%
  mean(dplyr::filter(dat_4, w_1==1)$u) # Male: 24%
  mean(dplyr::filter(dat_4, w_1==0)$u) # Female: 38%
  length(dplyr::filter(dat_4, year_00_07==1)$u) # Year 00-07: xx%
  length(dplyr::filter(dat_4, year_08_16==1)$u) # Year 08-16: xx%
  length(dplyr::filter(dat_4, year_17_22==1)$u) # Year 17-22: xx%
  
  nrow(dplyr::filter(dat_4, w_1==0 & age_below_50==1 & year_00_07==1))
  nrow(dplyr::filter(dat_4, w_1==1 & age_below_50==1 & year_00_07==1))
  nrow(dplyr::filter(dat_4, w_1==0 & age_below_50==1 & year_08_16==1))
  nrow(dplyr::filter(dat_4, w_1==1 & age_below_50==1 & year_08_16==1))
  nrow(dplyr::filter(dat_4, w_1==0 & age_below_50==1 & year_17_22==1))
  nrow(dplyr::filter(dat_4, w_1==1 & age_below_50==1 & year_17_22==1))
  
  chk(2, "DQA: END")
  
}



#########################.
##### Data analysis #####
#########################.

if (cfg2$run_analysis) {
  
  dat_i_names <- names(dat_objs[[1]]$dat_i)
  inds <- list(
    w = which(dat_i_names %in% c("w_1", "w_2")),
    spl = which(
      substr(dat_i_names, 1, 1)=="b" & (
        substr(dat_i_names, 3, 3)=="_" | substr(dat_i_names, 4, 4)=="_"
      )
    )
  )
  
  # Construct log likelihood function
  chk(3, "construct_negloglik: START")
  n <- attr(dat, "n")
  print(paste0("Using ", cfg$sim_n_cores, " cores."))
  if (cfg2$parallelize) {
    n_batches <- cfg$sim_n_cores
  } else {
    n_batches <- 2
  }
  folds <- cut(c(1:n), breaks=n_batches, labels=FALSE)
  batches <- lapply(c(1:n_batches), function(batch) { c(1:n)[folds==batch] })
  dat_objs_wrapper <- lapply(batches, function(i) { dat_objs[i] })
  
  if (cfg2$parallelize) {
    cl <- parallel::makeCluster(cfg$sim_n_cores)
    objs_to_export <- c("f_x", "f_y", "icll", "lik_fn2", "inds", "batches",
                        "uncompress")
    parallel::clusterExport(cl, objs_to_export, envir=.GlobalEnv)
    negloglik <- construct_negloglik(
      inds=inds, parallelize=T, model_version=cfg$model_version, use_counter=T
    )
  } else {
    negloglik <- construct_negloglik(
      inds=inds, parallelize=F, model_version=cfg$model_version, use_counter=T
    )
  }
  chk(3, "construct_negloglik: END")
  
  # Run optimizer
  chk(4, "optim: START")
  counter <- 0
  st <- system.time({ nll_init <- negloglik(par_init) })
  print(st)
  print(paste("negloglik(par_init):", nll_init))
  opt <- stats::optim(
    par = par_init,
    fn = negloglik,
    method = "Nelder-Mead",
    control = list(maxit=cfg2$opt_maxit, reltol=cfg2$opt_reltol)
  )
  print("optim() finished.")
  print(opt)
  
  if (F) {
    stats::optim(par=par_init, fn=negloglik, control=list(trace=6))
    library(optimParallel)
    opt <- stats::optim(par=par_init, fn=negloglik)
  }
  chk(4, "optim: END")
  
  # Compute Hessian
  chk(5, "hessian: START")
  hessian_est <- numDeriv::hessian(
    func = negloglik,
    x = opt$par,
    method = "Richardson", # "Richardson" "complex"
    method.args = list(
      eps = 0.0001,
      d = 0.0001, # d gives the fraction of x to use for the initial numerical approximation
      zero.tol = sqrt(.Machine$double.eps/7e-7),
      r = cfg2$opt_r, # r gives the number of Richardson improvement iterations (repetitions with successly smaller d
      v = 2 # v gives the reduction factor
    )
  )
  
  tryCatch(
    expr = {
      hessian_inv <- solve(hessian_est)
    },
    error = function(e) {
      print("Hessian not invertible.")
      print(hessian_est)
    }
  )
  
  chk(5, "hessian: END")
  parallel::stopCluster(cl)
  saveRDS(
    list(opt=opt, hessian_inv=hessian_inv),
    paste0("ests_", cfg$model_version,
           ifelse(cfg2$samp_size==0, "_full_", "_partial_"),
           substr(cfg$model_sex,1,1), "_", format(Sys.time(), "%Y%m%d"),
           ".rds")
  )
  
  # if (cfg2$parallelize) { stopCluster(cl) }
  
  # Parse results
  res <- data.frame(
    "par" = character(),
    "par_init" = double(),
    # "par_true" = double(),
    "est" = double(),
    "se" = double(),
    "ci_lo" = double(),
    "ci_up" = double()
  )
  for (i in c(1:length(par_init))) {
    est <- as.numeric(opt$par[i])
    se <- sqrt(diag(hessian_inv))[i]
    res[i,] <- c(
      par = names(par_init)[i],
      par_init = round(par_init[i],4),
      # par_true = round(par_true[i],4),
      est = round(est,4),
      se = round(se,4),
      ci_lo = round(est-1.96*se,4), # Maybe do CI transformation later
      ci_up = round(est + 1.96*se,4) # Maybe do CI transformation later
    )
  }
  print(res)

  # !!!!! temp
  .t_end <- Sys.time()
  print("Total Runtime:")
  print(.t_end-.t_start)
  # saveRDS(
  #   list(
  #     runtime = .t_end - .t_start,
  #     n_cores = n_cores,
  #     opt_maxit = cfg2$opt_maxit,
  #     opt_reltol = cfg2$opt_reltol,
  #     samp_size = cfg2$samp_size,
  #     opt_r = cfg2$opt_r,
  #     res = res
  #   ),
  #   file = paste0("res_", runif(1), ".rds")
  # )
  
}

chk(6, "END")
