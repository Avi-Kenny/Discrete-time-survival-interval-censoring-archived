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
  
  # Set initial parameter estimates
  if (cfg$model_version==7) {
    par_init <- c(a_x=-6.2967, g_x1=-0.1535, g_x2=0.9796, t_x=0.5343, a_s=-2.3111, g_s1=-0.5649, g_s2=0.6198, t_s=0.4245, beta_x=1.401, a_y=-5.5786, g_y1=0.3278, g_y2=4.2046, t_y=-0.7198)
  } else if (cfg$model_version==29) {
    par_init <- c(a_x=-6.087, g_x1=1.283, g_x2=-0.983, g_x3=-1.032, g_x4=-0.888, t_x1=0.602, t_x2=-1.874, t_x3=-1.472, t_x4=-0.530, a_s=-3.209, g_s1=3.658, g_s2=2.351, g_s3=4.186, g_s4=0.942, beta_x1=1.981, beta_x2=-0.784, beta_x3=0.0149, beta_x4=-1.884, beta_x5=-2.042, beta_x6=-0.565, a_y=-8.153, g_y1=2.255, g_y2=1.771, g_y3=4.729, g_y4=2.989, t_y1=-0.174, t_y2=-0.103, t_y3=-0.288, t_y4=0.126)
  } else if (cfg$model_version==30) {
    par_init <- c(a_x=-6.087, g_x1=1.283, g_x2=-0.983, g_x3=-1.032, g_x4=-0.888, t_x1=0.602, t_x2=-1.874, t_x3=-1.472, t_x4=-0.530, a_s=-3.209, g_s1=3.658, g_s2=2.351, g_s3=4.186, g_s4=0.942, beta_x1=0, beta_x2=0, beta_x3=0, beta_x4=0, a_y=-8.153, g_y1=2.255, g_y2=1.771, g_y3=4.729, g_y4=2.989, t_y1=-0.174, t_y2=-0.103, t_y3=-0.288, t_y4=0.126)
  } else if (cfg$model_version==31) {
    par_init <- c(a_x=-5.564, g_x1=0.569, g_x2=-1.992, g_x3=1.610, g_x4=-1.229, t_x1=1.113, t_x2=-1.773, t_x3=-1.267, t_x4=-0.453, a_s=-3.198, g_s1=3.563, g_s2=2.150, g_s3=4.358, g_s4=0.510, beta_x1=1.632, beta_x2=-0.406, beta_x3=-0.511, beta_x4=-0.525, beta_x5=-1.337, beta_x6=-0.472, a_y=-8.075, g_y1=2.393, g_y2=1.799, g_y3=5.165, g_y4=2.703, t_y1=-0.393, t_y2=-0.460, t_y3=-0.785, t_y4=0.233)
  } else if (cfg$model_version==32) {
    par_init <- c(a_x=-5.564, g_x1=0.569, g_x2=-1.992, g_x3=1.610, g_x4=-1.229, t_x1=1.113, t_x2=-1.773, t_x3=-1.267, t_x4=-0.453, a_s=-3.198, g_s1=3.563, g_s2=2.150, g_s3=4.358, g_s4=0.510, beta_x1=0, beta_x2=0, beta_x3=0, beta_x4=0, beta_x5=0, beta_x6=0, beta_x7=0, beta_x8=0, beta_x9=0, a_y=-8.075, g_y1=2.393, g_y2=1.799, g_y3=5.165, g_y4=2.703, t_y1=-0.393, t_y2=-0.460, t_y3=-0.785, t_y4=0.233)
  } else if (cfg$model_version==33) {
    par_init <- c(a_x=-5.564, g_x1=0.569, g_x2=-1.992, g_x3=1.610, g_x4=-1.229, t_x1=1.113, t_x2=-1.773, t_x3=-1.267, t_x4=-0.453, a_s=-3.198, g_s1=3.563, g_s2=2.150, g_s3=4.358, g_s4=0.510, beta_x1=0, beta_x2=0, beta_x3=0, beta_x4=0, beta_x5=0, beta_x6=0, beta_x7=0, beta_x8=0, beta_x9=0, a_y=-8.075, g_y1=2.393, g_y2=1.799, g_y3=5.165, g_y4=2.703, t_y1=-0.393, t_y2=-0.460, t_y3=-0.785, t_y4=0.233)
  } else if (cfg$model_version==34) {
    par_init <- c(a_x=-5.564, g_x1=0.569, g_x2=-1.992, g_x3=1.610, g_x4=-1.229, t_x1=1.113, t_x2=-1.773, t_x3=-1.267, t_x4=-0.453, a_s=-3.198, g_s1=3.563, g_s2=2.150, g_s3=4.358, g_s4=0.510, beta_x1=0, beta_x2=0, beta_x3=0, beta_x4=0, a_y=-8.075, g_y1=2.393, g_y2=1.799, g_y3=5.165, g_y4=2.703, t_y1=-0.393, t_y2=-0.460, t_y3=-0.785, t_y4=0.233)
  } else if (cfg$model_version==35) {
    par_init <- c(a_x=-5.564, g_x1=0.569, g_x2=-1.992, g_x3=1.610, g_x4=-1.229, t_x1=0, a_s=-3.198, g_s1=3.563, g_s2=2.150, g_s3=4.358, g_s4=0.510, beta_x1=0, beta_x2=0, beta_x3=0, beta_x4=0, a_y=-8.075, g_y1=2.393, g_y2=1.799, g_y3=5.165, g_y4=2.703, t_y1=-0.393, t_y2=-0.460, t_y3=-0.785, t_y4=0.233)
  } else if (cfg$model_version==36) {
    # par_init <- c(a_x=-5.2, g_x1=0.6, g_x2=-1.8, g_x3=0.3, g_x4=-2.0, t_x1=0.8, t_x2=-1.1, t_x3=-0.8, t_x4=-0.6, a_s=-3.1, g_s1=3.5, g_s2=2.0, g_s3=4.3, g_s4=0.6, beta_x1=2.0, beta_x2=0, beta_x3=-0.8, beta_x4=-0.1, a_y=-7.9, g_y1=2.0, g_y2=1.8, g_y3=4.6, g_y4=2.8, t_y1=-0.1, t_y2=-0.4, t_y3=-0.5, t_y4=0.3)
    if (cfg2$model_sex=="Female") { par_init <- c(a_x=-5.6, g_x1=0.69, g_x2=-1.27, g_x3=1.04, g_x4=-2.38, t_x1=1.21, t_x2=-1.64, t_x3=-1.04, t_x4=-0.51, a_s=-3.18, g_s1=3.51, g_s2=2.13, g_s3=4.35, g_s4=0.57, beta_x1=2.11, beta_x2=-0.03, beta_x3=-0.39, beta_x4=-0.15, a_y=-8.02, g_y1=2.13, g_y2=1.66, g_y3=4.37, g_y4=2.8, t_y1=-0.17, t_y2=-0.17, t_y3=-0.37, t_y4=0.06) }
    if (cfg2$model_sex=="Male") { par_init <- c(a_x=-6.61, g_x1=2.41, g_x2=-1.89, g_x3=-0.42, g_x4=-1.73, t_x1=0.48, t_x2=-3.19, t_x3=-0.52, t_x4=-0.43, a_s=-3.75, g_s1=3.84, g_s2=3.03, g_s3=4.03, g_s4=2.09, beta_x1=2.22, beta_x2=-0.03, beta_x3=-1.15, beta_x4=-0.17, a_y=-7.18, g_y1=1.91, g_y2=1.9, g_y3=3.95, g_y4=2.71, t_y1=-0.27, t_y2=-0.03, t_y3=-0.12, t_y4=0.02) }
  } else if (cfg$model_version==37) {
    if (cfg2$model_sex=="Female") { par_init <- c(a_x=-5.6, g_x1=25, g_x2=-32, g_x3=-5, g_x4=7, t_x1=0.02, t_x2=0.2, t_x3=-0.8, t_x4=0.7, a_s=-3.18, g_s1=30, g_s2=-17, g_s3=-18, g_s4=0, beta_x1=2.11, beta_x2=-0.03, beta_x3=-0.39, beta_x4=-0.15, a_y=-8.02, g_y1=20, g_y2=-13, g_y3=-5, g_y4=4, t_y1=-0.07, t_y2=0.075, t_y3=-0.01, t_y4=0.04) }
    if (cfg2$model_sex=="Male") { par_init <- c(a_x=-6.61, g_x1=13, g_x2=-3, g_x3=-30, g_x4=13, t_x1=0.2, t_x2=-0.4, t_x3=-0.5, t_x4=1.3, a_s=-3.75, g_s1=15, g_s2=10, g_s3=-25, g_s4=-4, beta_x1=2.22, beta_x2=-0.03, beta_x3=-1.15, beta_x4=-0.17, a_y=-7.18, g_y1=15, g_y2=-7, g_y3=-4, g_y4=1, t_y1=-0.04, t_y2=0.01, t_y3=0.07, t_y4=-0.02) }
  } else if (cfg$model_version==38) {
    # if (cfg2$model_sex=="Female") { par_init <- c(a_x=-5.6, g_x1=10, g_x2=-28, g_x3=13, t_x1=0.02, t_x2=0.2, t_x3=-0.8, t_x4=0.7, a_s=-3.18, g_s1=30, g_s2=-17, g_s3=-18, g_s4=0, beta_x1=2.11, beta_x2=-0.03, beta_x3=-0.39, beta_x4=-0.15, a_y=-8.02, g_y1=20, g_y2=-13, g_y3=-5, g_y4=4, t_y1=-0.07, t_y2=0.075, t_y3=-0.01, t_y4=0.04) }
    # if (cfg2$model_sex=="Male") { par_init <- c(a_x=-6.61, g_x1=11, g_x2=-30, g_x3=12, t_x1=0.2, t_x2=-0.4, t_x3=-0.5, t_x4=1.3, a_s=-3.75, g_s1=15, g_s2=10, g_s3=-25, g_s4=-4, beta_x1=2.22, beta_x2=-0.03, beta_x3=-1.15, beta_x4=-0.17, a_y=-7.18, g_y1=15, g_y2=-7, g_y3=-4, g_y4=1, t_y1=-0.04, t_y2=0.01, t_y3=0.07, t_y4=-0.02) }
    if (cfg2$model_sex=="Female") { par_init <- c(a_x=-5.47, g_x1=10.01, g_x2=-27.66, g_x3=13.78, t_x1=0.01, t_x2=0.2, t_x3=-0.81, t_x4=0.7, a_s=-3.17, g_s1=30.03, g_s2=-16.89, g_s3=-17.81, g_s4=0.04, beta_x1=2.1, beta_x2=-0.03, beta_x3=-0.4, beta_x4=-0.15, a_y=-8.03, g_y1=19.99, g_y2=-12.99, g_y3=-4.94, g_y4=4.24, t_y1=-0.07, t_y2=0.08, t_y3=-0.02, t_y4=0.03) }
    if (cfg2$model_sex=="Male") { par_init <- c(a_x=-6.59, g_x1=11.58, g_x2=-30.05, g_x3=13.13, t_x1=0.2, t_x2=-0.4, t_x3=-0.5, t_x4=1.32, a_s=-3.75, g_s1=14.95, g_s2=9.95, g_s3=-24.96, g_s4=-3.94, beta_x1=2.22, beta_x2=-0.03, beta_x3=-1.18, beta_x4=-0.18, a_y=-7.18, g_y1=15, g_y2=-7.04, g_y3=-4.05, g_y4=1.23, t_y1=-0.04, t_y2=0.01, t_y3=0.07, t_y4=-0.01) }
  }
  
  # Construct param init vector
  if (F) {
    ests <- readRDS("objs/ests_38_full_F_20241209.rds")
    str <- "par_init <- c("
    for (par in names(ests$opt$par)) {
      str <- paste0(str, par, "=", round(ests$opt$par[[par]], 2), ", ")
    }
    str <- paste0(substr(str, 1, nchar(str)-2), ")")
    print(str)
  }
  
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
           substr(cfg2$model_sex,1,1), "_", format(Sys.time(), "%Y%m%d"),
           ".rds")
  )
  
  # if (cfg2$parallelize) { stopCluster(cl) }
  
  # Parse results
  res <- data.frame(
    "param" = character(),
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
      param = names(par_init)[i],
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
