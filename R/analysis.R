##################.
##### Config #####
##################.

for (pkg in c(cfg$pkgs,cfg$pkgs_nocluster)) {
  suppressMessages({ do.call("library", list(pkg)) })
}

# !!!!! Testing: slurm_apply
if (F) {
  
  sim.pi <- function(iterations=1000) {
    print("hey 1")
    x.pos <- runif(iterations, min=-1, max=1)
    y.pos <- runif(iterations, min=-1, max=1)
    draw.pos <- ifelse(x.pos^2 + y.pos^2 <= 1, TRUE, FALSE)
    draws.in <- length(which(draw.pos == TRUE))
    print("hey 2")
    return(data.frame(iterations,draws.in))
  }
  
  params <- data.frame(iterations = rep(1000,100))
  
  sjob1 <- slurm_apply(
    f = sim.pi,
    params = params,
    jobname = "rslurm-pi-example",
    nodes = 10,
    cpus_per_node = 1
  )
  
  if (T) {
    slr_job <- sjob1
    res_files <- paste0("results_", 0:(slr_job$nodes - 1), ".RDS")
    print("res_files")
    print(res_files)
    tmpdir <- paste0("_rslurm_", slr_job$jobname)
    print("tmpdir")
    print(tmpdir)
  }
  
  res <- get_slurm_out(sjob1, outtype="table")
  
  my_pi <- 4/(sum(res$iterations)/sum(res$draws.in))
  cat("\n... done\n")
  cat(paste0(
    "pi estimated to ", my_pi, " over ", sum(res$iterations), " iterations\n"
  ))
  
  cleanup_files(sjob1)
  
}

# !!!!! Testing: slurm_map
if (F) {
  
  sim.pi <- function(iterations=1000) {
    x.pos <- runif(iterations, min=-1, max=1)
    y.pos <- runif(iterations, min=-1, max=1)
    draw.pos <- ifelse(x.pos^2 + y.pos^2 <= 1, TRUE, FALSE)
    draws.in <- length(which(draw.pos == TRUE))
    return(data.frame(iterations,draws.in))
  }
  
  params <- c(1:100)
  
  sjob1 <- slurm_map(
    x = params,
    f = sim.pi,
    jobname = "rslurm-pi-example",
    nodes = 10,
    cpus_per_node = 1
  )
  
  res <- get_slurm_out(sjob1, outtype="table")
  
  my_pi <- 4/(sum(res$iterations)/sum(res$draws.in))
  cat("\n... done\n")
  cat(paste0(
    "pi estimated to ", my_pi, " over ", sum(res$iterations), " iterations\n"
  ))
  
  cleanup_files(sjob1)
  
}

# stop("Stopping examples.")

chk(0, "START")
.t_start <- Sys.time()
cfg2 <- list(
  process_data = T,
  save_data = F,
  run_dqa = F,
  run_analysis = T,
  parallelize = T,
  use_simulated_dataset = F
)

# !!!!! TEMPORARY: testing different optim config options
if (T) {
  avi_maxit <- as.integer(Sys.getenv("avi_maxit"))
  avi_reltol <- as.numeric(Sys.getenv("avi_reltol"))
  avi_r <- as.numeric(Sys.getenv("avi_r"))
  # avi_maxit <- 200
  # avi_reltol <- 1e-5
  # avi_r <- 2
  print("CONFIG")
  print("------------")
  print(paste("maxit:", avi_maxit))
  print(paste("reltol:", avi_reltol))
  print(paste("r:", avi_r))
  print(paste("sample size:", as.integer(Sys.getenv("avi_samp_size"))))
  print("------------")
}



###########################.
##### Data processing #####
###########################.

chk(1, "Data reading/processing: START")
if (cfg2$use_simulated_dataset) {
  
  # Generate dataset
  # n <- 4000
  n <- 1000
  par_true_full <- lapply(list(
    a_x=0.005, a_y=0.003, a_v=0.7, a_z=0.01, g_x=c(1.3,1.2), g_y=c(1.2,1.1),
    g_v=c(1.2,1.1), g_z=c(1.2,1.1), beta_x=1.5, beta_z=0.7, a_s=0.05,
    g_s=c(2,1.5), t_s=1, t_x=1, t_y=1
  ), log)
  dat <- generate_data(n=n, max_time=70, params=par_true_full)
  print(paste("n:",n)) # !!!!!
  print(paste("rows in dataset:", nrow(dat))) # !!!!!
  
  # Manipulate dataset
  dat %<>% arrange(id, t_end)
  dat$x <- NULL # !!!!! New code: make this null and comment out below
  attr(dat, "T_plus") <- NULL
  attr(dat, "T_minus") <- NULL
  attr(dat, "case") <- NULL
  attr(dat, "max_time") <- NULL
  attr(dat, "params") <- NULL
  
  # if (cfg2$parallelize) { saveRDS(dat, "dat_sim_10000.rds") }
  
} else {
  
  if (cfg2$process_data) {
    
    # Read in data; see generate_data() for expected columns
    # dat_prc <- dat_raw <- read.csv("../../Data/data_raw_1000.csv")
    # dat_prc <- dat_raw <- read.csv("../Data/data_raw_1000_v2.csv")
    # dat_prc <- dat_raw <- read.csv("../Data/data_raw_1000_v3.csv")
    # dat_prc <- dat_raw <- read.csv("../Data/data_raw_10000.csv")
    
    set.seed(1)
    dat_prc <- dat_raw <- read.csv("../Data/data_raw_full.csv")
    iintids <- unique(dat_prc$iintid)
    samp_size <- 20000
    # samp_size <- as.integer(Sys.getenv("avi_samp_size")) # !!!!!
    iintids_sample <- sample(iintids, size=samp_size)
    dat_prc %<>% dplyr::filter(iintid %in% iintids_sample)
    
    # Drop unnecessary columns
    # dat_prc %<>% subset(select=-c(enter, exit, X_t, X_t0))
    dat_prc %<>% subset(select=-c(enter, exit))
    # !!!!! Drop hiv_result_fill, earliestartinitdate ?????
    
    # Start time
    window_start <- min(dat_prc$year)
    
    # Function to convert dates
    date_start <- as.Date("01jan1960", format="%d%b%Y")
    convert_date <- memoise(function(date) {
      date <- as.Date(date, format="%d%b%Y")
      # return(lubridate::interval(start=date_start, end=date) %/% years(1))
      return(lubridate::year(date))
    })
    
    # Data wrangling
    dat_prc %<>% dplyr::mutate(
      id = iintid,
      sex = ifelse(sex=="Male", 1, 0),
      died = ifelse(is.na(died), 0, died),
      dob = convert_date(dob),
      dod = convert_date(dod),
      first_hiv_pos_dt = convert_date(first_hiv_pos_dt),
      last_hiv_neg_dt = convert_date(last_hiv_neg_dt),
      earliestartinitdate = convert_date(earliestartinitdate),
      onart = ifelse(is.na(onart), 0, onart),
      year_prev = year - 1
    )
    
    # Filter out children with tests before age 13
    nrow(dat_prc)
    children_with_tests <- dplyr::filter(
      dat_prc, age_end<=13 & resultdate!=""
    )$iintid
    dat_prc %<>% dplyr::filter(!(iintid %in% children_with_tests))
    nrow(dat_prc)
    
    # Remove all data before 13th birthday
    nrow(dat_prc)
    dat_prc %<>% dplyr::filter(age_end>=13)
    nrow(dat_prc)
    
    # Filter out adults with tests after age 90
    # !!!!! Temporary
    nrow(dat_prc)
    adults90_with_tests <- dplyr::filter(
      dat_prc, age_end>90 & resultdate!=""
    )$iintid
    if (length(adults90_with_tests)>0) {
      dat_prc %<>% dplyr::filter(!(iintid %in% adults90_with_tests))
    }
    nrow(dat_prc)
    
    # Remove all data after 90th birthday
    nrow(dat_prc)
    dat_prc %<>% dplyr::filter(age_end<=90)
    nrow(dat_prc)
    
    # Filter out records with a negative test after a positive test
    nrow(dat_prc)
    dat_prc %<>% dplyr::filter(
      is.na(first_hiv_pos_dt) | is.na(last_hiv_neg_dt) |
        first_hiv_pos_dt>last_hiv_neg_dt
    )
    nrow(dat_prc)
    
    # Rearrange columns
    dat_prc %<>% dplyr::relocate(year_prev, .before=year)
    
    # Rename columns
    dat_prc %<>% dplyr::rename(
      "w_1" = sex,
      "y" = died,
      "z" = onart,
      "t_start" = year_prev,
      "t_end" = year
    )
    
    # Sort dataframe
    dat_prc %<>% dplyr::arrange(id,t_start)
    
    # Create grouped dataset
    dat_grp <- dat_prc %>% dplyr::group_by(id) %>%
      dplyr::summarize(
        count = n(),
        T_minus = last_hiv_neg_dt[1],
        T_plus = first_hiv_pos_dt[1],
        s_i = min(t_end),
        t_i = max(t_end)
      )
    dat_grp %<>% dplyr::mutate(
      case = case_when(
        is.na(T_minus) & is.na(T_plus) ~ 1,
        !is.na(T_minus) & is.na(T_plus) ~ 2,
        !is.na(T_minus) & !is.na(T_plus) ~ 3,
        is.na(T_minus) & !is.na(T_plus) ~ 4,
        TRUE ~ 999
      )
    )
    
    # !!!!! Check that `first_hiv_pos_dt` and `last_hiv_neg_dt` lie within `t_start` and `t_end`
    
    # Make sure this step is done after dropping rows
    dat_prc %<>% dplyr::mutate(id=as.integer(factor(id)))
    
    # Set data attributes
    attr(dat_prc, "n") <- as.integer(max(dat_prc$id))
    # attr(dat_prc, "T_minus") <- dat_grp$T_minus
    # attr(dat_prc, "T_plus") <- dat_grp$T_plus
    attr(dat_prc, "s_i") <- dat_grp$s_i
    attr(dat_prc, "t_i") <- dat_grp$t_i
    # attr(dat_prc, "case") <- dat_grp$case
    
    # Drop rows with duplicate time
    # !!!!! Maybe move this above
    nrow(dat_prc)
    dupe_time_rows <- which(
      dat_prc$id==c(NA,dat_prc$id[1:length(dat_prc$id)-1]) &
        dat_prc$t_end==c(NA,dat_prc$t_end[1:length(dat_prc$t_end)-1])
    )
    if (length(dupe_time_rows)>0) { dat_prc <- dat_prc[-dupe_time_rows,] }
    nrow(dat_prc)
    
    # Drop rows with DOB > t_start
    # !!!!! TO DO

    # Generate delta column
    delta <- rep(NA, nrow(dat_prc))
    for (id in c(1:attr(dat_prc, "n"))) {
      rows_i <- which(dat_prc$id==id)
      # print(paste("id:", id)) # For debugging
      delta_i <- g_delta(
        case = dat_grp$case[id],
        s_i = dat_grp$s_i[id],
        t_i = dat_grp$t_i[id],
        T_minus = dat_grp$T_minus[id],
        T_plus = dat_grp$T_plus[id]
      )
      if (length(rows_i)!=length(delta_i)) {
        stop(paste0("Error with computation of delta for ID ", id, "."))
      }
      delta[rows_i] <- delta_i
    }
    dat_prc$delta <- delta
    
    # Create V (testing) and U (positive/known) indicators
    dat_prc %<>% dplyr::mutate(
      v = In(!is.na(dat_prc$yearoftest)),
      u = In(!is.na(first_hiv_pos_dt) & t_end>=first_hiv_pos_dt)
    )
    
    # Rescale time variable to start at 1
    dat_prc %<>% dplyr::mutate(
      dob = (dob - window_start) + 1,
      dod = (dod - window_start) + 1,
      t_start = (t_start - window_start) + 1,
      t_end = (t_end - window_start) + 1,
      yearoftest = (yearoftest - window_start) + 1,
      first_hiv_pos_dt = (first_hiv_pos_dt - window_start) + 1,
      last_hiv_neg_dt = (last_hiv_neg_dt - window_start) + 1,
      earliestartinitdate = (earliestartinitdate - window_start) + 1
    )
    attr(dat_prc, "s_i") <- (attr(dat_prc, "s_i") - window_start) + 1
    attr(dat_prc, "t_i") <- (attr(dat_prc, "t_i") - window_start) + 1
    
    # Create (scaled) baseline age variable
    # !!!!! Check this later; some with age -1 = -.01
    # Consider rounding w_2 as well (for memoising)
    dat_prc %<>% dplyr::mutate(w_2 = (t_start-dob)/100)
    
    # DQA
    if (F) {
      
      set.seed(1)
      case_1_ids <- sample(dplyr::filter(dat_grp, case==1)$id, size=5)
      case_2_ids <- sample(dplyr::filter(dat_grp, case==2)$id, size=5)
      case_3_ids <- sample(dplyr::filter(dat_grp, case==3)$id, size=5)
      case_4_ids <- sample(dplyr::filter(dat_grp, case==4)$id, size=5)
      
    }
    
    # !!!!! Check all variables for a handful of each case type
    
    # !!!!! Drop variables no longer needed (after QA)
    
    # TO DO: change "baseline age" to "age" (i.e. make it time-varying)
    
    # # TEMP: filter out some observations for individuals who migrated out
    # if (T) {
    #   
    #   # !!!!! TO DO
    #   
    # }
    # 
    # if (cfg2$save_data) { saveRDS(dat, "dat.rds") }
    
    # !!!!! These checks should be done after dataset processing is complete
    id_1 <- c(as.integer(min(dat_prc$id)):as.integer(max(dat_prc$id)))
    id_2 <- sort(unique(dat_prc$id))
    if (!identical(id_1, id_2)) { stop("Issue with ID numbers.") }
    
    dat <- dat_prc
    
    # # Save datasets for validation
    # write.table(dat, file="dat.csv", sep=",", row.names=FALSE)
    # write.table(dat_prc, file="dat_prc.csv", sep=",", row.names=FALSE)
    # write.table(dat_raw, file="dat_raw.csv", sep=",", row.names=FALSE)
    
    cols_to_drop <- c(
      "iintid", "dob", "dod", "age_start", "age_end", "yearoftest", "resultdate",
      "hivresult", "first_hiv_pos_dt", "last_hiv_neg_dt", "hiv_result_fill",
      "earliestartinitdate"
    )
    for (col in cols_to_drop) { dat[[col]] <- NULL }
    
    rm(dat_raw,dat_prc)
    
    # Check estimates for model 10 against Cox model estimates
    # !!!!! Move this code elsewhere
    if (F) {
      
      dat$j <- dat$t_end/10 # Create calendar time variable
      dat$w_2b <- dat$w_2+0.01 # Create time variable

      library(survival)
      model <- coxph(
        formula = Surv(w_2, w_2b, y)~z+j+w_1,
        data = dat
      )
      summary(model)
      
      # BEFORE REMOVING 90+ YEAR OLDS
      # n= 141817, number of events= 1721 
      # 
      #         coef exp(coef) se(coef)       z Pr(>|z|)    
      # z    0.45797   1.58086  0.09468   4.837 1.32e-06 ***
      # j   -0.61830   0.53886  0.03810 -16.228  < 2e-16 ***
      # w_1  0.34421   1.41088  0.04923   6.992 2.72e-12 ***
      # ---
      # 
      # AFTER REMOVING 90+ YEAR OLDS
      # n= 141134, number of events= 1681 
      # 
      #         coef exp(coef) se(coef)       z Pr(>|z|)    
      # z    0.47341   1.60546  0.09483   4.992 5.97e-07 ***
      # j   -0.64619   0.52404  0.03859 -16.746  < 2e-16 ***
      # w_1  0.34264   1.40866  0.04971   6.893 5.46e-12 ***
      # ---

      bh <- survival::basehaz(model, centered=FALSE)
      
      ggplot(bh, aes(x=time, y=hazard)) + geom_line()
      
      # Death rates by age
      dat2 <- dat %>% dplyr::group_by(w_2) %>% dplyr::summarize(
        num_persontime = n(),
        num_deaths = sum(y),
        death_rate = round(mean(y),3)
      )
      ggplot(dat2, aes(x=w_2, y=death_rate)) + geom_point()
      
    }
    
  } else {
    
    dat <- readRDS("dat.rds")
    
  }
  
}
chk(1, "Data reading/processing: END")



###############################.
##### Data quality checks #####
###############################.

if (cfg2$run_dqa) {
  
  chk(2, "DQA: START")
  
  # Setup
  dqa <- function(test) { if (test==F) { stop("DQA Error") } }
  dat_grp2 <- dat_raw %>% dplyr::group_by(iintid) %>%
    dplyr::summarize(
      unique_T_plus = (function(vec){ length(unique(vec)) })(first_hiv_pos_dt),
      unique_T_minus = (function(vec){ length(unique(vec)) })(last_hiv_neg_dt)
    )
  
  # Tests on raw data
  dqa(identical(sort(unique(dat_raw$sex)), c("Female", "Male")))
  dqa(sum(dat_raw$died==0, na.rm=T) + sum(dat_raw$died==1, na.rm=T)
      == length(unique(dat_raw$iintid)))
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
  
  chk(2, "DQA: END")
  
}



#########################.
##### Data analysis #####
#########################.

if (cfg2$run_analysis) {
  
  # Construct log likelihood function
  chk(3, "construct_negloglik: START")
  if (cfg2$parallelize) {
    n_cores <- cfg$sim_n_cores
    print(paste0("Using ", n_cores, " cores."))
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, ls(.GlobalEnv))
    negloglik <- construct_negloglik(dat, parallelize=T, cl=cl,
                                     cfg$model_version)
  } else {
    negloglik <- construct_negloglik(dat, parallelize=F, cl=NULL,
                                     cfg$model_version)
  }
  chk(3, "construct_negloglik: END")
  
  # Set initial parameter estimates
  if (cfg$model_version==0) {
    # This version not yet working
    par_init <- c(a_x=-5.603, g_x1=0, g_x2=-0.3655, a_y=-6.020, g_y1=0, g_y2=4.282, beta_x=1.401, beta_z=0.0004, t_x=0.9609, t_y=-3.906, a_s=-1.740, t_s=-2.100, g_s1=0, g_s2=1.271) # Model iteration 0
  } else if (cfg$model_version==1) {
    par_init <- c(a_x=-5.651, a_y=-4.942, beta_x=1.423, beta_z=0.2235, a_s=-2.007)
  } else if (cfg$model_version==2) {
    par_init <- c(a_x=-5.651, a_y=-4.942, beta_x=1.423, beta_z=0.2235, a_s=-2.007, g_x1=0, g_y1=0, g_s1=0)
  } else if (cfg$model_version %in% c(3,4)) {
    par_init <- c(a_x=-5.8039, g_x1=-0.5518, g_x2=0.6733, a_y=-6.4375, g_y1=0.3011, g_y2=4.0686, beta_x=1.6762, beta_z=0.8045, a_s=-2.116, g_s1=-0.3937, g_s2=1.1987)
  } else if (cfg$model_version==5) {
    par_init <- c(a_x=-5.8039, g_x1=-0.5518, g_x2=0.6733, a_y=-6.4375, g_y1=0.3011, g_y2=4.0686, beta_x=1.6762, beta_z=0.8045, t_y=0, a_s=-2.116, g_s1=-0.3937, g_s2=1.1987)
  } else if (cfg$model_version==6) {
    par_init <- c(a_x=-4.6185, g_x1=-1.3117, g_x2=-0.3883, a_y=-6.1581, g_y1=0.3794, g_y2=4.4762, beta_x=1.8497, beta_z=1.708, t_x=0, t_y=-0.45141, a_s=-2.495, g_s1=-0.0177, g_s2=1.684)
  } else if (cfg$model_version==7) {
    par_init <- c(a_x=-6.2967, g_x1=-0.1535, g_x2=0.9796, a_y=-5.5786, g_y1=0.3278, g_y2=4.2046, beta_x=1.401, beta_z=1.3177, t_x=0.5343, t_y=-0.7198, a_s=-2.3111, t_s=0.4245, g_s1=-0.5649, g_s2=0.6198)
  } else if (cfg$model_version==8) {
    par_init <- c(a_x=-3.5607, g_x1=-0.3244, g_x2=-0.2809, a_y=-5.7446, g_y1=0.3544, g_y2=4.4057, g_y3=0, g_y4=0, beta_x=1.8096, beta_z=1.8153, t_x=-0.786, t_y=-0.7826, a_s=-2.87, t_s=0.6349, g_s1=-0.3768, g_s2=0.6409)
  } else if (cfg$model_version==9) {
    par_init <- c(a_x=-2.2308, g_x1=-0.4977, g_x2=-0.9101, a_y=-6.3404, g_y1=0.5996, g_y2=2.4098, g_y3=2.8139, g_y4=7.0955, g_y5=6.0127, beta_x=1.295, beta_z=1.2856, t_x=-1.366, t_y=-0.7141, a_s=-2.0978, t_s=0.3321, g_s1=-0.8771, g_s2=0.8316)
  } else if (cfg$model_version==10) {
    par_init <- c(beta_z=0.3, a_y=-9.4955, g_y1=0.3209, g_y2=5.7549, g_y3=5.2759, g_y4=13.7284, g_y5=5.2979, t_y=-0.6637)
  } else if (cfg$model_version==11) {
    par_init <- c(a_x=-2.2308, g_x1=-0.4977, g_x2=0, g_x3=0, g_x4=0, g_x5=0, t_x=-1.366, a_s=-2.0978, g_s1=-0.8771, g_s2=0.8316, t_s=0.3321, beta_x=1.295, beta_z=1.2856, a_y=-6.3404, g_y1=0.5996, g_y2=2.4098, g_y3=2.8139, g_y4=7.0955, g_y5=6.0127, t_y=-0.7141)
  }
  
  # par_true <- c(
  #   par_true_full$a_x, par_true_full$g_x[1], par_true_full$g_x[2],
  #   par_true_full$a_y, par_true_full$g_y[1], par_true_full$g_y[2],
  #   par_true_full$beta_x, par_true_full$beta_z, par_true_full$t_x,
  #   par_true_full$t_y, par_true_full$a_s, par_true_full$t_s,
  #   par_true_full$g_s[1], par_true_full$g_s[2]
  # )
  
  # Run optimizer
  chk(4, "optim: START")
  opt <- stats::optim(
    par = par_init,
    fn = negloglik,
    method = "Nelder-Mead",
    control = list(maxit=avi_maxit,
                   reltol=avi_reltol)) # !!!!!
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
      r = avi_r, # r gives the number of Richardson improvement iterations (repetitions with successly smaller d
      v = 2 # v gives the reduction factor
    )
  )
  print(hessian_est)
  hessian_inv <- solve(hessian_est)
  chk(5, "hessian: END")
  
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
  #     avi_maxit = avi_maxit,
  #     avi_reltol = avi_reltol,
  #     avi_samp_size = samp_size,
  #     avi_r = avi_r,
  #     res = res
  #   ),
  #   file = paste0("res_", runif(1), ".rds")
  # )
  
}

chk(6, "END")
