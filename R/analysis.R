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
cfg2 <- list(
  process_data = T,
  save_data = T,
  run_dqa = F,
  run_analysis = T,
  parallelize = T,
  use_simulated_dataset = F
)

# !!!!! TEMPORARY: testing different optim config options
if (T) {
  avi_maxit <- as.integer(Sys.getenv("avi_maxit"))
  avi_reltol <- as.numeric(Sys.getenv("avi_reltol"))
  print("OPTIM CONFIG")
  print("------------")
  print(paste("maxit:", avi_maxit))
  print(paste("reltol:", avi_reltol))
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
    dat_prc <- dat_raw <- read.csv("../Data/data_raw_1000_v2.csv")
    
    # Drop unnecessary columns
    # dat_prc %<>% subset(select=-c(enter, exit, X_t, X_t0))
    dat_prc %<>% subset(select=-c(enter, exit))
    # !!!!! Drop hiv_result_fill, earliestartinitdate ?????
    
    # Start time
    window_start <- min(dat_prc$month)
    
    # Function to convert dates
    date_start <- as.Date("01jan1960", format="%d%b%Y")
    convert_date <- memoise(function(date) {
      date <- as.Date(date, format="%d%b%Y")
      # Subtract window_start ?????
      return(lubridate::interval(start=date_start, end=date) %/% months(1))
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
      month_prev = month - 1
    )
    
    # Filter out records with a negative test after a positive test
    print(nrow(dat_prc))
    dat_prc %<>% dplyr::filter(
      is.na(first_hiv_pos_dt) | is.na(last_hiv_neg_dt) |
        first_hiv_pos_dt>last_hiv_neg_dt
    )
    print(nrow(dat_prc))
    
    # Rearrange columns
    dat_prc %<>% dplyr::relocate(month_prev, .before=month)
    
    # Rename columns
    dat_prc %<>% dplyr::rename(
      "w_1" = sex,
      "w_2" = dob,
      "y" = died,
      "z" = onart,
      "t_start" = month_prev,
      "t_end" = month
    )
    
    # !!!!! Standardize all variables to have [0,1] range before adding to model
    
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
    
    # !!!!! Remove all data before 13th birthday; check for any testing before this
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
    print(nrow(dat_prc))
    dupe_time_rows <- which(
      dat_prc$id==c(NA,dat_prc$id[1:length(dat_prc$id)-1]) &
        dat_prc$t_end==c(NA,dat_prc$t_end[1:length(dat_prc$t_end)-1])
    )
    if (length(dupe_time_rows)>0) { dat_prc <- dat_prc[-dupe_time_rows,] }
    print(nrow(dat_prc))
    
    # Generate delta column
    delta <- rep(NA, nrow(dat_prc))
    for (id in c(1:attr(dat_prc, "n"))) {
      rows_i <- which(dat_prc$id==id)
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
      v = In(!is.na(dat_prc$monthoftest)),
      u = In(!is.na(first_hiv_pos_dt) & t_end<=first_hiv_pos_dt)
    )
    
    # # Save for validation
    # write.table(dat_prc, file="dat_prc.csv", sep=",", row.names=FALSE)
    # write.table(dat_raw, file="dat_raw.csv", sep=",", row.names=FALSE)
    
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
    cols_to_drop <- c(
      "iintid", "dod", "age_start", "age_end", "monthoftest", "resultdate",
      "hivresult", "first_hiv_pos_dt", "last_hiv_neg_dt", "hiv_result_fill",
      "earliestartinitdate"
    )
    for (col in cols_to_drop) { dat[[col]] <- NULL }
    rm(dat_raw,dat_prc)
    
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
  chk(3, "construct_negloglik_miss: START")
  if (cfg2$parallelize) {
    n_cores <- parallel::detectCores() - 1
    n_cores <- 10 # !!!!!
    # assign("fn_calls", 0, envir=.GlobalEnv) # !!!!!
    print(paste0("Using ", n_cores, " cores."))
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, ls(.GlobalEnv))
    negloglik_miss <- construct_negloglik_miss(dat, parallelize=T, cl=cl)
  } else {
    negloglik_miss <- construct_negloglik_miss(dat, parallelize=F, cl=NULL)
  }
  chk(3, "construct_negloglik_miss: END")
  
  # Set initial parameter estimates
  par_init <- log(c(
    a_x=0.003, g_x1=1.2, g_x2=1.1, a_y=0.002, g_y1=1.3, g_y2=1, beta_x=1.3,
    beta_z=0.8, t_x=0.999, t_y=1.001, a_s=0.03, t_s=1.001, g_s1=1.8, g_s2=1.7
  ))
  # par_true <- c(
  #   par_true_full$a_x, par_true_full$g_x[1], par_true_full$g_x[2],
  #   par_true_full$a_y, par_true_full$g_y[1], par_true_full$g_y[2],
  #   par_true_full$beta_x, par_true_full$beta_z, par_true_full$t_x,
  #   par_true_full$t_y, par_true_full$a_s, par_true_full$t_s,
  #   par_true_full$g_s[1], par_true_full$g_s[2]
  # )
  
  # Run optimizer
  chk(4, "optim: START")
  # opt_miss <- stats::optim(par=par_init, fn=negloglik_miss) # !!!!!
  
  opt_miss <- stats::optim(
    par = par_init,
    fn = negloglik_miss,
    method = "Nelder-Mead",
    control = list(maxit=avi_maxit,
                   reltol=avi_reltol)) # !!!!!
  
  # print(paste0("objective function calls (optim): ", fn_calls)) # !!!!!
  if (F) {
    stats::optim(par=par_init, fn=negloglik_miss, control=list(trace=6))
    library(optimParallel)
    opt_miss <- stats::optim(par=par_init, fn=negloglik_miss)
  }
  chk(4, "optim: END")
  
  # Compute Hessian
  chk(5, "hessian: START")
  # hessian_miss <- numDeriv::hessian(func=negloglik_miss, x=opt_miss$par)
  hessian_miss <- numDeriv::hessian(
    func = negloglik_miss,
    x = opt_miss$par,
    method = "Richardson", # "Richardson" "complex"
    method.args = list(
      eps = 1e-4,
      d = 0.1,
      zero.tol = sqrt(.Machine$double.eps/7e-7),
      r = 4,
      v = 2,
      show.details = F
    )
  )
  hessian_inv <- solve(hessian_miss)
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
    est <- as.numeric(opt_miss$par[i])
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
  # print(paste0("objective function calls (total): ", fn_calls)) # !!!!!
  
  # !!!!! temp
  saveRDS(
    list(
      avi_maxit = avi_maxit,
      avi_reltol = avi_reltol,
      res = res
    ),
    file = paste0("res_", runif(1))
  )
  
}

chk(6, "END")
