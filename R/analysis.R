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
  process_data = TRUE,
  save_data = TRUE,
  run_dqa = TRUE,
  run_analysis = TRUE,
  parallelize = TRUE
)

###########################.
##### Data processing #####
###########################.

chk(1, "Data reading/processing: START")
if (cfg2$process_data) {
  
  if (F) {
    # Generate dataset
    dat <- generate_data(
      n = 200,
      # n = 2000,
      # n = 10000,
      max_time = 70,
      params = list(
        a_x=log(0.005), a_y=log(0.003), a_v=log(0.7), a_z=log(0.01),
        g_x=c(log(1.3),log(1.2)), g_y=c(log(1.2),log(1.1)),
        g_v=c(log(1.2),log(1.1)), g_z=c(log(1.2),log(1.1)),
        beta_x=log(1.5), beta_z=log(0.7),
        a_s=log(0.05), g_s=c(log(2),log(1.5))
      )
    )
    
    # Add x_prev column
    dat %<>% arrange(id, t_end)
    dat$x_prev <- ifelse(
      dat$id==c(0,dat$id[c(1:(length(dat$id)-1))]),
      c(0,dat$x[c(1:(length(dat$x)-1))]),
      0
    )
  }
  
  # Read in data; see generate_data() for expected columns
  dat_prc <- dat_raw <- read.csv("data_raw.csv")
  
  # Drop unnecessary columns
  dat_prc %<>% subset(select=-c(enter, exit, X_t, X_t0))
  # !!!!! Drop hiv_result_fill, earliestartinitdate ?????
  
  # Start time
  window_start <- min(dat_prc$month)
  
  # Function to convert dates
  date_start <- as.Date("01jan1960", format="%d%b%Y")
  # date_vars <- c("dob", "dod", "enter", "exit", "first_hiv_pos_dt",
  #                "last_hiv_neg_dt", "earliestartinitdate")
  convert_date <- memoise(function(date) {
    date <- as.Date(date, format="%d%b%Y")
    # Subtract window_start ?????
    return(lubridate::interval(start=date_start, end=date) %/% months(1))
  })
  
  # Data wrangling
  dat_prc %<>% dplyr::mutate(
    # iintid = as.numeric(factor(iintid)),
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
  
  # !!!!! Standardize all variables to have [0,1] range
  
  # Rearrange columns
  dat_prc %<>% dplyr::relocate(month_prev, .before=month)
  
  # Rename columns
  dat_prc %<>% dplyr::rename(
    "id" = iintid,
    "w_1" = sex,
    "w_2" = dob,
    "y" = died,
    "z" = onart,
    "t_start" = month_prev,
    "t_end" = month
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
  
  # !!!!! Remove all data before 13th birthday; check for any testing before this
  
  # !!!!! Temporary
  if (T) {
    bad_ids <- which(dat_grp$T_plus-dat_grp$T_minus<0)
    dat_prc %<>% dplyr::filter(!(id %in% bad_ids))
    dat_grp %<>% dplyr::filter(!(id %in% bad_ids))
  }
  
  # !!!!! Check that `first_hiv_pos_dt` and `last_hiv_neg_dt` lie within `t_start` and `t_end`
  
  # Set data attributes
  attr(dat_prc, "n") <- max(dat_prc$id)
  attr(dat_prc, "T_minus") <- dat_grp$T_minus
  attr(dat_prc, "T_plus") <- dat_grp$T_plus
  attr(dat_prc, "s_i") <- dat_grp$s_i
  attr(dat_prc, "t_i") <- dat_grp$t_i
  attr(dat_prc, "case") <- dat_grp$case
  
  # Generate delta column
  delta <- rep(NA, nrow(dat_prc))
  for (id in c(1:attr(dat_prc, "n"))) {
    index_first <- NA
    index_last <- NA
    delta_i <- g_delta(
      case = dat_grp$case[i],
      s_i = dat_grp$s_i[i],
      t_i = dat_grp$t_i[i],
      T_minus = dat_grp$T_minus[i],
      T_plus = dat_grp$plus[i]
    )
    delta[c(index_first:index_last)] <- delta_i
  }
  
  # !!!!!
  print(dat_prc[c(630:636),])
  
  #
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # TO DO: change "baseline age" to "age" (i.e. make it time-varying)
  

  # # Add x_prev column
  # dat %<>% arrange(id, t_end)
  # dat$x_prev <- ifelse(
  #   dat$id==c(0,dat$id[c(1:(length(dat$id)-1))]),
  #   c(0,dat$x[c(1:(length(dat$x)-1))]),
  #   0
  # )
  # 
  # # TEMP: filter out individuals who migrated out
  # if (T) {
  #   
  #   # !!!!! TO DO
  #   
  # }
  # 
  # if (cfg2$save_data) { saveRDS(dat, "dat.rds") }
  
} else {
  
  dat <- readRDS("dat.rds")
  
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
    print(paste0("Using ", n_cores, " cores."))
    # cl <- parallel::makeCluster(n_cores) # !!!!!
    cl <- parallel::makeCluster(35) # !!!!!
    parallel::clusterExport(cl, ls(.GlobalEnv))
    # cl <- NULL # !!!!!
    negloglik_miss <- construct_negloglik_miss(dat, parallelize=T, cl=cl)
  } else {
    negloglik_miss <- construct_negloglik_miss(dat, parallelize=F, cl=NULL)
  }
  chk(3, "construct_negloglik_miss: END")
  
  # Set initial parameter estimates
  par <- log(c(
    a_x=0.003, g_x1=1.2, g_x2=1.1, a_y=0.002, g_y1=1.3, g_y2=1, beta_x=1.3,
    beta_z=0.8, t_x=0.999, t_y=1.001, a_s=0.03, t_s=1.001, g_s1=1.8, g_s2=1.7
  ))
  
  # Run optimizer
  chk(4, "optim: START")
  opt_miss <- optim(par=par, fn=negloglik_miss)
  if (F) {
    optim(par=par, fn=negloglik_miss, control=list(trace=6))
    library(optimParallel)
    opt_miss <- optim(par=par, fn=negloglik_miss)
  }
  chk(4, "optim: END")
  
  # Compute Hessian
  chk(5, "hessian: START")
  hessian_miss <- hessian(func=negloglik_miss, x=opt_miss$par)
  hessian_inv <- solve(hessian_miss)
  chk(5, "hessian: END")
  
  # if (cfg2$parallelize) { stopCluster(cl) }
  
  # Parse results
  res <- data.frame(
    "param" = character(),
    "est" = double(),
    "se" = double(),
    "ci_lo" = double(),
    "ci_up" = double()
  )
  for (i in c(1:length(par))) {
    est <- as.numeric(opt_miss$par[i])
    se <- sqrt(diag(hessian_inv))[i]
    res[i,] <- c(
      param = par[i],
      est = est,
      se = se,
      ci_lo = est - 1.96*se, # Maybe do transformations later
      ci_up = est + 1.96*se # Maybe do transformations later
    )
  }
  print(res)
  
}


chk(6, "END")
