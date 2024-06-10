##################.
##### Config #####
##################.

for (pkg in c(cfg$pkgs,cfg$pkgs_nocluster)) {
  suppressMessages({ do.call("library", list(pkg)) })
}

chk(0, "START")
.t_start <- Sys.time()
cfg2 <- list(
  process_data = T, # !!!!!
  save_data = F,
  run_dqa = F,
  run_analysis = T,
  parallelize = T,
  use_simulated_dataset = F,
  samp_size = 20000, # as.integer(Sys.getenv("samp_size"))
  opt_maxit = 5000,
  opt_r = 2,
  opt_reltol = 1e-5
)

# !!!!! TEMPORARY: testing different optim config options
if (T) {
  print("CONFIG")
  print("------------")
  print(paste("model_version:", cfg$model_version))
  print(paste("maxit:", cfg2$opt_maxit))
  print(paste("reltol:", cfg2$opt_reltol))
  print(paste("r:", cfg2$opt_r))
  print(paste("sample size:", cfg2$samp_size))
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
    # cfg2$samp_size <- 5000 # !!!!! Rapid testing
    # cfg2$samp_size <- 150000 # !!!!! Testing code speed
    iintids_sample <- sample(iintids, size=cfg2$samp_size)
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
    
    # Remove tests with an "S" result
    rows_with_s <- which(dat_prc$hivresult=="S")
    dat_prc[rows_with_s, "yearoftest"] <- NA
    dat_prc[rows_with_s, "resultdate"] <- ""
    dat_prc[rows_with_s, "hivresult"] <- ""
    
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
    
    # # !!!!! TEMPORARY
    # if (Sys.getenv("sex")=="M") {
    #   dat_prc %<>% dplyr::filter(w_1==1)
    # } else if (Sys.getenv("sex")=="F") {
    #   dat_prc %<>% dplyr::filter(w_1==0)
    # }
    
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
    
    # Create transformed dataset object
    dat_objs <- transform_dataset(dat, model_version=cfg$model_version)
    
    saveRDS(dat, "dat.rds")
    saveRDS(dat_objs, "dat_objs.rds")
    
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
    dat_objs <- readRDS("dat_objs.rds")
    
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
  
  
  
  # Distribution of "year at entry"
  dat_prc %>%
    group_by(iintid) %>%
    mutate(min_year=min(year)) %>%
    filter(year==min_year) %>%
    xtabs(~year, data=.)
  
  # Distribution of "year of test"
  xtabs(~yearoftest, data=dat_prc)
  
  # Distribution of "year of first test"
  dat_prc %>%
    group_by(iintid) %>%
    filter(!is.na(yearoftest)) %>%
    mutate(min_year=min(year)) %>%
    filter(year==min_year) %>%
    xtabs(~year, data=.)
  
  # Distribution of "year of POS test"
  dat_prc %>%
    group_by(iintid) %>%
    filter(!is.na(yearoftest) & hivresult=="P") %>%
    xtabs(~year, data=.)
  
  # Distribution of "year of first POS test"
  dat_prc %>%
    group_by(iintid) %>%
    filter(!is.na(yearoftest) & hivresult=="P") %>%
    mutate(min_year=min(year)) %>%
    filter(year==min_year) %>%
    xtabs(~year, data=.)
  
  # Out of "year at entry" years, which years had HIV tests?
  dat_prc %>%
    group_by(iintid) %>%
    mutate(min_year=min(year)) %>%
    filter(year==min_year) %>%
    filter(!is.na(yearoftest)) %>%
    xtabs(~year, data=.)
  
  
  # !!!!! New DQA: checking "initial sero" model
  dat_2 <- dat %>%
    group_by(id) %>%
    mutate(min_time=min(t_start), max_time=max(t_start))
  dat_3 <- dplyr::filter(dat_2, t_start==min_time)
  dat_4 <- dplyr::filter(dat_3, v==1)
  dat_4 %<>% dplyr::mutate(
    age_below_50 = In(age_start<=50),
    year_00_07 = In(t_start<=7),
    year_08_16 = In(t_start>8 & t_start<=16),
    year_17_22 =  In(t_start>16)
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
    spl = which(substr(dat_i_names, 1, 1)=="b" & substr(dat_i_names, 3, 3)=="_")
  )
  
  # Construct log likelihood function
  chk(3, "construct_negloglik: START")
  if (cfg2$parallelize) {
    print(paste0("Using ", cfg$sim_n_cores, " cores."))
    n <- attr(dat, "n")
    folds <- cut(c(1:n), breaks=cfg$sim_n_cores, labels=FALSE)
    batches <- lapply(c(1:cfg$sim_n_cores), function(batch) {
      c(1:n)[folds==batch]
    })
    dat_objs_wrapper <- lapply(batches, function(i) { dat_objs[i] }) # !!!!! New
    cl <- parallel::makeCluster(cfg$sim_n_cores)
    objs_to_export <- c("f_x", "f_y", "icll", "lik_fn", "lik_fn2", "inds",
                        "batches")
    parallel::clusterExport(cl, objs_to_export, envir=.GlobalEnv)
    
    # parallel::clusterExport(cl, ls(.GlobalEnv))
    # negloglik <- construct_negloglik(parallelize=T, cl=cl,
    #                                  cfg$model_version)
    negloglik <- construct_negloglik(parallelize=T, cfg$model_version)
  } else {
    # negloglik <- construct_negloglik(parallelize=F, cl=NULL,
    #                                  cfg$model_version)
    negloglik <- construct_negloglik(parallelize=F, cfg$model_version)
  }
  chk(3, "construct_negloglik: END")
  
  # !!!!! Code profiling
  if (F) {
    
    if (F) {
      
      cl <- parallel::makeCluster(cfg$sim_n_cores)
      parallel::clusterEvalQ(cl, sink(paste0("C:/Users/ak811/Desktop/Avi/Research/HIVMI/output", Sys.getpid(), ".txt"))) # !!!!!
      parallel::parLapply(cl, c(1:5), function(i) {
        print("Check 1")
        print(pryr::mem_used())
      })
      parallel::clusterExport(cl, c("f_x", "f_y", "icll"), envir=.GlobalEnv)
      parallel::parLapply(cl, c(1:5), function(i) {
        print("Check 2")
        print(pryr::mem_used())
      })
      parallel::clusterExport(cl, c("lik_fn"), envir=.GlobalEnv)
      parallel::parLapply(cl, c(1:5), function(i) {
        print("Check 3")
        print(pryr::mem_used())
      })
      parallel::clusterExport(cl, c("inds", "batches"), envir=.GlobalEnv)
      parallel::parLapply(cl, c(1:5), function(i) {
        print("Check 4")
        print(pryr::mem_used())
      })
      print("Check 4.2")
      print("pryr::object_size(dat_objs)")
      print(pryr::object_size(dat_objs))
      print("pryr::object_size(get('dat_objs', envir=.GlobalEnv))")
      print(pryr::object_size(get("dat_objs", envir=.GlobalEnv)))
      parallel::clusterExport(cl, c("dat_objs"), envir=.GlobalEnv)
      parallel::parLapply(cl, c(1:5), function(i) {
        print("Check 5")
        print(pryr::mem_used())
        print("pryr::object_size(dat_objs)")
        print(pryr::object_size(dat_objs))
        # dat_objs[[i]]
      })
      
      parallel::clusterExport(cl, c("dat"), envir=.GlobalEnv)
      parallel::parLapply(cl, c(1:5), function(i) {
        pryr::object_size(dat)
        # pryr::mem_used()
      })
      
      
    }
    
    par_init <- c(a_x=-6.535, g_x1=-0.6737, g_x2=3.636, g_x3=0.2734, g_x4=0.4366, g_x5=-8.512, t_x1=-1.578, t_x2=-0.2818, t_x3=0.3750, t_x4=-2.160, a_s=-3.369, g_s1=-0.5761, g_s2=0.8899, t_s1=0, t_s2=0, t_s3=0, t_s4=0, beta_x=1, beta_z=0.5072, a_y=-5.944, g_y1=0.3940, g_y2=1.871, g_y3=2.923, g_y4=6.809, g_y5=3.004, t_y=-0.6077)
    
    negloglik <- construct_negloglik(parallelize=T, cfg$model_version)
    system.time({ qqq1 <- negloglik(par_init) })
    print(qqq1)
    
    negloglik <- construct_negloglik(parallelize=F, cfg$model_version)
    system.time({ qqq2 <- negloglik(par_init) })
    print(qqq2)
    
    # BEFORE: qqq = 17333.29
    # AFTER: qqq = 17333.29
  }
  
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
  } else if (cfg$model_version==12) {
    par_init <- c(a_x=-4.667, g_x1=-0.7444, g_x2=3.658, g_x3=-1.970, g_x4=-0.8989, g_x5=-6.690, t_x=-1.169, a_s=-3.299, g_s1=-0.6594, g_s2=0.8443, t_s=1.051, beta_x=0.4, beta_z=0.4, a_y=-5.654, g_y1=0.2733, g_y2=1.759, g_y3=2.802, g_y4=6.180, g_y5=2.649, t_y=-0.6332)
  } else if (cfg$model_version==13) {
    par_init <- c(a_x=-6.535, g_x1=-0.6737, g_x2=3.636, g_x3=0.2734, g_x4=0.4366, g_x5=-8.512, t_x1=-1.578, t_x2=-0.2818, t_x3=0.3750, t_x4=-2.160, a_s=-3.369, g_s1=-0.5761, g_s2=0.8899, t_s=1.051, beta_x=1, beta_z=0.5072, a_y=-5.944, g_y1=0.3940, g_y2=1.871, g_y3=2.923, g_y4=6.809, g_y5=3.004, t_y=-0.6077)
  } else if (cfg$model_version==14) {
    par_init <- c(a_x=-6.3567, g_x1=-0.7201, g_x2=3.5877, g_x3=0.8982, g_x4=1.2344, g_x5=-7.8761, t_x1=-2.1962, t_x2=-0.396, t_x3=-0.0074, t_x4=-2.2542, a_s=-3.0662, g_s1=-0.6222, g_s2=0.8341, t_s=0.8966, beta_x=0.9409, beta_z=0.7236, a_y=-6.5048, g_y1=0.4056, g_y2=1.9911, g_y3=3.0964, g_y4=6.9362, g_y5=3.4698, t_y1=-0.2982, t_y2=-0.8746, t_y3=-0.5772, t_y4=-1.0723)
  } else if (cfg$model_version==15) {
    par_init <- c(a_x=-6.67, g_x1=4.86, g_x2=1.33, g_x3=0.871, g_x4=-6.53, g_x5=3.66, g_x6=1.07, g_x7=2.60, g_x8=-7.59, t_x1=-1.82, t_x2=-0.728, t_x3=-0.593, t_x4=-2.06, a_s=-3.08, g_s1=-0.742, g_s2=0.753, t_s=0.936, beta_x=1.12, beta_z=1.00, a_y=-6.56, g_y1=0.431, g_y2=1.95, g_y3=3.29, g_y4=7.18, g_y5=3.60, t_y1=-0.319, t_y2=-1.00, t_y3=-0.934, t_y4=-1.12)
  } else if (cfg$model_version==16) {
    par_init <- c(a_x=-6.67, g_x1=4.86, g_x2=1.33, g_x3=0.871, g_x4=-6.53, g_x5=3.66, g_x6=1.07, g_x7=2.60, g_x8=-7.59, t_x1=-1.82, t_x2=-0.728, t_x3=-0.593, t_x4=-2.06, a_s=-3.08, g_s1=-0.742, g_s2=0.753, t_s1=0, t_s2=0, t_s3=0, t_s4=0, beta_x=1.12, beta_z=1.00, a_y=-6.56, g_y1=0.431, g_y2=1.95, g_y3=3.29, g_y4=7.18, g_y5=3.60, t_y1=-0.319, t_y2=-1.00, t_y3=-0.934, t_y4=-1.12)
  } else if (cfg$model_version==17) {
    par_init <- c(a_x=-6.67, g_x1=4.86, g_x2=1.33, g_x3=0.871, g_x4=-6.53, g_x5=3.66, g_x6=1.07, g_x7=2.60, g_x8=-7.59, t_x1=-1.82, t_x2=-0.728, t_x3=-0.593, t_x4=-2.06, a_s=-3.08, g_s1=-0.742, g_s2=0, g_s3=0, g_s4=0, g_s5=0, beta_x=1.12, beta_z=1.00, a_y=-6.56, g_y1=0.431, g_y2=1.95, g_y3=3.29, g_y4=7.18, g_y5=3.60, t_y1=-0.319, t_y2=-1.00, t_y3=-0.934, t_y4=-1.12)
  } else if (cfg$model_version==18) {
    par_init <- c(a_x=-6.9521, g_x1=4.2808, g_x2=-0.1476, g_x3=1.4616, g_x4=-5.7527, g_x5=2.4626, g_x6=0.1639, g_x7=3.6117, g_x8=-7.3954, t_x1=-1.093, t_x2=-0.888, t_x3=-1.3534, t_x4=-1.4756, a_s=-2.8117, g_s1=-0.3086, g_s2=0.8134, g_s3=-0.1185, g_s4=3.2074, g_s5=-1.4434, beta_x1=1.369, beta_x2=0, beta_z1=1.1501, beta_z2=0, a_y=-6.5127, g_y1=0.4162, g_y2=1.7359, g_y3=3.1651, g_y4=6.5974, g_y5=4.0223, t_y1=-0.1488, t_y2=-0.8869, t_y3=-0.6838, t_y4=-1.0119)
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
  counter <- 0
  nll_init <- negloglik(par_init)
  print(paste("negloglik(par_init):", nll_init))
  opt <- stats::optim(
    par = par_init,
    fn = negloglik,
    method = "Nelder-Mead",
    control = list(maxit=cfg2$opt_maxit,
                   reltol=cfg2$opt_reltol)) # !!!!!
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
