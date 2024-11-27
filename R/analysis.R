##################.
##### Config #####
##################.

set.seed(1)

for (pkg in c(cfg$pkgs,cfg$pkgs_nocluster)) {
  suppressMessages({ do.call("library", list(pkg)) })
}

chk(0, "START")
.t_start <- Sys.time()
cfg2 <- list(
  process_data = T,
  save_data = T,
  run_dqa = F,
  run_analysis = T,
  parallelize = T,
  use_simulated_dataset = F,
  # samp_size = 20000,
  # samp_size = as.integer(Sys.getenv("samp_size")),
  samp_size = 0, # Full sample
  opt_maxit = 5000,
  opt_r = 2,
  opt_reltol = 1e-5,
  window_start = 2010, # Spline bases for year must reflect this
  window_end = 2022, # Spline bases for year must reflect this
  age_end = 60, # Spline bases for age must reflect this
  model_sex = Sys.getenv("model_sex")
  # temp = T
)

if (T) {
  print("CONFIG")
  print("------------")
  print(paste("model_version:", cfg$model_version))
  print(paste("maxit:", cfg2$opt_maxit))
  print(paste("reltol:", cfg2$opt_reltol))
  print(paste("r:", cfg2$opt_r))
  print(paste("sample size:", cfg2$samp_size))
  print(paste("sex:", cfg2$model_sex))
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
    
    set.seed(1)
    
    # Set up data frame to track filtering
    dat_log <- data.frame("note"=character(), "num_rows"=integer())
    log_note <- function(note, num_rows) {
      dat_log[nrow(dat_log)+1,] <<- list(note, num_rows)
    }
    
    # Read in data
    # dat_prc <- dat_raw <- read.csv("../Data/data_raw_full_v2.csv")
    dat_prc <- dat_raw <- read.csv("../Data/pip_hiv_missingness_Avi_2024-11-21_mod.csv")
    log_note("# rows, original", nrow(dat_prc))
    log_note("# individuals, original", length(unique(dat_prc$IIntId)))
    
    # Take sample from dataset (for model development)
    if (cfg2$samp_size!=0) {
      iintids <- unique(dat_prc$IIntId)
      iintids_sample <- sample(iintids, size=cfg2$samp_size)
      dat_prc %<>% dplyr::filter(IIntId %in% iintids_sample)
    }
    log_note("# rows, subsample", nrow(dat_prc))
    
    # Rename columns
    dat_prc %<>% dplyr::rename(
      "id" = IIntId,
      "year" = Year
    )
    
    # Drop unnecessary columns
    dat_prc %<>% subset(select=-c(X, LocationId, EarliestARTInitDate, HIV_update, age_start,
                                  age_end, first_hiv_pos, last_hiv_neg))

    # Function to convert dates
    convert_date <- memoise(function(date) {
      # date <- as.Date(date, format="%d%b%Y")
      date <- as.Date(date, format="%Y-%m-%d")
      return(lubridate::year(date))
    })
    
    # Filter out obs with sex=="Unknown"
    rows_pre <- nrow(dat_prc)
    dat_prc %<>% dplyr::filter(sex!="Unknown")
    log_note("# rows dropped, sex=='Unknown'", rows_pre-nrow(dat_prc))
    
    # Filter out obs with missing `PIPSA`
    rows_pre <- nrow(dat_prc)
    dat_prc %<>% dplyr::filter(!is.na(PIPSA))
    log_note("# rows dropped, missing `PIPSA`", rows_pre-nrow(dat_prc))
    
    # Filter dataset based on sex
    rows_pre <- nrow(dat_prc)
    dat_prc %<>% dplyr::filter(sex==cfg2$model_sex)
    if (nrow(dat_prc)==0) { stop("`model_sex` incorrectly specified.") }
    log_note("# rows after filtering by sex", rows_pre-nrow(dat_prc))
    
    # Misc data wrangling
    dat_prc %<>% dplyr::mutate(
      died = ifelse(is.na(died), 0, died),
      dob = convert_date(dob),
      age = year - dob,
      ART_update = ifelse(is.na(ART_update), 0, ART_update)
    )
    
    # Filter out children with tests before age 12
    rows_pre <- nrow(dat_prc)
    children_with_tests <- dplyr::filter(
      dat_prc, age<=13 & ResultDate!=""
    )$id
    dat_prc %<>% dplyr::filter(!(id %in% children_with_tests))
    log_note("# rows dropped, children with tests before age 12",
             rows_pre-nrow(dat_prc))

    # Remove all data before 13th birthday
    rows_pre <- nrow(dat_prc)
    dat_prc %<>% dplyr::filter(age>=13)
    log_note("# rows dropped, age<13", rows_pre-nrow(dat_prc))
    
    # Remove all data after age cfg2$age_end
    rows_pre <- nrow(dat_prc)
    dat_prc %<>% dplyr::filter(age<cfg2$age_end)
    log_note("# rows dropped, age>=cfg2$age_end", rows_pre-nrow(dat_prc))
    
    # Remove tests with an "S" result
    rows_with_s <- which(dat_prc$HIVResult=="S")
    dat_prc[rows_with_s, "ResultDate"] <- ""
    dat_prc[rows_with_s, "HIVResult"] <- ""
    
    # Create flag for individuals whose status is positive/known at window_start
    dat_pos_at_start <- dat_prc %>%
      dplyr::filter(year<=cfg2$window_start) %>%
      mutate(pos=ifelse(is.na(HIVResult), 0, In(HIVResult=="P"))) %>%
      dplyr::group_by(id) %>%
      dplyr::summarize(pos_at_start = max(pos)) %>%
      dplyr::filter(pos_at_start==1)
    ids_pos_at_start <- dat_pos_at_start$id
    rm(dat_pos_at_start)
    
    # Remove all data prior to window_start
    rows_pre <- nrow(dat_prc)
    dat_prc %<>% dplyr::filter(year>=cfg2$window_start)
    log_note("# rows dropped, year<cfg2$window_start", rows_pre-nrow(dat_prc))
    
    # Set V=1 if status is known at window_start
    rows <- which(
      dat_prc$id %in% ids_pos_at_start & dat_prc$year==cfg2$window_start
    )
    for (row in rows) {
      hiv_res <- dat_prc[row,"HIVResult"]
      if (is.na(hiv_res) || hiv_res!="P") {
        dat_prc[row,"HIVResult"] <- "P"
        dat_prc[row,"ResultDate"] <- as.Date(
          paste0(cfg2$window_start, "-01-01"),
          format = "%Y-%m-%d"
        )
      }
    }
    
    # Remove all data following window_end
    rows_pre <- nrow(dat_prc)
    dat_prc %<>% dplyr::filter(year<=cfg2$window_end)
    log_note("# rows dropped, year>cfg2$window_end", rows_pre-nrow(dat_prc))
    
    # Add `first_hiv_pos_dt` and `last_hiv_neg_dt`
    dat_prc %<>% dplyr::mutate(
      hiv_pos_dts = ifelse(HIVResult=="P", year, NA),
      hiv_pos_dts = ifelse(is.na(hiv_pos_dts), 9999, hiv_pos_dts),
      hiv_neg_dts = ifelse(HIVResult=="N", year, NA),
      hiv_neg_dts = ifelse(is.na(hiv_neg_dts), 0, hiv_neg_dts),
      art_pos_dts = ifelse(ART_update==1, year, NA),
      art_pos_dts = ifelse(is.na(art_pos_dts), 9999, art_pos_dts)
    )
    dat_prc %<>% dplyr::mutate(
      first_hiv_pos_dt = min(hiv_pos_dts),
      first_hiv_pos_dt = ifelse(first_hiv_pos_dt==9999, NA, first_hiv_pos_dt),
      last_hiv_neg_dt = max(hiv_neg_dts),
      last_hiv_neg_dt = ifelse(last_hiv_neg_dt==0, NA, last_hiv_neg_dt),
      first_art_pos_dt = min(art_pos_dts),
      .by = id
    )
    dat_prc[["hiv_pos_dts"]] <- NULL
    dat_prc[["hiv_neg_dts"]] <- NULL
    dat_prc[["art_pos_dts"]] <- NULL
    
    # Create ART status variable
    dat_prc %<>% dplyr::mutate(
      ART_status_new = ifelse(year>=first_art_pos_dt, 1, 0),
      .by = id
    )
    
    # Filter out records with a negative test after a positive test
    rows_pre <- nrow(dat_prc)
    dat_prc %<>% dplyr::filter(
      is.na(first_hiv_pos_dt) | is.na(last_hiv_neg_dt) |
        first_hiv_pos_dt>last_hiv_neg_dt
    )
    log_note("# rows dropped, NEG test after POS test", rows_pre-nrow(dat_prc))

    # Drop rows with duplicate person-time
    rows_pre <- nrow(dat_prc)
    dupe_time_rows <- which(
      dat_prc$id==c(NA,dat_prc$id[1:length(dat_prc$id)-1]) &
        dat_prc$year==c(NA,dat_prc$year[1:length(dat_prc$year)-1])
    )
    if (length(dupe_time_rows)>0) { dat_prc <- dat_prc[-dupe_time_rows,] }
    log_note("# rows dropped, duplicate person-time", rows_pre-nrow(dat_prc))
    log_note("# rows, final", nrow(dat_prc))
    log_note("# individuals, final", length(unique(dat_prc$id)))
    log_note("# deaths, final", sum(dat_prc$died))
    
    # Print data log
    print(dat_log)
    
    # Sort dataframe
    dat_prc %<>% dplyr::arrange(id,year)
    
    # Rename columns
    dat_prc %<>% dplyr::rename(
      "y" = died,
      "z" = ART_status_new,
      "t_end" = year
    )
    
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
    
    # # !!!!! TEMPORARY
    # if (cfg2$temp) {
    #   case1_ids <- dplyr::filter(dat_grp, case==1)$id
    #   dat_grp %<>% dplyr::filter(case!=1)
    #   dat_prc %<>% dplyr::filter(!(id %in% case1_ids))
    # }
    
    # !!!!! Check that `first_hiv_pos_dt` and `last_hiv_neg_dt` lie within `t_start` and `t_end`
    
    # Drop rows with DOB > t_start
    # !!!!! TO DO

    # Renumber IDs
    dat_prc %<>% dplyr::mutate(
      id_orig = id,
      id = as.integer(factor(id))
    )
    
    # Set data attributes
    attr(dat_prc, "n") <- as.integer(max(dat_prc$id))
    attr(dat_prc, "s_i") <- dat_grp$s_i
    attr(dat_prc, "t_i") <- dat_grp$t_i
    
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
      v = In(!is.na(dat_prc$HIVResult)),
      u = In(!is.na(first_hiv_pos_dt) & t_end>=first_hiv_pos_dt)
    )
    
    # Rescale time variables via scale_time function
    dat_prc %<>% dplyr::mutate(
      dob = scale_time(dob, st=cfg2$window_start),
      t_end = scale_time(t_end, st=cfg2$window_start),
      first_hiv_pos_dt = scale_time(first_hiv_pos_dt, st=cfg2$window_start),
      last_hiv_neg_dt = scale_time(last_hiv_neg_dt, st=cfg2$window_start)
    )
    attr(dat_prc, "s_i") <- scale_time(attr(dat_prc, "s_i"),
                                       st=cfg2$window_start)
    attr(dat_prc, "t_i") <- scale_time(attr(dat_prc, "t_i"),
                                       st=cfg2$window_start)
    
    # Create scaled age variable (w_1) and geography covariate (w_2)
    # Geography covariate: 1="PIPSA North", 0="PIPSA South"
    dat_prc %<>% dplyr::mutate(
      w_1 = scale_age(t_end-dob),
      w_2 = In(PIPSA=="N")
    )
    
    # DQA
    if (F) {
      
      set.seed(1)
      case_1_ids <- sample(dplyr::filter(dat_grp, case==1)$id, size=5)
      case_2_ids <- sample(dplyr::filter(dat_grp, case==2)$id, size=5)
      case_3_ids <- sample(dplyr::filter(dat_grp, case==3)$id, size=5)
      case_4_ids <- sample(dplyr::filter(dat_grp, case==4)$id, size=5)
      
      id_1 <- c(as.integer(min(dat_prc$id)):as.integer(max(dat_prc$id)))
      id_2 <- sort(unique(dat_prc$id))
      if (!identical(id_1, id_2)) { stop("Issue with ID numbers.") }
      
    }
    
    dat <- dat_prc
    
    # # Save datasets for validation
    # write.table(dat, file="dat.csv", sep=",", row.names=FALSE)
    # write.table(dat_prc, file="dat_prc.csv", sep=",", row.names=FALSE)
    # write.table(dat_raw, file="dat_raw.csv", sep=",", row.names=FALSE)
    
    cols_to_drop <- c(
      "DoB", "dob", "age", "ResultDate", "HIVResult", "hiv_result_fill",
      "VisitDate", "ReceivedHIVTestResult", "CurrentlyOnART", "HadPosHIVResult",
      "first_hiv_pos_dt", "last_hiv_neg_dt", "ART_update", "first_art_pos_dt",
      "sex", "id_orig", "PIPSA", "year_begin", "year_end"
    )
    for (col in cols_to_drop) { dat[[col]] <- NULL }
    
    rm(dat_raw,dat_prc)
    
    # Create transformed dataset object
    dat_objs <- transform_dataset(
      dat = dat,
      model_version = cfg$model_version,
      window_start = cfg2$window_start,
      window_end = cfg2$window_end
    )
    
    if (cfg2$save_data) {
      saveRDS(dat, paste0("../Data/dat_", cfg2$model_sex, ".rds"))
      saveRDS(dat_objs, paste0("../Data/dat_objs_", cfg2$model_sex, ".rds"))
    }
    
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
    
    dat <- readRDS("../Data/dat.rds")
    dat_objs <- readRDS("../Data/dat_objs.rds")

  }
  
}
chk(1, "Data reading/processing: END")



###############################.
##### Data quality checks #####
###############################.

if (cfg2$run_dqa) {
  
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
    # par_init <- c(a_x=-6.0, g_x1=0.8, g_x2=-1.7, g_x3=2.1, g_x4=-1.5, g_x5=0, t_x1=0.9, t_x2=-1.5, t_x3=-1.3, t_x4=-0.8, a_s=-3.1, g_s1=3.5, g_s2=2.1, g_s3=4.3, g_s4=0.5, g_s5=0, beta_x1=1.6, beta_x2=-0.1, beta_x3=0.4, beta_x4=-1.6, a_y=-7.9, g_y1=2.2, g_y2=1.6, g_y3=4.6, g_y4=2.7, g_y5=0, t_y1=-0.3, t_y2=-0.3, t_y3=-0.7, t_y4=-0.1)
    par_init <- c(a_x=-5.2, g_x1=0.6, g_x2=-1.8, g_x3=0.3, g_x4=-2.0, t_x1=0.8, t_x2=-1.1, t_x3=-0.8, t_x4=-0.6, a_s=-3.1, g_s1=3.5, g_s2=2.0, g_s3=4.3, g_s4=0.6, beta_x1=2.0, beta_x2=0, beta_x3=-0.8, beta_x4=-0.1, a_y=-7.9, g_y1=2.0, g_y2=1.8, g_y3=4.6, g_y4=2.8, t_y1=-0.1, t_y2=-0.4, t_y3=-0.5, t_y4=0.3)
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
