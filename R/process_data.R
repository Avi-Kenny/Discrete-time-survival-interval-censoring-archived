##################.
##### Config #####
##################.

# Set configutation
chk(0, "START")
.t_start <- Sys.time()
cfg2 <- list(
  process_data = T,
  save_data = F,
  run_analysis = T,
  parallelize = T,
  use_simulated_dataset = F,
  # samp_size = 20000,
  samp_size = 0, # Full sample
  window_start = 2010,
  window_end = 2022,
  age_end = 60
)

# !!!!! DEBUGGING
if (T) {
  cfg2$parallelize <- F
  cfg2$samp_size <- 20000
}

# Load packages
for (pkg in c(cfg$pkgs,cfg$pkgs_nocluster)) {
  suppressMessages({ do.call("library", list(pkg)) })
}

# Print config info
print("CONFIG")
print("------------")
print(paste("model_version:", cfg$model_version))
print(paste("maxit:", cfg2$opt_maxit))
print(paste("reltol:", cfg2$opt_reltol))
print(paste("r:", cfg2$opt_r))
print(paste("sample size:", cfg2$samp_size))
print(paste("sex:", cfg$model_sex))
print("------------")



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
  dat <- generate_data(n=n, max_time=70, par=par_true_full)
  print(paste("n:",n))
  print(paste("rows in dataset:", nrow(dat)))
  
  # Manipulate dataset
  dat %<>% arrange(id, t_end)
  dat$x <- NULL
  attr(dat, "T_plus") <- NULL
  attr(dat, "T_minus") <- NULL
  attr(dat, "case") <- NULL
  attr(dat, "max_time") <- NULL
  attr(dat, "par") <- NULL
  
} else {
  
  if (cfg2$process_data) {
    
    set.seed(1)
    
    # Set up data frame to track filtering
    dat_log <- data.frame("note"=character(), "num_rows"=integer())
    log_note <- function(note, num_rows) {
      dat_log[nrow(dat_log)+1,] <<- list(note, num_rows)
    }
    
    # Read in data
    dat_prc <- read.csv("../Data/data_raw_full_v2.csv") # !!!!! Temporarily reverting to this version until Stephen adds 2023 back in
    dat_prc$PIPSA <- 1 # !!!!!! TEMP
    dat_prc$LocationId <- 1 # !!!!!! TEMP
    dat_prc$HIV_update <- 1 # !!!!!! TEMP
    dat_prc$age_start <- 1 # !!!!! TEMP
    dat_prc$age_end <- 1 # !!!!! TEMP
    dat_prc$first_hiv_pos <- 1 # !!!!! TEMP
    dat_prc$last_hiv_neg <- 1 # !!!!! TEMP
    dat_prc$HIV_status <- NULL # !!!!! TEMP
    dat_prc$dod <- NULL # !!!!! TEMP
    dat_prc %<>% dplyr::rename("ART_update"=ART_status) # !!!!! TEMP
    
    # dat_prc <- read.csv("../Data/pip_hiv_missingness_Avi_2024-11-21_mod.csv")
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
    dat_prc %<>% dplyr::filter(sex==cfg$model_sex)
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
    
    cols_to_drop <- c(
      "DoB", "dob", "age", "ResultDate", "HIVResult", "hiv_result_fill",
      "VisitDate", "ReceivedHIVTestResult", "CurrentlyOnART", "HadPosHIVResult",
      "first_hiv_pos_dt", "last_hiv_neg_dt", "ART_update", "first_art_pos_dt",
      "sex", "id_orig", "PIPSA", "year_begin", "year_end"
    )
    for (col in cols_to_drop) { dat[[col]] <- NULL }
    
    rm(dat_prc)
    
    # Create transformed dataset object
    dat_objs <- transform_dataset(
      dat = dat,
      model_version = cfg$model_version,
      window_start = cfg2$window_start,
      window_end = cfg2$window_end
    )
    
    if (cfg2$save_data) {
      saveRDS(dat, paste0("../Data/dat_", cfg$model_sex, ".rds"))
      saveRDS(dat_objs, paste0("../Data/dat_objs_", cfg$model_sex, ".rds"))
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
