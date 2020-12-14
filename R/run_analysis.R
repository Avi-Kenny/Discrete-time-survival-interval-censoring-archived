#' Run Cox PH analysis
#'
#' @param dataset A dataset returned by either generate_dataset() or
#'     perform_imputation()
#' @param method One of c("ideal", "censor", "mi")
#' @return A list containing the following:
#'     est: point estimate of exposure coefficient
#'     se: standard error of exposure coefficient estimate

run_analysis <- function(dataset, method) {
  
  # Transform dataset
  {
    # Data cleaning
    dataset$start_year <- attributes(dataset)$start_year
    dataset$death_year <- replace_na(dataset$death_year,
                                     attributes(dataset)$end_year)
    dataset %<>% rename("end_year" = death_year)
    
    # Put dataset into "counting process" format
    dataset_cp <- survSplit(
      formula = Surv(start_year, end_year, died) ~.,
      data = dataset,
      cut = c((attributes(dataset)$start_year):(attributes(dataset)$end_year))
    )
    
    # Create time-varying covariates
    dataset_cp %<>% mutate(
      casc_status = ifelse(
        is.na(sero_year), "HIV-", ifelse(
          (start_year>=sero_year) & (start_year<replace_na(art_init,9999)),
          "HIV+ART-", "HIV+ART+")
      ),
      age = start_year - birth_year
    )
    
    # Create time-varying age bucket variable
    dataset_cp %<>% mutate(
      age_bin = factor(case_when(
        age %in% c(1:9) ~ 0,
        age %in% c(10:19) ~ 1,
        age %in% c(20:29) ~ 2,
        age %in% c(30:39) ~ 3,
        age %in% c(40:49) ~ 4,
        age %in% c(50:59) ~ 5,
        age %in% c(60:69) ~ 6,
        age %in% c(70:79) ~ 7,
        age %in% c(80:89) ~ 8,
        age %in% c(90:99) ~ 9,
        age %in% c(100:109) ~ 10
      ), levels=c(5,0,1,2,3,4,6,7,8,9,10) # Reference group is 50-59
    ))
    
  }
  
  # Create "censor" dataset
  if (type=="censor") {
    
    # !!!!! Clean this up, given changes to generate_dataset()
    
    # Set up censoring variable: c0.x == don't censor, c1.x == censor
    dataset_cp$censor <- "c0.0"
    
    # Don't censor any observations for baseline HIV+ART+ patients
    # !!!!! This should be redundant now
    dataset_cp %<>% mutate(
      censor = ifelse(baseline_status=="HIV+ART+", "c0.1", censor)
    )
    # xtabs(~censor, data=dataset_cp, addNA=T)
    
    # Censor all obs for patients with no testing data
    dataset_cp %<>% mutate(
      censor = ifelse(baseline_status!="HIV+ART+" & is.na(first_test),
                      "c1.1", censor)
    )
    
    # For baseline HIV- or HIV+ART- patients, censor obs before first test
    dataset_cp %<>% mutate(
      censor = ifelse(
        baseline_status!="HIV+ART+" & !is.na(first_test) &
          start_year<first_test, "c1.2", censor)
    )
    
    # For patients who only had negative tests, censor obs after last test
    dataset_cp %<>% mutate(
      censor = ifelse(
        baseline_status!="HIV+ART+" & !is.na(last_neg_test) &
          is.na(first_pos_test) & start_year>last_neg_test, "c1.3", censor)
    )

    # !!!!! "c0.2", "c0.3", etc.
    dataset_cp %<>% filter(censor %in% c("c0.0","c0.1"))
    
    # xtabs(~casc_status, data=dataset_cp_censor, addNA=T)

    # !!!!! Check censoring manually
  }
  
  # Fit time-varying Cox model ("ideal")
  fit <- coxph(
    Surv(start_year, end_year, died) ~ factor(casc_status) + age_bin +
      cluster(patient_id),
    data = dataset_cp
  )
  summ <- summary(fit)$coefficients
  
  results <- list(
    est_hiv = summ["factor(casc_status)HIV+ART-","coef"],
    se_hiv = summ["factor(casc_status)HIV+ART-","se(coef)"],
    est_art = summ["factor(casc_status)HIV+ART+","coef"],
    se_art = summ["factor(casc_status)HIV+ART+","se(coef)"]
  )
  
  return(results)
  
}
