#' Transform dataset into "counting process" format
#'
#' @param dataset A dataset returned by either generate_dataset() or
#'     perform_imputation()
#' @return A dataset in "counting process" format for Cox PH analysis

transform_dataset <- function(dataset) {
  
  # Data cleaning
  dataset$start_year <- attributes(dataset)$start_year
  dataset$death_year <- replace_na(
    dataset$death_year+1,
    attributes(dataset)$end_year
  )
  dataset %<>% rename("end_year" = death_year)
  
  # Put dataset into "counting process" format
  # !!!!! Make sure time is being split properly
  dataset_cp <- survSplit(
    formula = Surv(start_year, end_year, died) ~.,
    data = dataset,
    cut = c((attributes(dataset)$start_year):(attributes(dataset)$end_year))
  )
  
  # Create time-varying covariates
  # d1 <- dataset_cp %>% mutate(
  dataset_cp %<>% mutate(
    casc_status = case_when(
      is.na(sero_year) | start_year<sero_year ~ "HIV-",
      start_year>=sero_year & start_year<replace_na(art_init,9999) ~ "HIV+ART-",
      start_year>=sero_year & start_year>=art_init ~ "HIV+ART+",
      TRUE ~ "error"
    ),
    age = start_year - birth_year
  )
  # d2 <- dataset_cp %>% mutate(
  #   # dataset_cp %<>% mutate(
  #   casc_status = case_when(
  #     is.na(sero_year) | start_year<sero_year ~ "HIV-",
  #     start_year>=sero_year & start_year<replace_na(art_init,9999) ~ "HIV+ART-",
  #     start_year>=sero_year & start_year>=art_init ~ "HIV+ART+",
  #     TRUE ~ "error"
  #   ),
  #   age = start_year - birth_year
  # )
  
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
  
  return(dataset_cp)
  
}
