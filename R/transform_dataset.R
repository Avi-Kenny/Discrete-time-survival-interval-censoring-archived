#' Transform dataset into "counting process" format
#'
#' @param dat_baseline A dataset returned by generate_data_baseline()
#' @param dat_events A dataset returned by generate_data_events()
#' @return A dataset in "counting process" format for Cox PH analysis

transform_dataset <- function(dat_baseline, dat_events) {
  
  # Extract variables
  start_year <- attr(dat_baseline, "start_year")
  end_year <- attr(dat_events, "end_year")
  I <- nrow(dat_baseline)
  
  # Data transformation
  dat_baseline$start_time <- 0
  dat_baseline$end_time <- sapply(dat_events, function(d) { d$T_i })
  dat_baseline$y_copy <- sapply(dat_events, function(d) { d$y[d$T_i] })
  
  # Put dataset into "counting process" format
  dat_cp <- survSplit(
    formula = Surv(start_time, end_time, y_copy) ~.,
    data = dat_baseline,
    cut = c(1:(12*(end_year-start_year)))
  )
  
  # Convert dat_events to a dataframe and attach to dat_cp
  # cbind is functioning as an inner join since both dataframes are sorted
  df_ev <- as.data.frame(rbindlist(dat_events))
  df_ev %<>% filter(y!=9)
  df_ev %<>% subset(select=-id)
  dat_cp %<>% cbind(df_ev)
  
  # Create exposure variable
  dat_cp %<>% mutate(
    casc_status = case_when(
      x==0 ~ "HIV-",
      x==1 & z==0 ~ "HIV+ART-",
      x==1 & z==1 ~ "HIV+ART+"
    ),
    age = b_age + start_time/12
  )
  
  return(dat_cp)
  
}
