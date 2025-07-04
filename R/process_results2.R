#####################################.
##### Setup and data processing #####
#####################################.

# Set up data frame for tracking
dat_log <- data.frame("note"=character(), "num_rows"=integer())
log_note <- function(note, num_rows) {
  dat_log[nrow(dat_log)+1,] <<- list(note, num_rows)
}

# Config
c_date <- format(Sys.time(), "%Y-%m-%d")
status_vec <- c("HIV-", "HIV+")
sources <- c("Data (prc)", "Model (HIV-)", "Model (HIV+)")
colors_1 <- c("brown3", "cyan3", "cyan4")
# colors_2 <- c("cyan3", "cyan4", "brown3", "orange")
colors_2 <- colors_1
if (cfg$w_start==2010) {
  breaks_year <- c(2010,2015,2020)
  years_plot <- c(2010,2013,2016,2019,2022)
} else if (cfg$w_start==2017) {
  breaks_year <- c(2017,2022)
  years_plot <- c(2017,2018,2020,2021,2022)
}

# Load datasets
{
  # # dat_raw is raw/original dataset
  # dat_raw <- readstata13::read.dta13("../Data/Raw data/SurveillanceEpisodesHIV.dta")
  
  # dat_prc is the processed dataset via process_data.R
  # To save a copy, need to uncomment the line marked "Save datasets for
  # validation" and run this once for each sex
  dat_prc_M <- readRDS("../Data/dat_prc_M_20250205.rds")
  dat_prc_F <- readRDS("../Data/dat_prc_F_20250205.rds")
  dat_prc_M$sex <- 1
  dat_prc_F$sex <- 0
  dat_prc <- rbind(dat_prc_M, dat_prc_F)
}

# dat_prc: Rescale and rename variables
dat_prc %<>% dplyr::mutate(
  t_end = round(unscale_time(t_end, st=cfg$w_start)),
  age = unscale_age(w_1)
)

# # dat_raw: Generate age-at-death variable
# dat_raw$AoD <- round((dat_raw$DoD-dat_raw$DoB)/365.25)

# Summary stats
# log_note("rows, dat_raw", nrow(dat_raw))
# log_note("# people, dat_raw", length(unique(dat_raw$IIntId)))
log_note("rows, dat_prc", nrow(dat_prc))
log_note("# people, dat_prc", length(unique(dat_prc$id)))



#########################################################################.
##### Compare person-time exposure and death rates between datasets #####
#########################################################################.

# Prep work
age_bins <- list(
  c(0,12), c(13,20), c(21,30), c(31,40), c(41,50), c(51,60), c(61,999)
)
df_summ <- data.frame(
  "year" = integer(),
  "sex" = character(),
  "age_bin" = character(),
  "n_deaths" = integer(),
  "n_py" = double(),
  "rate" = double(),
  "source" = character()
)

# Populate df_summ
for (year_ in c(cfg$w_start:cfg$w_end)) {
  
  yr_s <- as.Date(x=paste0(year_,"-01-01"))
  yr_e <- as.Date(x=paste0(year_,"-12-31"))
  
  for (sex_ in c("Male", "Female")) {
    for (age_bin in age_bins) {
      
      # # Summary stats from raw dataset
      # dat_raw$age_at_start <- (
      #   as.Date(x=paste0(year_,"-01-01")) - dat_raw$DoB
      # ) / 365.25
      
      # dat_raw$n_days_exp <- pmax(
      #   pmin(dat_raw$EndDate,yr_e) - pmax(dat_raw$StartDate,yr_s) + 1,
      #   0
      # )
      # n_deaths_raw <- sum(
      #   dat_raw$Died=="Yes" &
      #     as.integer(substr(dat_raw$DoD,1,4))==year_ &
      #     dat_raw$Sex==sex_ &
      #     round(dat_raw$AoD) >= age_bin[1] &
      #     round(dat_raw$AoD) <= age_bin[2],
      #   na.rm=T
      # )
      # n_py_raw <- sum(
      #   as.numeric(dat_raw$n_days_exp/365.25) * (
      #     dat_raw$Sex==sex_ &
      #       round(dat_raw$age_at_start) >= age_bin[1] &
      #       round(dat_raw$age_at_start) <= age_bin[2]
      #   ),
      #   na.rm=T
      # )
      
      # Summary stats from processed dataset
      dat_prc_filt <- dplyr::filter(
        dat_prc,
        t_end==year_ & sex==as.integer(sex_=="Male") &
          age >= age_bin[1] & age <= age_bin[2]
      )
      n_deaths_prc <- sum(dat_prc_filt$y)
      n_py_prc <- nrow(dat_prc_filt)
      
      # # Summary from raw dataset
      # df_summ[nrow(df_summ)+1,] <- list(
      #   year_, sex_, paste0(age_bin, collapse="-"), n_deaths_raw, n_py_raw,
      #   round(1000*(n_deaths_raw/n_py_raw), 1), "raw"
      # )
      
      # Summary from PRC dataset
      df_summ[nrow(df_summ)+1,] <- list(
        year_, sex_, paste0(age_bin, collapse="-"), n_deaths_prc, n_py_prc,
        round(1000*(n_deaths_prc/n_py_prc), 1), "Data (prc)"
      )
      
      # Summary from model (age midpoint)
      # Note: if this code is uncommented, the code in sections 1 and 3 of
      #       process.R needs to be run first
      for (status in status_vec) {
        age_ <- mean(age_bin)
        type <- paste0("mort (", status, ")")
        rate <- round(1000 * prob(
          type=type, j=year_, w_1=age_, w_2=cfg$w_2, w_3=as.integer(sex_=="Male"),
          year_start=cfg$w_start, which="est"
        ), 1)
        df_summ[nrow(df_summ)+1,] <- list(
          year_, sex_, paste0(age_bin, collapse="-"), NA, NA,
          rate, paste0("Model (", status, ")")
        )
      }
      
    }
  }
}

# Plot death rates by year, filtering out highest age bin
# Export: 10 x 5
plot <- ggplot(
  dplyr::filter(df_summ, !(age_bin %in% c("0-12", "61-999"))),
  aes(x=year, y=rate, color=factor(source))
) +
  geom_line() +
  facet_grid(rows=dplyr::vars(sex), cols=dplyr::vars(age_bin)) +
  scale_x_continuous(breaks=breaks_year) +
  scale_color_manual(values=colors_1) +
  labs(color="Source", y="Deaths per 1,000 person-years")
ggsave(
  filename = paste0("../Figures + Tables/", c_date, " death_rates_by_year - ",
                    "model ", cfg$model_version, ".pdf"),
  plot = plot, device="pdf", width=10, height=5
)

# # Plot death counts, filtering out highest age bin
# ggplot(
#   dplyr::filter(df_summ, age_bin!="61-999"),
#   aes(x=year, y=n_deaths, color=factor(source))
# ) +
#   geom_line() +
#   facet_grid(rows=dplyr::vars(sex), cols=dplyr::vars(age_bin)) +
#   scale_x_continuous(breaks=breaks_year) +
#   labs(color="Source", y="# Deaths total")

# # Plot person-years, filtering out highest age bin
# ggplot(
#   dplyr::filter(df_summ, age_bin!="61-999"),
#   aes(x=year, y=n_py, color=factor(source))
# ) +
#   geom_line() +
#   facet_grid(rows=dplyr::vars(sex), cols=dplyr::vars(age_bin)) +
#   scale_x_continuous(breaks=breaks_year) +
#   labs(color="Source", y="# person-years total")

# Summary stats: person-time
# log_note("Person-time, based on dat_raw$Days", sum(dat_raw$Days)/365)
# log_note("Person-time, based on dat_raw$EndDate and dat_raw$StartDate",
#          sum(as.numeric(dat_raw$EndDate-dat_raw$StartDate)/365.25))
log_note("Person-time, based on df_summ$n_py", sum(df_summ$n_py))
log_note("Person-time, based on nrow(dat_prc)", nrow(dat_prc))


# Summary stats, deaths
# log_note("Deaths, based on dat_raw$Died", sum(dat_raw$Died=="Yes"))
log_note("Deaths, based on dat_prc$y", sum(dat_prc$y))



#######################################################.
##### Generate death rates by age and 45q15 plots #####
#######################################################.

# Prep work
df_summ2 <- data.frame(
  "year" = integer(),
  "sex" = character(),
  "age" = integer(),
  "n_deaths" = integer(),
  "n_py" = double(),
  "rate" = double(),
  "source" = character()
)

# Populate df_summ2
for (year_ in c(cfg$w_start:cfg$w_end)) {
  
  yr_s <- as.Date(x=paste0(year_,"-01-01"))
  yr_e <- as.Date(x=paste0(year_,"-12-31"))
  
  for (sex_ in c("Male", "Female")) {
    for (age_ in c(15:59)) {
      
      # # Summary stats from raw dataset
      # dat_raw$age_at_start <- (
      #   as.Date(x=paste0(year_,"-01-01")) - dat_raw$DoB
      # ) / 365.25
      
      # dat_raw$n_days_exp <- pmax(
      #   pmin(dat_raw$EndDate,yr_e) - pmax(dat_raw$StartDate,yr_s) + 1,
      #   0
      # )
      # n_deaths_raw <- sum(
      #   dat_raw$Died=="Yes" &
      #     as.integer(substr(dat_raw$DoD,1,4))==year_ &
      #     dat_raw$Sex==sex_ &
      #     round(dat_raw$AoD) == age_,
      #   na.rm=T
      # )
      # n_py_raw <- sum(
      #   as.numeric(dat_raw$n_days_exp/365.25) * (
      #     dat_raw$Sex==sex_ &
      #       round(dat_raw$age_at_start) == age_
      #   ),
      #   na.rm=T
      # )
      
      # Summary stats from processed dataset
      dat_prc_filt <- dplyr::filter(
        dat_prc,
        t_end==year_ & sex==as.integer(sex_=="Male") & age == age_
      )
      n_deaths_prc <- sum(dat_prc_filt$y)
      n_py_prc <- nrow(dat_prc_filt)
      
      # # Summary from raw dataset
      # df_summ2[nrow(df_summ2)+1,] <- list(
      #   year_, sex_, age_, n_deaths_raw, n_py_raw,
      #   round(1000*(n_deaths_raw/n_py_raw), 1), "raw"
      # )
      
      # Summary from PRC dataset
      df_summ2[nrow(df_summ2)+1,] <- list(
        year_, sex_, age_, n_deaths_prc, n_py_prc,
        round(1000*(n_deaths_prc/n_py_prc), 1), "Data (prc)"
      )
      
      # Summary from model (age midpoint)
      # Note: if this code is uncommented, the code in sections 1 and 3 of
      #       process.R needs to be run first
      for (status in status_vec) {
        type <- paste0("mort (", status, ")")
        rate <- round(1000 * prob(
          type=type, j=year_, w_1=age_, w_2=cfg$w_2, w_3=as.integer(sex_=="Male"),
          year_start=cfg$w_start, which="est"
        ), 1)
        df_summ2[nrow(df_summ2)+1,] <- list(
          year_, sex_, age_, NA, NA,
          rate, paste0("Model (", status, ")")
        )
      }
      
    }
  }
}

# # Calculate smoothed rates
# for (sex_ in c("Male", "Female")) {
#   for (year_ in c(cfg$w_start:cfg$w_end)) {
# 
#     source_ <- "Data (prc)"
#     df_filt <- dplyr::filter(
#       df_summ2,
#       year==year_ & sex==sex_ & source==source_
#     )
#     spl <- smooth.spline(x=df_filt$age, y=df_filt$rate, df=8)
# 
#     for (i in c(1:length(spl$x))) {
#       age_ <- spl$x[i]
#       df_summ2[nrow(df_summ2)+1,] <- list(
#         year_, sex_, age_, NA, NA, spl$y[i], "Data (smoothed)"
#       )
#     }
#     
#   }
# }

# Plot death rates by age, filtering out highest age bin
# Export: 10 x 5
plot <- ggplot(
  dplyr::filter(df_summ2, year %in% years_plot),
  aes(x=age, y=rate, color=factor(source))
) +
  geom_line() +
  facet_grid(rows=dplyr::vars(sex), cols=dplyr::vars(year)) +
  scale_x_continuous(breaks=seq(10,60,10)) +
  scale_color_manual(values=colors_2) +
  labs(color="Source", y="Deaths per 1,000 person-years")
ggsave(
  filename = paste0("../Figures + Tables/", c_date, " death_rates_by_age - ",
                    "model ", cfg$model_version, ".pdf"),
  plot = plot, device="pdf", width=10, height=5
)

# !!!!! New graph; eventually replace the one above with this one
# Excludes 2010 for conenience to form 3-year age bins
if (T) {
  
  df_summ3 <- df_summ2 %>%
    # dplyr::filter(year!=2010 & source!="Data (smoothed)") %>%
    dplyr::filter(year!=2022 & source!="Data (smoothed)") %>%
    dplyr::mutate(
      year_bin = dplyr::case_when(
        # year %in% seq(2011,2013) ~ "2011 - 2013",
        # year %in% seq(2014,2016) ~ "2014 - 2016",
        # year %in% seq(2017,2019) ~ "2017 - 2019",
        # year %in% seq(2020,2022) ~ "2020 - 2022",
        year %in% seq(2010,2012) ~ "2010 - 2012",
        year %in% seq(2013,2015) ~ "2013 - 2015",
        year %in% seq(2016,2018) ~ "2016 - 2018",
        year %in% seq(2019,2021) ~ "2019 - 2021",
        TRUE ~ NA
      )
    ) %>%
    dplyr::group_by(year_bin, sex, age, source) %>%
    dplyr::summarize(
      n_deaths = sum(n_deaths),
      n_py = sum(n_py),
      rate = mean(rate)
    ) %>% mutate(
      rate = ifelse(source %in% c("Data"),
                    round(1000*(n_deaths/n_py), 1),
                    rate)
    )
  
  plot <- ggplot(
    dplyr::filter(df_summ3),
    aes(x=age, y=rate, color=factor(source))
  ) +
    geom_line() +
    facet_grid(rows=dplyr::vars(sex), cols=dplyr::vars(year_bin)) +
    scale_x_continuous(breaks=seq(10,60,10)) +
    scale_color_manual(values=colors_2) +
    labs(color="Source", y="Deaths per 1,000 person-years")
  ggsave(
    filename = paste0("../Figures + Tables/", c_date, " death_rates_by_age_v2 - ",
                      "model ", cfg$model_version, ".pdf"),
    plot = plot, device="pdf", width=10, height=5
  )
  
}

# Calculate
df_45q15 <- data.frame(
  "year" = integer(),
  "sex" = character(),
  "rate" = double(),
  "source" = character()
)

for (sex_ in c("Male", "Female")) {
  for (year_ in c(cfg$w_start:cfg$w_end)) {
    for (source_ in sources) {
      rate_vec <- dplyr::filter(
        df_summ2,
        year==year_ & sex==sex_ & source==source_
      )$rate
      rate <- 1-prod(1-rate_vec/1000)
      df_45q15[nrow(df_45q15)+1,] <- list(year_, sex_, rate, source_)
    }
  }
}

# Add Thembisa lines
if (cfg$add_thembisa) {
  
  df_45q15_thembisa <- data.frame(
    year = rep(c(2010:2022), 4),
    sex = rep(c("Female", "Male"), each=26),
    rate = c(c(0.2,0.2,0.2,0.2,0.19,0.19,0.19,0.18,0.18,0.18,0.28,0.24,0.2), # Female HIV-
             c(0.77,0.7,0.62,0.59,0.55,0.52,0.49,0.46,0.42,0.39,0.46,0.39,0.37), # Female HIV+
             c(0.38,0.37,0.37,0.36,0.35,0.35,0.34,0.33,0.33,0.31,0.41,0.39,0.33), # Male HIV-
             c(0.89,0.85,0.80,0.78,0.77,0.75,0.73,0.71,0.7,0.68,0.71,0.68,0.67)), # Male HIV+
    source = rep(rep(c("Thembisa (HIV-)", "Thembisa (HIV+)"), 2), each=13)
  )
  df_45q15 <- rbind(df_45q15, df_45q15_thembisa)
  
  # Plot 45q15 death rates
  plot <- ggplot(
    df_45q15,
    aes(x=year, y=rate, color=factor(source), linetype=factor(source))
  ) +
    geom_line() +
    facet_grid(cols=dplyr::vars(sex)) +
    scale_x_continuous(breaks=breaks_year) +
    # scale_color_manual(values=c(colors_1,colors_1[1:2])) +
    scale_color_manual(values=c(colors_1,colors_1[2:3])) +
    # scale_linetype_manual(values=c(rep("solid",4), rep("dashed",2))) +
    scale_linetype_manual(values=c(rep("solid",3), rep("dashed",2))) +
    labs(color="Source", linetype="Source",
         y="Probability of death between ages 15-60")
  
} else {
  
  # Plot 45q15 death rates
  plot <- ggplot(
    df_45q15,
    aes(x=year, y=rate, color=factor(source))
  ) +
    geom_line() +
    facet_grid(cols=dplyr::vars(sex)) +
    scale_x_continuous(breaks=breaks_year) +
    scale_color_manual(values=colors_1) +
    labs(color="Source", y="Probability of death between ages 15-60")
  
}

ggsave(
  filename = paste0("../Figures + Tables/", c_date, " 45q15 - ",
                    "model ", cfg$model_version, ".pdf"),
  plot = plot, device="pdf", width=6, height=4
)
