# Title: "Modeling HIV seroconversion dates"
# Author: Avi Kenny



##################.
##### CONFIG #####
##################.

# Set global config
cfg <- list(
  level_set_which = "level_set_1",
  # run_or_update = "run",
  num_sim = 500,
  pkgs = c("dplyr", "survival", "data.table", "tidyr", "memoise", "Rsolnp",
           "numDeriv"), # "rjags", "rstan"
  pkgs_nocluster = c("ggplot2"),
  parallel = "none", # none outer
  stop_at_error = F
  # mcmc = list(n.adapt=1000, n.iter=1000, n.burn=1000, n.chains=2, thin=1)
)

# Set cluster config
cluster_config <- list(
  js = "slurm",
  dir = "/home/akenny/z.hivmi"
)



#################.
##### Setup #####
#################.

# Set local vs. cluster variables
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
  # Local
  setwd(paste0("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Resear",
               "ch/Mark Siedner/Project - HIVMI/HIV.multiple.imputation/R"))
  cfg$local <- TRUE
} else {
  # Cluster
  setwd("z.hivmi/R")
  cfg$local <- FALSE
}

# Load packages (if running locally)
if (cfg$local) {
  for (pkg in c(cfg$pkgs,cfg$pkgs_nocluster)) {
    suppressMessages({ do.call("library", list(pkg)) })
  }
}

# Load SimEngine + functions
{
  library(SimEngine)
  source("one_simulation.R")
  source("generate_data.R")
  source("likelihood_miss.R")
  source("likelihood_miss_nocovariates.R")
  source("helpers.R")
}



##########################################################.
##### MAIN: Set level sets for different simulations #####
##########################################################.

if (Sys.getenv("sim_run") %in% c("first", "")) {
  
  # L <- list(n=300,max_time=70,params=list(a_x=log(0.005),a_y=log(0.003),a_v=log(0.7),g_x=c(log(1.3),log(1.2)),g_y=c(log(1.2),log(1.1)),g_v=c(log(1.2),log(1.1)),beta=log(1.5)))
  
  # Simulation 1: basic
  # n=500,t=100 rep runs in 3.2 hrs # !!!!! outdated
  level_set_1 <- list(
    n = 500,
    # n = c(500,1000,2000),
    max_time = 70,
    # max_time = 100,
    params = list(
      # "10pct testing" = list(
      #   a_x=log(0.005), a_y=log(0.003), a_v=log(0.1),
      #   g_x=c(log(1.3),log(1.2)), g_y=c(log(1.2),log(1.1)),
      #   g_v=c(log(1.2),log(1.1)), beta=log(1.5)
      # ),
      # "70pct testing" = list(
      #   a_x=log(0.005), a_y=log(0.003), a_v=log(0.7),
      #   g_x=c(log(1.3),log(1.2)), g_y=c(log(1.2),log(1.1)),
      #   g_v=c(log(1.2),log(1.1)), beta=log(1.5)
      # )
      
      
      
      "70pct testing" = list( # !!!!! 2022-09-28
        a_x=log(0.005), a_y=log(0.003), a_v=log(0.7), # !!!!! 2022-09-28
        g_x=c(0,0), g_y=c(log(1.5),log(1.1)), # !!!!! 2022-09-28
        g_v=c(log(1.2),log(1.1)), beta=log(1.5) # !!!!! 2022-09-28
      ) # !!!!! 2022-09-28
    )
  )
  
  # Simulation 2: no covariates
  # n=500,t=100 rep runs in 3.2 hrs # !!!!! outdated
  level_set_2 <- list(
    n = 500,
    max_time = 70,
    params = list(
      "70pct testing" = list(
        a_x=log(0.005), a_y=log(0.003), a_v=log(0.7), g_x=c(0,0),
        g_y=c(0,0), g_v=c(0,0), beta=log(1.5)
      )
    )
  )
  
  level_set <- get(cfg$level_set_which)
  
}



#################################.
##### MAIN: Simulation code #####
#################################.

# Commands for job sumbission on Slurm:
# sbatch --export=sim_run='first',cluster='bionic',type='R',project='z.hivmi' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:11 --array=1-16 --export=sim_run='main',cluster='bionic',type='R',project='z.hivmi' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:12 --export=sim_run='last',cluster='bionic',type='R',project='z.hivmi' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh

run_on_cluster(
  
  first = {
    
    # Set up and configure simulation object
    sim <- new_sim()
    sim %<>% set_config(
      num_sim = cfg$num_sim,
      parallel = cfg$parallel,
      stop_at_error = cfg$stop_at_error,
      packages = cfg$pkgs
    )
    
    # Set simulation script
    sim %<>% set_script(one_simulation)
    
    # Set levels
    sim <- do.call(set_levels, c(list(sim), level_set))
    
  },
  
  main = { sim %<>% run() },
  
  last = { sim %>% summarize() %>% print() },
  
  cluster_config = cluster_config
  
)



#######################################.
##### VIZ: All params (one model) #####
#######################################.

if (F) {
  
  # True vals: a_x=log(0.005), a_y=log(0.003), a_v=log(0.7), g_x=c(log(1.3),log(1.2)), g_y=c(log(1.2),log(1.1)), g_v=c(log(1.2),log(1.1)), beta=log(1.5)
  
  # v <- c("lik_M_g_y1_est", "lik_M_beta_est")
  # true_vals <- log(c(1.5,1.5))
  v <- c("lik_M_a_x_est", "lik_M_g_x1_est", "lik_M_g_x2_est", "lik_M_a_y_est",
         "lik_M_g_y1_est", "lik_M_g_y2_est", "lik_M_beta_est")
  # true_vals <- log(c(0.005,1.3,1.2,0.003,1.5,1.1,1.5))
  true_vals <- log(c(0.005,1,1,0.003,1.5,1.1,1.5))
  r <- filter(sim$results, params=="70pct testing")
  x <- unlist(lapply(v, function(col) { r[,col] }))
  df_true <- data.frame(
    which = factor(v, levels=v),
    val = true_vals
  )
  df_plot <- data.frame(
    x = x,
    y = rep(0, length(x)),
    which = rep(factor(v, levels=v), each=round(length(x)/length(v)))
  )
  
  # Export 8" x 5"
  ggplot(df_plot, aes(x=x, y=y, color=which)) +
    geom_jitter(width=0, height=1, alpha=0.3, size=3) +
    geom_vline(aes(xintercept=val), df_true, alpha=0.5, linetype="dashed") +
    facet_wrap(~which, ncol=1, strip.position="left") + # scales="free_x"
    labs(y=NULL, title="70% testing") +
    ylim(-2,2) +
    # xlim(-3,3) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text.y.left = element_text(angle=0),
          legend.position="none")
  
  # Coverage
  summ <- sim %>% summarize(
    mean = list(
      list(name="sd_a_x_est", x="lik_M_a_x_se"),
      list(name="sd_g_x1_est", x="lik_M_g_x1_se"),
      list(name="sd_g_x2_est", x="lik_M_g_x2_se"),
      list(name="sd_a_y_est", x="lik_M_a_y_se"),
      list(name="sd_g_y1_est", x="lik_M_g_y1_se"),
      list(name="sd_g_y2_est", x="lik_M_g_y2_se"),
      list(name="sd_beta_est", x="lik_M_beta_se")
    ),
    sd = list(
      list(name="sd_a_x_actual", x="lik_M_a_x_est"),
      list(name="sd_g_x1_actual", x="lik_M_g_x1_est"),
      list(name="sd_g_x2_actual", x="lik_M_g_x2_est"),
      list(name="sd_a_y_actual", x="lik_M_a_y_est"),
      list(name="sd_g_y1_actual", x="lik_M_g_y1_est"),
      list(name="sd_g_y2_actual", x="lik_M_g_y2_est"),
      list(name="sd_beta_actual", x="lik_M_beta_est")
    ),
    coverage = list(
      list(name="cov_a_x", truth=log(0.005), estimate="lik_M_a_x_est", se="lik_M_a_x_se", na.rm=T),
      list(name="cov_g_x1", truth=log(1), estimate="lik_M_g_x1_est", se="lik_M_g_x1_se", na.rm=T),
      list(name="cov_g_x2", truth=log(1), estimate="lik_M_g_x2_est", se="lik_M_g_x2_se", na.rm=T),
      list(name="cov_a_y", truth=log(0.003), estimate="lik_M_a_y_est", se="lik_M_a_y_se", na.rm=T),
      list(name="cov_g_y1", truth=log(1.5), estimate="lik_M_g_y1_est", se="lik_M_g_y1_se", na.rm=T),
      list(name="cov_g_y2", truth=log(1.1), estimate="lik_M_g_y2_est", se="lik_M_g_y2_se", na.rm=T),
      list(name="cov_beta", truth=log(1.5), estimate="lik_M_beta_est", se="lik_M_beta_se", na.rm=T)
    )
  )
  
}



##############################.
##### VIZ: Others (temp) #####
##############################.

if (F) {
  
  
  
  # !!!!!
  {
    # True params: a_x=0.005, a_y=0.003, a_v=0.1/0.7, g_x=c(1.3,1.002), g_y=c(1.2,1.001), g_v=c(1.2,1.001), beta=1.5
    true_val <- 1.5
    x_window <- 0.2
    # v1 <- "cox_g_y2_est"
    # v2 <- "lik_F_g_y2_est"
    # v3 <- "lik_M_g_y2_est"
    v1 <- "cox_beta_est"
    v2 <- "lik_F_beta_est"
    v3 <- "lik_M_beta_est"
    r1 <- filter(sim$results, params=="10pct testing")
    r2 <- filter(sim$results, params=="70pct testing")
    # r1 <- filter(sim$results, n==1000 & max_time==100)
    x <- c(r1[[v1]], r1[[v2]], r1[[v3]],
           r2[[v1]], r2[[v2]], r2[[v3]])
    stats <- c("10pct, Cox", "10pct, Lik F", "10pct, Lik M",
               "70pct, Cox", "70pct, Lik F", "70pct, Lik M")
    df_plot <- data.frame(
      x = x,
      y = rep(0, length(x)),
      which = rep(factor(stats, levels=stats), each=round(length(x)/length(stats)))
    )
    ggplot(df_plot, aes(x=x, y=y, color=which)) +
      geom_vline(xintercept=true_val, alpha=0.5, linetype="dashed") +
      geom_jitter(width=0, height=1, alpha=0.3, size=3) +
      facet_wrap(~which, ncol=1, strip.position="left") + # scales="free_x"
      labs(y=NULL) +
      ylim(-2,2) +
      xlim(true_val-x_window,true_val+x_window) +
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            strip.text.y.left = element_text(angle=0),
            legend.position="none")
    
  }
  
  # r <- sim$results
  r500 <- filter(sim$results, n==500)
  r1000 <- filter(sim$results, n==1000)
  r2000 <- filter(sim$results, n==2000)
  ln <- c(nrow(r500), nrow(r1000), nrow(r2000))
  plot_data <- data.frame(
    x = c(r500$lik_beta_est, r1000$lik_beta_est, r2000$lik_beta_est),
    y = c(r500$cox_beta_est, r1000$cox_beta_est, r2000$cox_beta_est),
    n = c(rep(500, ln[1]), rep(1000, ln[2]), rep(2000, ln[3]))
  )
  ggplot(plot_data, aes(x=x, y=y)) +
    geom_abline(slope=1, intercept=0, color="orange") +
    geom_point(alpha=0.3) +
    geom_point(data=data.frame(x=1.5,y=1.5), color="forestgreen", size=3, alpha=0.7) +
    facet_wrap(~n) +
    labs(x="Likelihood model", y="Cox model")
  
  sim %>% summarize(
    coverage = list(
      list(name="cov_cox", truth=1.5, lower="cox_beta_ci_lo", upper="cox_beta_ci_hi"),
      list(name="cov_lik", truth=1.5, lower="lik_beta_ci_lo", upper="lik_beta_ci_hi")
    )
  )

  
  
  # # Read in simulation object
  # sim <- readRDS("../simba.out/sim_20210615.simba")
  
  # # Transform results
  # sim$results %<>% mutate(
  #   # est_hr_hiv = exp(est_hiv),
  #   # est_hr_art = exp(est_art),
  #   hr_hiv_lab = paste("HIV+ HR:",hr_hiv),
  #   hr_art_lab = paste("ART+ HR:",hr_art),
  #   method = ifelse(method=="ideal","Ideal",ifelse(method=="mi","MI","")),
  #   log_hr_hiv = log(hr_hiv),
  #   log_hr_art = log(hr_art)
  # )
  
  # # Plot of estimates (HIV)
  # # Export: 6" X 3"
  # ggplot(sim$results, aes(x=method, y=est_hiv, color=method)) +
  #   geom_point(alpha=0.2, size=2) +
  #   geom_hline(aes(yintercept=log(hr_hiv)), linetype="dotted") +
  #   facet_wrap(~hr_hiv_lab, ncol=2) +
  #   theme(legend.position="none")
  
  # # Plot of estimates (ART)
  # # Export: 6" X 3"
  # ggplot(sim$results, aes(x=method, y=est_art, color=method)) +
  #   geom_point(alpha=0.2, size=2) +
  #   geom_hline(aes(yintercept=log(hr_art)), linetype="dotted") +
  #   facet_wrap(~hr_art_lab, ncol=4) +
  #   theme(legend.position="none")
  
  
  
  
  
  # # HIV graph
  # # Export PDF 3x8
  # ggplot(sim$results, aes(x=method, y=est_hr_hiv, color=method)) +
  #   geom_point(alpha=0.2, size=2) +
  #   geom_hline(aes(yintercept=hr_hiv), linetype="dotted") +
  #   facet_wrap(~hr_hiv_lab, ncol=4) +
  #   theme(legend.position="none") +
  #   labs(title="Point estimates (50 simulation replicates per level)",
  #        x="Method", y="Estimated hazard ratio (HIV+ART-)")
  
  # # Coverage plots
  # summ <- sim %>% summary(
  #   coverage = list(
  #     name="cov_hiv", estimate="est_hiv", truth="log_hr_hiv", se="se_hiv"
  #   )
  # ) %>% mutate(
  #   hr_hiv_lab = paste("HIV+ HR:",hr_hiv),
  #   hr_art_lab = paste("ART+ HR:",hr_art),
  # )
  
  # ggplot(summ, aes(x=method, y=cov_hiv, color=method)) +
  #   geom_point(size=3) +
  #   geom_hline(aes(yintercept=0.95), linetype="dotted") +
  #   facet_wrap(~hr_hiv_lab, ncol=4) +
  #   theme(legend.position="none") +
  #   labs(title="Coverage (50 simulation replicates per level)",
  #        x="Method", y="95% CI Coverage (HIV+ART-)")
  
}



##############################.
##### Run dataset checks #####
##############################.

if (F) {
  
  # Setup
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(magrittr)
  library(data.table)
  library(survival)
  source("generate_dataset.R")
  source("perform_imputation.R")
  source("transform_dataset.R")
  source("run_analysis.R")
  source("helpers.R")
  p_sero_year <- convert_p_sero(list(male = list("1"=0.03, "2-10"=0, "11-15"=0, "16-20"=0.01, "21-25"=0.02,"26-30"=0.03, "31-35"=0.02, "36-40"=0.01, "41-45"=0.005, "46-50"=0.005),female = list("1"=0.03, "2-10"=0, "11-15"=0.005, "16-20"=0.02, "21-25"=0.03,"26-30"=0.02, "31-35"=0.015, "36-40"=0.01, "41-45"=0.005, "46-50"=0)))
  
  # Create "original" dataset
  dat_orig <- generate_dataset(
    num_patients = 10000,
    start_year = 2000,
    end_year = 2020,
    hazard_ratios = list("hiv"=1.7, "art"=1.4),
    p_sero_year = p_sero_year,
    p_death_year = p_death_year(1),
    u_mult = list("sero"=1, "death"=1)
  )
  
  # Create "imputed" dataset
  dat_imp <- perform_imputation(dat_orig, p_sero_year)
  
  # Create transformed datasets
  dat_cp_orig <- transform_dataset(dat_orig)
  dat_cp_imp <- transform_dataset(dat_imp)
  
  # !!!!! MICE TESTING: START
  # !!!!! TO DO
  # !!!!! MICE TESTING: END
  
  # Check 1: Compare overall seroconversion rates
  print(1-sum(is.na(dat_orig$sero_year))/nrow(dat_orig))
  print(1-sum(is.na(dat_imp$sero_year))/nrow(dat_imp))
  dat_orig %>% filter(sero_year<2000) %>% nrow()
  dat_imp %>% filter(sero_year<2000) %>% nrow()
  dat_orig %>% filter(sero_year>=2000) %>% nrow()
  dat_imp %>% filter(sero_year>=2000) %>% nrow()
  
  # Check 2: Check to see if distribution of seroconversion years is similar
  # !!!!! orig is getting much higher pre-2000 seroconversion rate
  ggplot(data.frame(
    sero_year = c(dat_orig$sero_year,dat_imp$sero_year),
    which = rep(c("original","imputed"),each=nrow(dat_orig))
  ), aes(x=sero_year, group=which, fill=factor(which))) +
    geom_histogram(color="white", bins=100) +
    geom_vline(xintercept = 2000, linetype="dotted") +
    facet_wrap(~which, ncol=2)
  
  # Check 3: Look at seroconversion years by "testing case"
  d3_orig <- dat_orig %>% group_by(case, sero_year) %>% summarize(count=n())
  d3_imp <- dat_imp %>% group_by(case, sero_year) %>% summarize(count=n())
  d3_orig$which <- "orig"
  d3_imp$which <- "imp"
  ggplot(
    rbind(d3_orig,d3_imp),
    aes(x=sero_year, y=count, group=which, fill=factor(case))
  ) +
    geom_bar(stat="identity", width=0.4) +
    geom_vline(xintercept = 2000, linetype="dotted") +
    facet_wrap(~case+which, ncol=2)
  
  # Check 4: Plot cascade status over time
  d4_orig <- dat_cp_orig %>% mutate(
    case = case_when(
      case==1 ~ "1. No testing data",
      case==2 ~ "2. Last test was neg",
      case==3 ~ "3. Neg test then pos test",
      case==4 ~ "4. First test was pos"
    )
  ) %>% group_by(start_year, casc_status, case) %>% summarize(num=n())
  ggplot(d4_orig, aes(x=start_year, y=num, color=casc_status)) +
    geom_line() + facet_wrap(~case, ncol=2)
  d4_imp <- dat_cp_imp %>% mutate(
    case = case_when(
      case==1 ~ "1. No testing data",
      case==2 ~ "2. Last test was neg",
      case==3 ~ "3. Neg test then pos test",
      case==4 ~ "4. First test was pos"
    )
  ) %>% group_by(start_year, casc_status, case) %>% summarize(num=n())
  ggplot(d4_imp, aes(x=start_year, y=num, color=casc_status)) +
    geom_line() + facet_wrap(~case, ncol=2)
  
  # Check 5: examine cascade status by case
  d5 <- dat_cp_orig %>% group_by(case, casc_status) %>%
    summarize(count=n())
  d5_imp <- dat_cp_imp %>% group_by(case, casc_status) %>% summarize(count=n())
  d5$imp <- d5_imp$count
  d5 %<>% rename(orig="count")
  d5b <- d5 %>% group_by(casc_status) %>% summarize(orig=sum(orig),imp=sum(imp))
  print(d5)
  print(d5b)
  print(d5 %>% filter(casc_status=="HIV+ART-") %>% mutate(diff=imp-orig))
  
  # Check 6: run analysis on both datasets; should yield comparable SEs
  run_analysis(dat_cp_orig, method="ideal")
  run_analysis(dat_cp_imp, method="ideal")
  
  # !!!!!!
  run_analysis(dat_cp_orig, method="censor")
  run_analysis(dat_cp_imp, method="censor")
  
  # Check 7: Examine case counts
  xtabs(~case, data=dat_orig)

  # Check 8: Examine death counts
  xtabs(~died, data=dat_orig)
  
  # Check 9: Compare death counts by casc_status
  # !!!!! This illustrates the source of the bias
  xtabs(~casc_status, data=filter(dat_cp_orig, died==1))
  xtabs(~casc_status, data=filter(dat_cp_imp, died==1))
  
  # Check 10: Compare death counts by case+casc_status
  xtabs(~case+casc_status, data=filter(dat_cp_orig, died==1))
  xtabs(~case+casc_status, data=filter(dat_cp_imp, died==1))
  
  # Check 11: examine datasets manually for abnormalities
  # write.table(dat_orig, file="dat_orig.csv", sep=",", row.names=FALSE)
  # write.table(dat_imp, file="dat_imp.csv", sep=",", row.names=FALSE)
  # write.table(dat_cp_orig, file="dat_cp_orig.csv", sep=",", row.names=FALSE)
  # write.table(dat_cp_imp, file="dat_cp_imp.csv", sep=",", row.names=FALSE)
  
}



###########################.
##### Mini-simulation #####
###########################.

if (F) {
  
  n <- 100000
  
  # x is our binary exposure variable
  prob_x <- 0.2
  x <- rbinom(n=n, size=1, prob=prob_x)
  
  # y is our binary outcome, correlated with x
  # rho_xy, the correlation between x and y, is the outcome of interest
  rho_xy <- 0.8
  y <- ifelse(runif(n)<rho_xy, x, rbinom(n=n, size=1, prob=prob_x))
  print(cor(x,y))
  
  est_corr <- c()
  rho_zx_vec <- seq(0,1,0.1)
  for (rho_zx in rho_zx_vec) {
    
    # z is a covariate correlated directly with x
    # Within this loop, we test different values of rho_zx, the correlation
    #     between x and z
    z <- ifelse(runif(n)<rho_zx, x, rbinom(n=n, size=1, prob=prob_x))
    
    # Impose missingness by deleting (100*k)% of the x values
    k <- 0.8
    x_trunc <- c(x[round(1:(n*(1-k)))],rep(NA,(n*k)))
    
    # Perform multiple imputation of x based on z; then estimate rho_xy
    est_cor_mi <- c()
    for (m in 1:10) {
      z_miss <- z[round((n*(1-k)+1):n)]
      x_imp <- ifelse(runif(n*k)<rho_zx, z_miss,
                      rbinom(n=(n*k), size=1, prob=prob_x))
      x_new <- c(x[round(1:(n*(1-k)))],x_imp)
      est_cor_mi <- c(est_cor_mi,cor(x_new,y))
    }
    est_corr <- c(est_corr, mean(est_cor_mi))

  }
  
  library(ggplot2)
  ggplot(
    data.frame(x=rho_zx_vec,y=est_corr),
    aes(x=x,y=y)) +
    geom_point() +
    geom_hline(yintercept=0.8, linetype="dotted") +
    labs(x="Correlation between x and z", y="Estimated rho_xy")
  
}
