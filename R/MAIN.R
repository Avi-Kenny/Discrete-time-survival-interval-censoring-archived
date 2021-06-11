# Title: "Modeling HIV seroconversion dates"
# Author: Avi Kenny



##################.
##### CONFIG #####
##################.

# Set global config
cfg <- list(
  level_set_which = "level_set_1",
  # run_or_update = "run",
  num_sim = 1,
  pkgs = c("dplyr", "survival", "data.table", "tidyr", "rjags", "rstan",
           "tidyr"),
  pkgs_nocluster = c(),
  parallel = "none",
  stop_at_error = FALSE,
  mcmc = list(n.adapt=1000, n.iter=1000, n.burn=1000, n.chains=2, thin=1),
  local = FALSE
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
}

# Load packages (if running locally)
if (cfg$local) {
  for (pkg in c(cfg$pkgs,cfg$pkgs_nocluster)) {
    do.call("library", list(pkg))
  }
}

# Load simba + functions
library(simba) # devtools::install_github(repo="Avi-Kenny/simba")
source("generate_data_baseline.R")
source("generate_data_events.R")
source("run_analysis.R")
source("perform_imputation.R")
source("fit_jags.R")
source("transform_jags.R")
source("transform_dataset.R")
source("one_simulation.R")
source("helpers.R")



##########################################################.
##### MAIN: Set level sets for different simulations #####
##########################################################.

if (Sys.getenv("run") %in% c("first", "")) {
  
  # Simulation 1: !!!!!
  level_set_1 <- list(
    method = c("ideal", "mi"), # c("ideal", "mi", "censor")
    hr_hiv = c(0.6,1.0,1.4,1.8),
    hr_art = 0.7,
    p_death_mult = 1,
    u_mult_sero = 1,
    u_mult_death = 1
  )
  
  level_set <- eval(as.name(cfg$level_set_which))
  
}



#################################.
##### MAIN: Simulation code #####
#################################.

# Commands for job sumbission on Slurm:
# sbatch --export=run='first',cluster='bionic',type='R',project='z.hivmi' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:11 --array=1-16 --export=run='main',cluster='bionic',type='R',project='z.hivmi' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:12 --export=run='last',cluster='bionic',type='R',project='z.hivmi' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh

run_on_cluster(
  
  first = {
    
    # Set up and configure simulation object
    sim <- new_sim()
    sim %<>% set_config(
      num_sim = cfg$num_sim, # !!!!!
      parallel = cfg$parallel,
      stop_at_error = cfg$stop_at_error,
      packages = cfg$pkgs
    )
    
    # Specify seroconversion conditional probabilities (discrete hazards)
    # For an individual of (exactly) age X-1, the number in the corresponding
    #     age bin represents the probability that the individual will
    #     seroconvert sometime in their Xth full year of life, given that they
    #     did not seroconvert by their (X-1)th full year of life.
    # !!!!! Change the actual numbers later based on AHRI cohort data
    p_sero_year <- convert_p_sero(list(
      male = list(
        "1"=0.03, "2-10"=0, "11-15"=0, "16-20"=0.01, "21-25"=0.02,
        "26-30"=0.03, "31-35"=0.02, "36-40"=0.01, "41-45"=0.005, "46-50"=0.005
      ),
      female = list(
        "1"=0.03, "2-10"=0, "11-15"=0.005, "16-20"=0.02, "21-25"=0.03,
        "26-30"=0.02, "31-35"=0.015, "36-40"=0.01, "41-45"=0.005, "46-50"=0
      )
    ))
    
    sim %<>% add_constants(
      p_sero_year = p_sero_year,
      m = 5, # Number of MI replicates
      num_patients = 5000,
      start_year = 2000,
      end_year = 2020
    )
    
    # Add functions to simulation object
    sim %<>% add_creator(generate_dataset)
    sim %<>% set_script(one_simulation)
    sim %<>% add_method(expit)
    sim %<>% add_method(convert_p_sero)
    sim %<>% add_method(construct_m_probs)
    sim %<>% add_method(run_analysis)
    sim %<>% add_method(perform_imputation)
    sim %<>% add_method(transform_dataset)
    sim %<>% add_method(p_death_year)
    
    # Set levels
    sim <- do.call(set_levels, c(list(sim), level_set))
    
  },
  
  main = {
    sim %<>% run()
  },
  
  last = {
    
    sim %>% summary() %>% print()
    
  },
  
  cluster_config = cluster_config
  
)



###########################.
##### Process results #####
###########################.
if (FALSE) {
  
  library(ggplot2)
  
  # !!!!! Calculate bias, variance inflation, power
  

  sim$results %<>% mutate(
    est_hr_hiv = exp(est_hiv),
    est_hr_art = exp(est_art),
    hr_hiv_lab = paste("HIV+ HR:",hr_hiv),
    hr_art_lab = paste("ART+ HR:",hr_art),
    method = ifelse(method=="ideal", "Ideal",
                    ifelse(method=="mi","MI","error")),
    log_hr_hiv = log(hr_hiv),
    log_hr_art = log(hr_art)
  )
  
  # HIV graph
  # Export PDF 3x8
  ggplot(sim$results, aes(x=method, y=est_hr_hiv, color=method)) +
    geom_point(alpha=0.2, size=2) +
    geom_hline(aes(yintercept=hr_hiv), linetype="dotted") +
    facet_wrap(~hr_hiv_lab, ncol=4) +
    theme(legend.position="none") +
    labs(title="Point estimates (50 simulation replicates per level)",
         x="Method", y="Estimated hazard ratio (HIV+ART-)")
  
  # Coverage plots
  summ <- sim %>% summary(
    coverage = list(
      name="cov_hiv", estimate="est_hiv", truth="log_hr_hiv", se="se_hiv"
    )
  ) %>% mutate(
    hr_hiv_lab = paste("HIV+ HR:",hr_hiv),
    hr_art_lab = paste("ART+ HR:",hr_art),
  )
  
  ggplot(summ, aes(x=method, y=cov_hiv, color=method)) +
    geom_point(size=3) +
    geom_hline(aes(yintercept=0.95), linetype="dotted") +
    facet_wrap(~hr_hiv_lab, ncol=4) +
    theme(legend.position="none") +
    labs(title="Coverage (50 simulation replicates per level)",
         x="Method", y="95% CI Coverage (HIV+ART-)")
  
}



##############################.
##### Run dataset checks #####
##############################.

if (FALSE) {
  
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
  dataset_orig <- generate_dataset(
    num_patients = 10000,
    start_year = 2000,
    end_year = 2020,
    hazard_ratios = list("hiv"=1.7, "art"=1.4),
    p_sero_year = p_sero_year,
    p_death_year = p_death_year(1),
    u_mult = list("sero"=1, "death"=1)
  )
  
  # Create "imputed" dataset
  dataset_imp <- perform_imputation(dataset_orig, p_sero_year)
  
  # Create transformed datasets
  dataset_cp_orig <- transform_dataset(dataset_orig)
  dataset_cp_imp <- transform_dataset(dataset_imp)
  
  # !!!!! MICE TESTING: START
  # !!!!! TO DO
  # !!!!! MICE TESTING: END
  
  # Check 1: Compare overall seroconversion rates
  print(1-sum(is.na(dataset_orig$sero_year))/nrow(dataset_orig))
  print(1-sum(is.na(dataset_imp$sero_year))/nrow(dataset_imp))
  dataset_orig %>% filter(sero_year<2000) %>% nrow()
  dataset_imp %>% filter(sero_year<2000) %>% nrow()
  dataset_orig %>% filter(sero_year>=2000) %>% nrow()
  dataset_imp %>% filter(sero_year>=2000) %>% nrow()
  
  # Check 2: Check to see if distribution of seroconversion years is similar
  # !!!!! orig is getting much higher pre-2000 seroconversion rate
  ggplot(data.frame(
    sero_year = c(dataset_orig$sero_year,dataset_imp$sero_year),
    which = rep(c("original","imputed"),each=nrow(dataset_orig))
  ), aes(x=sero_year, group=which, fill=factor(which))) +
    geom_histogram(color="white", bins=100) +
    geom_vline(xintercept = 2000, linetype="dotted") +
    facet_wrap(~which, ncol=2)
  
  # Check 3: Look at seroconversion years by "testing case"
  d3_orig <- dataset_orig %>% group_by(case, sero_year) %>% summarize(count=n())
  d3_imp <- dataset_imp %>% group_by(case, sero_year) %>% summarize(count=n())
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
  d4_orig <- dataset_cp_orig %>% mutate(
    case = case_when(
      case==1 ~ "1. No testing data",
      case==2 ~ "2. Last test was neg",
      case==3 ~ "3. Neg test then pos test",
      case==4 ~ "4. First test was pos"
    )
  ) %>% group_by(start_year, casc_status, case) %>% summarize(num=n())
  ggplot(d4_orig, aes(x=start_year, y=num, color=casc_status)) +
    geom_line() + facet_wrap(~case, ncol=2)
  d4_imp <- dataset_cp_imp %>% mutate(
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
  d5 <- dataset_cp_orig %>% group_by(case, casc_status) %>%
    summarize(count=n())
  d5_imp <- dataset_cp_imp %>% group_by(case, casc_status) %>% summarize(count=n())
  d5$imp <- d5_imp$count
  d5 %<>% rename(orig="count")
  d5b <- d5 %>% group_by(casc_status) %>% summarize(orig=sum(orig),imp=sum(imp))
  print(d5)
  print(d5b)
  print(d5 %>% filter(casc_status=="HIV+ART-") %>% mutate(diff=imp-orig))
  
  # Check 6: run analysis on both datasets; should yield comparable SEs
  run_analysis(dataset_cp_orig, method="ideal")
  run_analysis(dataset_cp_imp, method="ideal")
  
  # !!!!!!
  run_analysis(dataset_cp_orig, method="censor")
  run_analysis(dataset_cp_imp, method="censor")
  
  # Check 7: Examine case counts
  xtabs(~case, data=dataset_orig)

  # Check 8: Examine death counts
  xtabs(~died, data=dataset_orig)
  
  # Check 9: Compare death counts by casc_status
  # !!!!! This illustrates the source of the bias
  xtabs(~casc_status, data=filter(dataset_cp_orig, died==1))
  xtabs(~casc_status, data=filter(dataset_cp_imp, died==1))
  
  # Check 10: Compare death counts by case+casc_status
  xtabs(~case+casc_status, data=filter(dataset_cp_orig, died==1))
  xtabs(~case+casc_status, data=filter(dataset_cp_imp, died==1))
  
  # Check 11: examine datasets manually for abnormalities
  # write.table(dataset_orig, file="dataset_orig.csv", sep=",", row.names=FALSE)
  # write.table(dataset_imp, file="dataset_imp.csv", sep=",", row.names=FALSE)
  # write.table(dataset_cp_orig, file="dataset_cp_orig.csv", sep=",", row.names=FALSE)
  # write.table(dataset_cp_imp, file="dataset_cp_imp.csv", sep=",", row.names=FALSE)
  
}



###########################.
##### Mini-simulation #####
###########################.

if (FALSE) {
  
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



###################.
##### Archive #####
###################.

if (FALSE) {
  
  # Misc MICE code
  {
    # Run analysis on each imputed dataset
    for (j in 1:m) {
      d_imputed <- mice::complete(imputation_object, j)
    }
    
    # Get degrees of freedom from a single analysis
    dfcom <- summary(cox_model_ideal)$waldtest[[2]]
    
    # Store list as mice::mira object and pool results
    imputation_results <- as.mira(analysis_list)
    pooled_model <- pool(imputation_results, dfcom = dfcom)
    
    # Get coefficient of hiv_status(3)
    coeff_ideal <- summary(cox_model_ideal)$coefficients[2,1]
    coeff_mi <- pooled_model$pooled[2,1]
  }
  
  # Convert yearly probabilities to monthly probabilities
  {
    convert_to_monthly_prob <- function(p) { 1 - (1-p)^(1/12) }
    psero <- list(
      mtct = psero_year$mtct,
      male = lapply(psero_year$male, convert_to_monthly_prob),
      female = lapply(psero_year$female, convert_to_monthly_prob)
    )
  }
  
  # Dataset checks
  {
    
    # Check 4: Look at cascade status by year and "testing case"
    d4_orig <- dataset_cp_orig %>% group_by(case, start_year) %>%
      summarize(num_hivpos_artneg=sum(casc_status=="HIV+ART-"))  
    d4_imp <- dataset_cp_imp %>% group_by(case, start_year) %>%
      summarize(num_hivpos_artneg=sum(casc_status=="HIV+ART-"))  
    d4_orig$which <- "orig"
    d4_imp$which <- "imp"
    ggplot(
      rbind(d4_orig,d4_imp),
      aes(x=start_year, y=num_hivpos_artneg, group=which, fill=factor(case))
    ) +
      geom_bar(stat="identity", width=0.4) +
      facet_grid(rows=vars(case), cols=vars(which), scales="free_y")
    
  }
  
}
