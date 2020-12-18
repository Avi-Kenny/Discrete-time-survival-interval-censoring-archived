# Title: "Modeling HIV seroconversion dates"
# Author: Avi Kenny
# Date: 2020-12-17



#################.
##### Setup #####
#################.

# Set working directory
# If running this locally in RStudio, skip this section and click "Session" >>
#     "Set working directory" >> "To Source File Location"
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
  setwd("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Research/Mark Siedner/Project - HIVMI/HIV.multiple.imputation/R")
} else {
  setwd("z.hivmi/R")
}

# Set code blocks to run
{
  run_main <- TRUE
  run_results <- FALSE
  run_checks <- FALSE
}



###########################.
##### Simulation code #####
###########################.

if (run_main) {
  
  library(simba)
  
  # Commands for submission on Bionic cluster:
  # sbatch --export=run='first',cluster='bionic',type='R',project='z.hivmi' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out run_r.sh
  # sbatch --depend=afterok:101 --array=1-32 --export=run='main',cluster='bionic',type='R',project='z.hivmi' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out run_r.sh
  # sbatch --depend=afterok:102 --export=run='last',cluster='bionic',type='R',project='z.hivmi' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out run_r.sh
  
  run_on_cluster(
    
    first = {
      
      # install.packages(
      #   pkgs = 'dplyr', # data.table, survival, tidyr
      #   lib = '/home/akenny/R_lib',
      #   repos = 'http://cran.us.r-project.org',
      #   dependencies = TRUE
      # )
      
      # Declare functions
      source("generate_dataset.R")
      source("run_analysis.R")
      source("perform_imputation.R")
      source("transform_dataset.R")
      source("one_simulation.R")
      source("helpers.R")
      
      # Set up and configure simulation object
      sim <- new_sim()
      sim %<>% set_config(
        num_sim = 1, # !!!!!
        # stop_at_error = TRUE,
        # parallel = "outer",
        packages = c("dplyr", "survival", "data.table", "tidyr")
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
      sim %<>% add_script(one_simulation)
      sim %<>% add_method(convert_p_sero)
      sim %<>% add_method(construct_m_probs)
      sim %<>% add_method(run_analysis)
      sim %<>% add_method(perform_imputation)
      sim %<>% add_method(transform_dataset)
      sim %<>% add_method(p_death_year)
      
      # Set levels
      sim %<>% set_levels(
        method = c("ideal", "mi"),
        # method = c("ideal", "censor", "mi"),
        hr_hiv = c(0.5,1.0,1.4,1.8),
        hr_art = c(0.5,1.0,1.4,1.8),
        p_death_mult = 1,
        include_no_testers = TRUE,
        # include_no_testers = c(TRUE,FALSE), # FALSE by default for "censor" method
        u_mult_sero = 1, # 1.5
        u_mult_death = 1 # 1.3
      )
      
    },
    
    main = {
      sim %<>% run("one_simulation")
    },
    
    last = {
      # !!!!! Calculate bias, variance inflation, power
      print(sim$results)
    },
    
    cluster_config = list(
      sim_var = "sim",
      js = "slurm",
      dir = "/home/akenny/z.hivmi"
    )
    
  )
  
}



###########################.
##### Process results #####
###########################.

if (run_results) {
  
  library(ggplot2)
  
  # !!!!! TO DO
  
}



##############################.
##### Run dataset checks #####
##############################.

if (run_checks) {
  
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
  
  # # !!!!!
  # dataset_orig2
  # dataset_imp2
  
  # Create "imputed" dataset
  dataset_imp <- perform_imputation(dataset_orig, p_sero_year)
  
  # Create transformed datasets
  dataset_cp_orig <- transform_dataset(dataset_orig)
  dataset_cp_imp <- transform_dataset(dataset_imp)
  
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
      # case==1 ~ "1. No testing data",
      # case>=2 ~ "2/3/4. Some testing data"
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
      # case==1 ~ "1. No testing data",
      # case>=2 ~ "2/3/4. Some testing data"
    )
  ) %>% group_by(start_year, casc_status, case) %>% summarize(num=n())
  ggplot(d4_imp, aes(x=start_year, y=num, color=casc_status)) +
    geom_line() + facet_wrap(~case, ncol=2)
  
  # Check 5: examine cascade status by case
  # !!!!! Should these actually be similar ?????
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
  run_analysis(dataset_cp_orig, method="ideal", include_no_testers=TRUE)
  run_analysis(dataset_cp_imp, method="ideal", include_no_testers=TRUE)
  
  # Check 7: Examine case counts
  xtabs(~case, data=dataset_orig)

  # Check 8: Examine death counts
  xtabs(~died, data=dataset_orig)
  
  # Check 9: Compare death counts by casc_status
  xtabs(~casc_status, data=filter(dataset_cp_orig, died==1))
  xtabs(~casc_status, data=filter(dataset_cp_imp, died==1))
  
  # Check 10: Compare death counts by case+casc_status
  xtabs(~case+casc_status, data=filter(dataset_cp_orig, died==1))
  xtabs(~case+casc_status, data=filter(dataset_cp_imp, died==1))
  
  # !!!!! Checks 10+11 collectively illustrate the problem
  
  # Check 12: examine datasets manually for abnormalities
  # write.table(dataset_orig, file="dataset_orig.csv", sep=",", row.names=FALSE)
  # write.table(dataset_imp, file="dataset_imp.csv", sep=",", row.names=FALSE)
  # write.table(dataset_cp_orig, file="dataset_cp_orig.csv", sep=",", row.names=FALSE)
  # write.table(dataset_cp_imp, file="dataset_cp_imp.csv", sep=",", row.names=FALSE)
  
  # !!!!! Look for weird age effects, especially among infants who seroconverted
  
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
