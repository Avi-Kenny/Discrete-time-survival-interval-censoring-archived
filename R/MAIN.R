# Title: "Modeling HIV seroconversion dates"
# Author: Avi Kenny
# Date: 2020-09-30



#################.
##### Setup #####
#################.

# Set working directory
# If running this locally in RStudio, skip this section and click "Session" >>
#     "Set working directory" >> "To Source File Location"
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
  setwd("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Research/Jim Hughes/Project - Stepped wedge/z.stepped.wedge/R")
} else {
  setwd("z.stepped.wedge/R")
}

# Set code blocks to run
{
  run_main <- TRUE
  run_results <- FALSE
}



###########################.
##### Simulation code #####
###########################.

if (run_main) {
  
  library(simba)
  library(magrittr)
  
  run_on_cluster(
    
    first = {
      
      # Set up and configure simulation object
      sim <- new_sim()
      sim %<>% set_config(
        num_sim = 10, # !!!!!
        parallel = "cluster",
        packages = c("dplyr", "mice", "survival")
      )
      # sim %<>% add_constant(alpha=log(0.1))
      
      # Add functions to simulation object
      source("generate_dataset.R")
      source("one_simulation.R")
      sim %<>% add_creator(generate_dataset)
      sim %<>% add_script(one_simulation)
      
      # Set levels
      sim %<>% set_levels(
        n_clusters = 24,
        n_time_points = 7
      )
      
    },
    
    main = { sim %<>% run("one_simulation") },
    
    last = {
      
      # Calculate bias, variance inflation, power
      
      sim %>% summary()
      
    },
    
    cluster_config = list(
      sim_var = "sim",
      js = "slurm",
      dir = "/home/akenny/HIV.multiple.imputation"
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
  
}
