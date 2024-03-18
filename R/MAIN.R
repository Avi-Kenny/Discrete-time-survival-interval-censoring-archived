# Main config
cfg <- list(
  run_sims = F,
  run_process = F,
  run_analysis = T,
  sim_level_set = "level_set_1",
  sim_run_or_update = "run",
  sim_num = 1000,
  sim_parallel = F,
  # sim_n_cores = 50,
  sim_n_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")),
  sim_stop_at_error = F,
  model_version = 7
)

# Secondary config
source("R/config.R", local=T)

# Load SimEngine + functions
{
  library(SimEngine)
  source("R/one_simulation.R", local=T)
  source("R/generate_data.R", local=T)
  source("R/likelihood_miss.R", local=T)
  source("R/helpers.R", local=T)
}

if (cfg$run_analysis) {
  
  set.seed(1)
  source("R/analysis.R", local=T)
  
} else if (cfg$run_sims) {
  
  # Set level sets
  source("R/levels.R", local=T)
  
  # Run simulation
  source("R/run.R", local=T)
  
} else if (cfg$run_process) {
  
  # Tables and figures
  source("R/process.R", local=T)
  
}
