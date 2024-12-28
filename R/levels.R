# # Set global constants
# C <- list(alpha_1=0.5)

# Set simulation levels
if (cfg$run_sims && Sys.getenv("sim_run") %in% c("first", "")) {
  
  # L <- list(n=300,max_time=70,par=list(a_x=log(0.005), a_y=log(0.003), a_v=log(0.7),a_z=log(0.004),g_x=c(log(1.3),log(1.2)), g_y=c(log(1.2),log(1.1)),g_v=c(log(1.2),log(1.1)), g_z=c(log(1.2),log(1.1)),beta_x=log(1.5), beta_z=log(0.6)))
  
  level_sets <- list()
  
  # Simulation 1: basic
  # Figures: ...
  par_10 <- par_20 <- par_40 <- list(
    a_x = -3,
    g_x = c(0.3,0.2),
    t_x = -0.1,
    a_s = -1.6,
    g_s = c(0.5,0.3),
    t_s = 0.1,
    beta_x = 0.4,
    a_y = -3.5,
    g_y = c(0.2,0.1),
    t_y = -0.1,
    a_v = -0.8,
    g_v = c(0.2,0.1)
  )
  par_20$a_v <- -1.6
  par_10$a_v <- -2.4
  level_sets[["level_set_1"]] <- list(
    n = 1000,
    max_time = 20,
    model_version = 7,
    par = list(
      # "40% testing" = par_40,
      # "20% testing" = par_20,
      "10% testing" = par_10
    )
  )
  
  level_set <- level_sets[[cfg$sim_level_set]]
  
  # if (cfg$sim_level_set=="asdf") { cfg$keep = c(1:3,7:9,16:18,22:24) }
  
}
