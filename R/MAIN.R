# Title: "Modeling HIV seroconversion dates"
# Author: Avi Kenny



##################.
##### CONFIG #####
##################.
# install.packages(
#   pkgs = "stats",
#   lib = "/home/akenny/R_lib",
#   repos = "http://cran.us.r-project.org",
#   dependencies = TRUE
# )
# Set global config
cfg <- list(
  level_set_which = "level_set_1",
  # num_sim = 1, # 1000
  num_sim = 500, # 1000
  pkgs = c("dplyr", "survival", "data.table", "tidyr", "memoise", "Rsolnp",
           "numDeriv"),
  pkgs_nocluster = c("ggplot2"),
  parallel = F,
  n_cores = 500,
  stop_at_error = F
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
  source("one_simulation.R", local=T)
  source("generate_data.R", local=T)
  source("likelihood_miss.R", local=T)
  source("helpers.R", local=T)
}



##########################################################.
##### MAIN: Set level sets for different simulations #####
##########################################################.

if (Sys.getenv("sim_run") %in% c("first", "")) {
  
  # L <- list(n=300,max_time=70,params=list(a_x=log(0.005), a_y=log(0.003), a_v=log(0.7),a_z=log(0.004),g_x=c(log(1.3),log(1.2)), g_y=c(log(1.2),log(1.1)),g_v=c(log(1.2),log(1.1)), g_z=c(log(1.2),log(1.1)),beta_x=log(1.5), beta_z=log(0.6)))
  
  # Simulation 1: basic
  # n=500,t=100 rep runs in 3.2 hrs # !!!!! outdated
  level_set_1 <- list(
    # n = 50,
    n = 500,
    # n = c(500,1000,2000),
    max_time = 70,
    # max_time = 100,
    params = list(
      # "10pct testing" = list( a_v=log(0.1) ),
      "70pct testing" = list(
        a_x=log(0.005), a_y=log(0.003), a_v=log(0.7), a_z=log(0.01),
        g_x=c(log(1.3),log(1.2)), g_y=c(log(1.2),log(1.1)),
        g_v=c(log(1.2),log(1.1)), g_z=c(log(1.2),log(1.1)),
        beta_x=log(1.5), beta_z=log(0.7)
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
      n_cores = cfg$n_cores,
      stop_at_error = cfg$stop_at_error,
      seed = 123,
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
  
  v <- c("lik_M_a_x_est", "lik_M_g_x1_est", "lik_M_g_x2_est", "lik_M_a_y_est",
         "lik_M_g_y1_est", "lik_M_g_y2_est", "lik_M_beta_x_est",
         "lik_M_beta_z_est", "lik_M_t_x_est", "lik_M_t_y_est")
  true_vals <- log(c(0.005,1.3,1.2,0.003,1.2,1.1,1.5,0.7,1,1))
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
  summ <- sim %>% SimEngine::summarize(
    list(stat="mean", name="a_x__sd_est", x="lik_M_a_x_se"),
    list(stat="mean", name="g_x1__sd_est", x="lik_M_g_x1_se"),
    list(stat="mean", name="g_x2__sd_est", x="lik_M_g_x2_se"),
    list(stat="mean", name="a_y__sd_est", x="lik_M_a_y_se"),
    list(stat="mean", name="g_y1__sd_est", x="lik_M_g_y1_se"),
    list(stat="mean", name="g_y2__sd_est", x="lik_M_g_y2_se"),
    list(stat="mean", name="beta_x__sd_est", x="lik_M_beta_x_se"),
    list(stat="mean", name="beta_z__sd_est", x="lik_M_beta_z_se"),
    list(stat="mean", name="t_x__sd_est", x="lik_M_t_x_se"),
    list(stat="mean", name="t_y__sd_est", x="lik_M_t_y_se"),
    list(stat="sd", name="a_x__sd_actual", x="lik_M_a_x_est"),
    list(stat="sd", name="g_x1__sd_actual", x="lik_M_g_x1_est"),
    list(stat="sd", name="g_x2__sd_actual", x="lik_M_g_x2_est"),
    list(stat="sd", name="a_y__sd_actual", x="lik_M_a_y_est"),
    list(stat="sd", name="g_y1__sd_actual", x="lik_M_g_y1_est"),
    list(stat="sd", name="g_y2__sd_actual", x="lik_M_g_y2_est"),
    list(stat="sd", name="beta_x__sd_actual", x="lik_M_beta_x_est"),
    list(stat="sd", name="beta_z__sd_actual", x="lik_M_beta_z_est"),
    list(stat="sd", name="t_x__sd_actual", x="lik_M_t_x_est"),
    list(stat="sd", name="t_y__sd_actual", x="lik_M_t_y_est"),
    list(stat="coverage", name="a_x__cov", truth=true_vals[1],
         estimate="lik_M_a_x_est", se="lik_M_a_x_se", na.rm=T),
    list(stat="coverage", name="g_x1__cov", truth=true_vals[2],
         estimate="lik_M_g_x1_est", se="lik_M_g_x1_se", na.rm=T),
    list(stat="coverage", name="g_x2__cov", truth=true_vals[3],
         estimate="lik_M_g_x2_est", se="lik_M_g_x2_se", na.rm=T),
    list(stat="coverage", name="a_y__cov", truth=true_vals[4],
         estimate="lik_M_a_y_est", se="lik_M_a_y_se", na.rm=T),
    list(stat="coverage", name="g_y1__cov", truth=true_vals[5],
         estimate="lik_M_g_y1_est", se="lik_M_g_y1_se", na.rm=T),
    list(stat="coverage", name="g_y2__cov", truth=true_vals[6],
         estimate="lik_M_g_y2_est", se="lik_M_g_y2_se", na.rm=T),
    list(stat="coverage", name="beta_x__cov", truth=true_vals[7],
         estimate="lik_M_beta_x_est", se="lik_M_beta_x_se", na.rm=T),
    list(stat="coverage", name="beta_z__cov", truth=true_vals[8],
         estimate="lik_M_beta_z_est", se="lik_M_beta_z_se", na.rm=T),
    list(stat="coverage", name="t_x__cov", truth=true_vals[9],
         estimate="lik_M_t_x_est", se="lik_M_t_x_se", na.rm=T),
    list(stat="coverage", name="t_y__cov", truth=true_vals[10],
         estimate="lik_M_t_y_est", se="lik_M_t_y_se", na.rm=T)
  )
  
  l_id <- 1
  summ2 <- summ[summ$level_id==l_id]
  df_results <- data.frame(
    "param" = character(),
    "sd_est" = double(),
    "sd_actual" = double(),
    "coverage" = double()
  )
  params <- c("a_x", "g_x1", "g_x2", "a_y", "g_y1", "g_y2", "beta_x", "beta_z")
  for (p in params) {
    df_results[nrow(df_results)+1,] <- c(
      p,
      round(summ[[paste0(p,"__sd_actual")]], 3),
      round(summ[[paste0(p,"__sd_est")]], 3),
      round(summ[[paste0(p,"__cov")]], 3)
    )
  }
  print(df_results)
  
}
