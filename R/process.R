#################.
##### Setup #####
#################.

cfg2 <- list(
  process_sims = F,
  process_analysis = T,
  # m = 24,
  m = 30,
  w_start = 2010,
  w_end = 2022,
  # ests = readRDS("objs/ests_28_full_20240926.rds")
  # ests_M = readRDS("objs/ests_30_full_M_20240927.rds"),
  # ests_F = readRDS("objs/ests_30_full_F_20240927.rds")
  ests_M = readRDS("objs/ests_31_full_M_20241015.rds"),
  ests_F = readRDS("objs/ests_31_full_F_20241015.rds")
)

# Construct spline bases
b1 <- construct_basis("age (0-100), 4DF")
b2 <- construct_basis("age (13,20,30,60,90)")
b3 <- construct_basis("age (13,30,60,75,90)")
b4 <- construct_basis("year (00,05,10,15,20)", window_start=cfg2$w_start)
b5 <- construct_basis("year (10,13,17,20,23)", window_start=cfg2$w_start)
b6 <- construct_basis("age (13,28,44,60,75)")
b7 <- construct_basis("year (10,13,17,20,23) +i", window_start=cfg2$w_start)
b8 <- construct_basis("age (13,28,44,60,75) +i")
b9 <- construct_basis("age (13,20,30,40,60)")
b10 <- construct_basis("year (10,13,16,19,22)", window_start=cfg2$w_start,
                       window_end=cfg2$w_end)
b11 <- construct_basis("year (10,13,16,19,22) +i", window_start=cfg2$w_start,
                       window_end=cfg2$w_end)

# Get current date
cfg2$d <- format(Sys.time(), "%Y-%m-%d")



#####################################################.
##### VIZ (simulations): all params (one model) #####
#####################################################.

if (cfg2$process_sims) {
  
  sim <- readRDS("SimEngine.out/sim_20240719.rds")
  
  p_set <- "10% testing"
  
  # p_names should match those used by negloglik()
  p_names <- c("a_x", "g_x1", "g_x2", "t_x",
               "a_s", "g_s1", "g_s2", "t_s",
               "beta_x",
               "a_y", "g_y1", "g_y2", "t_y")
  p <- sim$levels$params[[p_set]]
  true_vals <- c(p$a_x, p$g_x, p$t_x,
                 p$a_s, p$g_s, p$t_s,
                 p$beta_x,
                 p$a_y, p$g_y, p$t_y)
  names(true_vals) <- p_names
  
  # prm_sim <- sim$levels$params[[1]]
  # true_vals2 <- c(prm_sim$a_x, prm_sim$g_x[1], prm_sim$g_x[2], prm_sim$a_y,
  #                prm_sim$g_y[1], prm_sim$g_y[2], prm_sim$beta_x, prm_sim$beta_z,
  #                # prm_sim$t_x, prm_sim$t_y,
  #                prm_sim$a_s, prm_sim$t_s,
  #                prm_sim$g_s[1], prm_sim$g_s[2])
  
  v <- paste0("lik_M_",p_names,"_est")
  r <- dplyr::filter(sim$results, params==p_set)
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
  
  plot <- ggplot(df_plot, aes(x=x, y=y, color=which)) +
    geom_jitter(width=0, height=1, alpha=0.3, size=3) +
    geom_vline(aes(xintercept=val), df_true, alpha=0.5, linetype="dashed") +
    facet_wrap(~which, ncol=1, strip.position="left") + # scales="free_x"
    labs(y=NULL, title=p_set) +
    ylim(-2,2) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text.y.left = element_text(angle=0),
          legend.position="none")
  
  ggsave(
    filename=paste0("../Figures + Tables/", cfg2$d, " sim_est_scatterplot.pdf"),
    plot=plot, device="pdf", width=8, height=5
  )
  
  # Summary stats
  summ_mean <- summ_bias <- summ_mean2 <- summ_sd <- summ_cov <- list()
  for (i in c(1:length(p_names))) {
    p <- p_names[i]
    summ_mean[[i]] <- list(stat="mean",
                           name=paste0(p,"__est"),
                           x=paste0("lik_M_",p,"_est"), na.rm=T)
    summ_bias[[i]] <- list(stat="bias",
                           name=paste0(p,"__bias"),
                           estimate=paste0("lik_M_",p,"_est"),
                           truth=true_vals[[p]])
    summ_mean2[[i]] <- list(stat="mean",
                           name=paste0(p,"__sd_est"),
                           x=paste0("lik_M_",p,"_se"), na.rm=T)
    summ_sd[[i]] <- list(stat="sd", name=paste0(p,"__sd_actual"),
                         x=paste0("lik_M_",p,"_est"), na.rm=T)
    summ_cov[[i]] <- list(stat="coverage", name=paste0(p,"__cov"),
                          truth=true_vals[i],
                          estimate=paste0("lik_M_",p,"_est"),
                          se=paste0("lik_M_",p,"_se"), na.rm=T)
  }
  summ <- do.call(
    SimEngine::summarize,
    c(list(sim), summ_mean, summ_bias, summ_mean2, summ_sd, summ_cov)
  )
  l_id <- 1
  summ2 <- summ[summ$level_id==l_id]
  df_results <- data.frame(
    "param" = character(),
    "truth" = double(),
    "est" = double(),
    "bias" = double(),
    "sd_est" = double(),
    "sd_actual" = double(),
    "coverage" = double()
  )
  for (p in p_names) {
    df_results[nrow(df_results)+1,] <- c(
      p,
      true_vals[[p]],
      round(summ[[paste0(p,"__est")]], 3),
      round(summ[[paste0(p,"__bias")]], 3),
      round(summ[[paste0(p,"__sd_est")]], 3),
      round(summ[[paste0(p,"__sd_actual")]], 3),
      round(summ[[paste0(p,"__cov")]], 3)
    )
  }
  utils::write.table(
    x = df_results,
    file = paste0("../Figures + Tables/", cfg2$d, " sims_sd_and_coverage.csv"),
    sep = ",",
    row.names = FALSE
  )
  
}



########################################################.
##### Functions for plotting modeled probabilities #####
########################################################.

#' Return modeled probability
#' @param type One of c("sero", "init", "mort (HIV-)", "mort (HIV+)",
#'     "mort (HIV+ART-)", "mort (HIV+ART+)")
#' @param j Calendar time (in years, starting at 1, where 1 represents the
#'     first year in the dataset)
#' @param w_1 Age (in years)
#' @param w_2 Sex (0=female, 1=male)
#' @param m An integer representing the model version number
#' @param which One of c("est", "ci_lo", "ci_up")
#' @return Numeric probability
prob <- function(type, m, j, w_1, w_2, which="est") {
  
  j <- j/10
  w_1 <- w_1/100
  
  if (type=="sero") {
    
    if (m %in% c(29:31)) {
      
      A <- t(matrix(c(
        1, b10(j,1), b10(j,2), b10(j,3), b10(j,4), b9(w_1,1), b9(w_1,2),
        b9(w_1,3), b9(w_1,4)
      )))
      p2 <- c("a_x", "t_x1", "t_x2", "t_x3", "t_x4", "g_x1", "g_x2", "g_x3",
              "g_x4")
      
    }
    
  } else if (type=="init") {
    
    if (m %in% c(29:31)) {
      A <- t(matrix(c(
        1, b9(w_1,1), b9(w_1,2), b9(w_1,3), b9(w_1,4)
      )))
      p2 <- c("a_s", "g_s1", "g_s2", "g_s3", "g_s4")
    }
    
  } else {
    
    if (type=="mort (HIV-)") { x <- 0 }
    if (type=="mort (HIV+)") { x <- 1 }

    if (m==29) {
      
      A <- t(matrix(c(
        x, x*j, x*max(w_1-0.3,0), x*j*max(w_1-0.3,0), x*max(w_1-0.45,0),
        x*j*max(w_1-0.45,0), 1, b9(w_1,1), b9(w_1,2), b9(w_1,3), b9(w_1,4),
        b10(j,1), b10(j,2), b10(j,3), b10(j,4)
      )))
      p2 <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5", "beta_x6",
              "a_y", "g_y1", "g_y2", "g_y3", "g_y4", "t_y1", "t_y2", "t_y3",
              "t_y4")
      
    } else if (m==30) {
      
      A <- t(matrix(c(
        x, x*j, x*w_1, x*j*w_1, 1, b9(w_1,1), b9(w_1,2), b9(w_1,3), b9(w_1,4),
        b10(j,1), b10(j,2), b10(j,3), b10(j,4)
      )))
      p2 <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "a_y", "g_y1", "g_y2",
              "g_y3", "g_y4", "t_y1", "t_y2", "t_y3", "t_y4")
      
    } else if (m==31) {
      
      A <- t(matrix(c(
        x, x*j, x*max(w_1-0.3,0), x*j*max(w_1-0.3,0), x*max(w_1-0.5,0),
        x*j*max(w_1-0.5,0), 1, b9(w_1,1), b9(w_1,2), b9(w_1,3), b9(w_1,4),
        b10(j,1), b10(j,2), b10(j,3), b10(j,4)
      )))
      p2 <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5", "beta_x6",
              "a_y", "g_y1", "g_y2", "g_y3", "g_y4", "t_y1", "t_y2", "t_y3",
              "t_y4")
      
    }
    
  }
  
  if (w_2==1) { ests <- cfg2$ests_M } else { ests <- cfg2$ests_F }
  indices <- as.numeric(sapply(p2, function(p) {
    which(names(ests$opt$par)==p)
  }))
  beta <- matrix(ests$opt$par[indices])
  Sigma <- ests$hessian_inv[indices,indices]
  if (which=="est") {
    fac <- 0
  } else if (which=="ci_lo") {
    fac <- -1.96
  } else if (which=="ci_up") {
    fac <- 1.96
  }
  est <- c(A %*% beta)
  se <- c(sqrt(A %*% Sigma %*% t(A)))
  # se <- 0.1 # !!!!! For debugging
  prob <- icll(est+fac*se)
  
  return(prob)
  
}



#' Return plot of modeled probabilities
#' @param x_axis One of c("Year", "Age"); the variable to go on the X-axis
#' @param type One of c("sero", "init", "mort (HIV-)", "mort (HIV+)",
#'     "mort (HIV+ART-)", "mort (HIV+ART+)")
#'     
#' @param m An integer representing the model version number
#' @param w_start An integer representing the window start calendar year
#' @param w_end An integer representing the window end calendar year
#' @param y_max Maximum Y value for the plot
#' @return ggplot2 object
plot_mod <- function(x_axis, type, m, w_start, w_end, y_max) {
  
  if (w_start==2000) {
    j_vals <- c(1,11,21) # This was formerly incorrectly set to c(0,10,20)
    j_labs <- c("2000","2010","2020")
  } else if (w_start==2010) {
    j_vals <- c(1,6,11)
    j_labs <- c("2010","2015","2020")
  }
  
  if (x_axis=="Age") {
    grid <- seq(13,60,0.1)
    color <- "Year"
    plot_data <- data.frame(
      x = rep(grid,6),
      Probability = c(
        sapply(grid, function(w_1) { prob(type, m, j=j_vals[1], w_1, w_2=0) }),
        sapply(grid, function(w_1) { prob(type, m, j=j_vals[1], w_1, w_2=1) }),
        sapply(grid, function(w_1) { prob(type, m, j=j_vals[2], w_1, w_2=0) }),
        sapply(grid, function(w_1) { prob(type, m, j=j_vals[2], w_1, w_2=1) }),
        sapply(grid, function(w_1) { prob(type, m, j=j_vals[3], w_1, w_2=0) }),
        sapply(grid, function(w_1) { prob(type, m, j=j_vals[3], w_1, w_2=1) })
      ),
      Sex = rep(rep(c("Female", "Male"),3), each=length(grid)),
      color = rep(j_labs, each=2*length(grid))
    )
  } else if (x_axis=="Year") {
    grid <- seq(w_start,w_end,0.01) %>% (function(x) { x-(w_start-1) })
    color <- "Age"
    plot_data <- data.frame(
      x = rep(grid,6) + (w_start-1),
      Probability = c(
        sapply(grid, function(j) { prob(type, m, j, w_1=20, w_2=0) }),
        sapply(grid, function(j) { prob(type, m, j, w_1=20, w_2=1) }),
        sapply(grid, function(j) { prob(type, m, j, w_1=35, w_2=0) }),
        sapply(grid, function(j) { prob(type, m, j, w_1=35, w_2=1) }),
        sapply(grid, function(j) { prob(type, m, j, w_1=50, w_2=0) }),
        sapply(grid, function(j) { prob(type, m, j, w_1=50, w_2=1) })
      ),
      Sex = rep(rep(c("Female", "Male"),3), each=length(grid)),
      color = rep(c("20","35","50"), each=2*length(grid))
    )
  }
  plot_data %<>% dplyr::mutate(grp=factor(paste0(Sex,color)))
  
  # Seroconversion prob as function of age
  if (type=="sero") {
    title = "Seroconversion prob (per year); model"
  } else if (type=="init") {
    title = "Prob an individial's initial status is HIV+; model"
  } else if (type=="mort (HIV-)") {
    title = "Mortality prob among HIV- (per year); model"
  } else if (type=="mort (HIV+)") {
    title = "Mortality prob among HIV+ (per year); model"
  } else if (type=="mort (HIV+ART-)") {
    title = "Mortality prob among HIV+ART- (per year); model"
  } else if (type=="mort (HIV+ART+)") {
    title = "Mortality prob among HIV+ART+ (per year); model"
  }
  title <- paste0(title, " v", m)
  plot <- ggplot(
    plot_data,
    aes(x=x, y=Probability, color=color, group=grp,
        linetype=Sex)
  ) +
    geom_line() +
    coord_cartesian(ylim=c(0,y_max)) +
    labs(title=title, x=x_axis, color=color) +
    theme(plot.background = element_rect(color="black"))
  
  return(plot)
  
}



#' Return plot of modeled mortality probabilities (HIV- vs. HIV+) with CIs
#' @param x_axis One of c("Year", "Age"); the variable to go on the X-axis
#' @param m An integer representing the model version number
#' @param w_start An integer representing the window start calendar year
#' @param w_start An integer representing the window end calendar year
#' @param y_max Maximum Y value for the plot
#' @param title Boolean; if F, title is suppressed
#' @return ggplot2 object
plot_mort3 <- function(x_axis, m, w_start, w_end, y_max=NA, title=T) {
  
  if (x_axis=="Age") {
    
    grid <- seq(13,60,0.1)
    breaks <- seq(20,60, length.out=5)
    color <- "Year"
    if (w_start==2000) {
      outer <- c(1,11,21)
    } else if (w_start==2010) {
      outer <- c(1,6,11)
    }
    prob2 <- function(type, which, outer) {
      1000 * sapply(grid, function(inner) {
        prob(type=type, m=m, j=outer, w_1=inner, w_2=sex, which=which)
      })
    }
    
  } else if (x_axis=="Year") {
    
    grid <- seq(w_start,w_end,0.01) %>% (function(x) { x-(w_start-1) })
    breaks <- seq(w_start,w_end, length.out=5)
    color <- "Age"
    outer <- c(20,35,50)
    prob2 <- function(type, which, outer) {
      1000 * sapply(grid, function(inner) {
        prob(type=type, m=m, j=inner, w_1=outer, w_2=sex, which=which)
      })
    }
    
  }
  
  init <- T
  
  for (o in outer) {
    for (sex in c(0,1)) {
      
      plot_data <- data.frame(
        x = rep(grid,2),
        Rate = c(prob2(type="mort (HIV-)", which="est", outer=o),
                 prob2(type="mort (HIV+)", which="est", outer=o)),
        ci_lo = c(prob2(type="mort (HIV-)", which="ci_lo", outer=o),
                  prob2(type="mort (HIV+)", which="ci_lo", outer=o)),
        ci_up = c(prob2(type="mort (HIV-)", which="ci_up", outer=o),
                  prob2(type="mort (HIV+)", which="ci_up", outer=o)),
        color = rep(c("HIV-","HIV+"), each=length(grid)),
        outer = o,
        sex = ifelse(sex, "Male", "Female")
      )
      
      if (x_axis=="Year") {
        plot_data$x <- plot_data$x + (w_start-1)
        plot_data$outer <- paste0("Age: ", plot_data$outer)
      } else if (x_axis=="Age") {
        plot_data$outer <- plot_data$outer + (w_start-1)
      }
      
      if (init) {
        plot_data2 <- plot_data
        init <- F
      } else {
        plot_data2 <- rbind(plot_data2, plot_data)
      }
      
    }
  }
  
  if (title) {
    title <- "Conditional mortality rates by HIV status, over calendar time"
  } else {
    title <- NULL
  }
  
  plot <- ggplot(
    plot_data2,
    aes(x=x, y=Rate, color=color)
  ) +
    geom_line() +
    geom_ribbon(
      aes(ymin=ci_lo, ymax=ci_up, fill=color),
      alpha = 0.2,
      linetype = "blank"
    ) +
    facet_grid(rows=dplyr::vars(sex), cols=dplyr::vars(outer)) +
    coord_cartesian(ylim=c(0,y_max)) +
    labs(
      title = title,
      x = x_axis,
      color = "HIV Status",
      fill = "HIV Status",
      y = "Deaths per 1,000 person-years"
    ) +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=breaks) +
    scale_color_manual(values=c("forestgreen", "#56B4E9")) +
    scale_fill_manual(values=c("forestgreen", "#56B4E9"))
  
  return(plot)
  
}



#' Return plot of modeled mortality probabilities (HIV- vs. HIV+) with CIs
#' @param m An integer representing the model version number
#' @param w_start An integer representing the window start calendar year
#' @param y_max Maximum Y value for the plot
#' @param title Boolean; if F, title is suppressed
#' @return ggplot2 object
plot_sero3 <- function(m, w_start, y_max=NA, title=T) {
  
  grid <- seq(13,60,0.1)
  breaks <- seq(20,60, length.out=5)
  color <- "Year"
  if (w_start==2000) {
    outer <- c(1,11,21)
  } else if (w_start==2010) {
    outer <- c(1,6,11)
  }
  prob2 <- function(which, outer) {
    sapply(grid, function(inner) {
      prob(type="sero", m=m, j=outer, w_1=inner, w_2=sex, which=which)
    })
  }
  
  init <- T
  
  for (o in outer) {
    for (sex in c(0,1)) {
      
      plot_data <- data.frame(
        x = rep(grid,2),
        Rate = prob2(which="est", outer=o),
        ci_lo = prob2(which="ci_lo", outer=o),
        ci_up = prob2(which="ci_up", outer=o),
        outer = o,
        sex = ifelse(sex, "Male", "Female")
      )
      
      plot_data$outer <- plot_data$outer + (w_start-1)
      
      if (init) {
        plot_data2 <- plot_data
        init <- F
      } else {
        plot_data2 <- rbind(plot_data2, plot_data)
      }
      
    }
  }
  
  if (title) {
    title <- "Probability of seroconversion, by age"
  } else {
    title <- NULL
  }
  
  plot <- ggplot(
    plot_data2,
    aes(x=x, y=Rate)
  ) +
    geom_line(color="forestgreen") +
    geom_ribbon(
      aes(ymin=ci_lo, ymax=ci_up),
      alpha = 0.2,
      fill = "forestgreen"
    ) +
    facet_grid(rows=dplyr::vars(sex), cols=dplyr::vars(outer)) +
    coord_cartesian(ylim=c(0,y_max)) +
    labs(
      title = title,
      x = "Age",
      y = "Probability of seroconversion (in one year)"
    ) +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=breaks)

  return(plot)
  
}



###############################################.
##### VIZ: Plotting modeled probabilities #####
###############################################.

if (cfg2$process_analysis) {
  
  y_max <- c(0.06, 0.8, 0.2, 0.05)
  
  # Seroconversion prob as a function of age
  p01 <- plot_mod(x_axis="Age", type="sero", m=cfg2$m, w_start=cfg2$w_start, w_end=cfg2$w_end, y_max=y_max[1])
  
  # Seroconversion prob as a function of calendar time
  p02 <- plot_mod(x_axis="Year", type="sero", m=cfg2$m, w_start=cfg2$w_start, w_end=cfg2$w_end, y_max=y_max[1])
  
  # HIV+ initial status as a function of age
  p03 <- plot_mod(x_axis="Age", type="init", m=cfg2$m, w_start=cfg2$w_start, w_end=cfg2$w_end, y_max=y_max[2])
  
  # HIV+ initial status as a function of calendar time
  p04 <- plot_mod(x_axis="Year", type="init", m=cfg2$m, w_start=cfg2$w_start, w_end=cfg2$w_end, y_max=y_max[2])
  
  # Mortality prob as a function of age
  p05 <- plot_mod(x_axis="Age", type="mort (HIV-)", m=cfg2$m, w_start=cfg2$w_start, w_end=cfg2$w_end, y_max=y_max[3])
  p06 <- plot_mod(x_axis="Age", type="mort (HIV+)", m=cfg2$m, w_start=cfg2$w_start, w_end=cfg2$w_end, y_max=y_max[3])
  
  # Mortality prob as a function of calendar time
  p08 <- plot_mod(x_axis="Year", type="mort (HIV-)", m=cfg2$m, w_start=cfg2$w_start, w_end=cfg2$w_end, y_max=y_max[4])
  p09 <- plot_mod(x_axis="Year", type="mort (HIV+)", m=cfg2$m, w_start=cfg2$w_start, w_end=cfg2$w_end, y_max=y_max[4])
  
  # Old ART plots
  if (F) {
    # Mortality prob as a function of age
    p05 <- plot_mod(x_axis="Age", type="mort (HIV-)", m=cfg2$m, w_start=cfg2$w_start, w_end=cfg2$w_end, y_max=y_max[3])
    p06 <- plot_mod(x_axis="Age", type="mort (HIV+ART-)", m=cfg2$m, w_start=cfg2$w_start, w_end=cfg2$w_end, y_max=y_max[3])
    p07 <- plot_mod(x_axis="Age", type="mort (HIV+ART+)", m=cfg2$m, w_start=cfg2$w_start, w_end=cfg2$w_end, y_max=y_max[3])
    
    # Mortality prob as a function of calendar time
    p08 <- plot_mod(x_axis="Year", type="mort (HIV-)", m=cfg2$m, w_start=cfg2$w_start, w_end=cfg2$w_end, y_max=y_max[4])
    p09 <- plot_mod(x_axis="Year", type="mort (HIV+ART-)", m=cfg2$m, w_start=cfg2$w_start, w_end=cfg2$w_end, y_max=y_max[4])
    p10 <- plot_mod(x_axis="Year", type="mort (HIV+ART+)", m=cfg2$m, w_start=cfg2$w_start, w_end=cfg2$w_end, y_max=y_max[4])
    
    # Combined plot
    plot_04 <- ggpubr::ggarrange(p05, p06, p07, p08, p09, p10, ncol=3, nrow=2)
  }
  
  # Plots without CIs
  plot_01 <- ggpubr::ggarrange(p01, p02)
  plot_02 <- ggpubr::ggarrange(p03, p04)
  plot_03 <- ggpubr::ggarrange(p05, p08, p06, p09)

  ggsave(
    filename = paste0("../Figures + Tables/", cfg2$d, " p1 (seroconversion) - ",
                      "model ", cfg2$m, ".pdf"),
    plot = plot_01, device="pdf", width=10, height=5
  )
  ggsave(
    filename = paste0("../Figures + Tables/", cfg2$d, " p2 (HIV initial status",
                      ") - model ", cfg2$m, ".pdf"),
    plot = plot_02, device="pdf", width=10, height=5
  )
  ggsave(
    filename = paste0("../Figures + Tables/", cfg2$d, " p3 (mortality) - model ",
                      cfg2$m, ".pdf"),
    plot = plot_03, device="pdf", width=10, height=10
  )
  
  # Mortality by calendar time, with CIs
  plot_05 <- plot_mort3(x_axis="Year", m=cfg2$m, w_start=cfg2$w_start,
                        w_end=cfg2$w_end, y_max=65, title=F)
  ggsave(
    filename = paste0("../Figures + Tables/", cfg2$d, " p5 (mort CIs, by year)",
                      " - model ", cfg2$m, ".pdf"),
    plot = plot_05, device="pdf", width=8, height=5
  )
  
  # Mortality by age, with CIs
  plot_06 <- plot_mort3(x_axis="Age", m=cfg2$m, w_start=cfg2$w_start,
                        w_end=cfg2$w_end, y_max=120, title=F)
  ggsave(
    filename = paste0("../Figures + Tables/", cfg2$d, " p6 (mort CIs, by age) ",
                      "- model ", cfg2$m, ".pdf"),
    plot = plot_06, device="pdf", width=8, height=5
  )
  
  # Seroconversion by age, with CIs
  plot_07 <- plot_sero3(m=cfg2$m, w_start=cfg2$w_start, y_max=0.05, title=F)
  ggsave(
    filename = paste0("../Figures + Tables/", cfg2$d, " p7 (sero CIs, by age) ",
                      "- model ", cfg2$m, ".pdf"),
    plot = plot_07, device="pdf", width=8, height=5
  )
  
}



#######################################################.
##### VIZ: Marginalized modeled probability curve #####
#######################################################.

# !!!!! Note: deprioritizing for now

# if (cfg2$process_analysis) {
#   
#   if (cfg2$m==29) {
#     dat_M <- readRDS("../Data/dat_29_full_M_20240926.rds")
#     dat_F <- readRDS("../Data/dat_29_full_F_20240926.rds")
#   } else if (cfg2$m==30) {
#     dat_M <- readRDS("../Data/dat_30_full_M_20240927.rds")
#     dat_F <- readRDS("../Data/dat_30_full_F_20240927.rds")
#   } else {
#     stop("")
#   }
#   
#   w_2_mean <- nrow(dat_M)/(nrow(dat_M)+nrow(dat_F))
#   # !!!!! COntinue here
#   
#   b9_1_mean <- mean(sapply(dat$w_1, function(w_1) { b9(w_1,1) }))
#   b9_2_mean <- mean(sapply(dat$w_1, function(w_1) { b9(w_1,2) }))
#   b9_3_mean <- mean(sapply(dat$w_1, function(w_1) { b9(w_1,3) }))
#   b9_4_mean <- mean(sapply(dat$w_1, function(w_1) { b9(w_1,4) }))
#   
#   prob_m <- function(type, m, j, which="est") {
#     
#     j <- j/10
#     
#     if (type=="mort (HIV-)") { x <- 0; z <- NA; }
#     if (type=="mort (HIV+)") { x <- 1; z <- NA; }
#     
#     if (m==25) {
#       A <- function(j) { t(matrix(c(
#         x*b7(j,1), x*b7(j,2), x*b7(j,3), x*b7(j,4), x*b7(j,5), 1, w_1_mean,
#         b9_1_mean, b9_2_mean, b9_3_mean, b9_4_mean, b5(j,1), b5(j,2), b5(j,3),
#         b5(j,4)
#       ))) }
#       p2 <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5", "a_y",
#               "g_y1", "g_y2", "g_y3", "g_y4", "g_y5", "t_y1", "t_y2", "t_y3",
#               "t_y4")
#     } else if (m==26) {
#       A <- function(j) { t(matrix(c(
#         x*b11(j,1), x*b11(j,2), x*b11(j,3), x*b11(j,4), x*b11(j,5), 1, w_1_mean,
#         b9_1_mean, b9_2_mean, b9_3_mean, b9_4_mean, b10(j,1), b10(j,2),
#         b10(j,3), b10(j,4)
#       ))) }
#       p2 <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5", "a_y",
#               "g_y1", "g_y2", "g_y3", "g_y4", "g_y5", "t_y1", "t_y2", "t_y3",
#               "t_y4")
#     }
#     
#     indices <- as.numeric(sapply(p2, function(p) {
#       which(names(cfg2$ests$opt$par)==p)
#     }))
#     beta <- matrix(cfg2$ests$opt$par[indices])
#     Sigma <- cfg2$ests$hessian_inv[indices,indices]
#     if (which=="est") {
#       fac <- 0
#     } else if (which=="ci_lo") {
#       fac <- -1.96
#     } else if (which=="ci_up") {
#       fac <- 1.96
#     }
#     est <- c(A(j) %*% beta)
#     se <- c(sqrt(A(j) %*% Sigma %*% t(A(j))))
#     prob <- icll(est+fac*se)
#     
#     return(prob)
#     
#   }
#   
#   grid <- seq(cfg2$w_start,cfg2$w_end,0.01) %>% (function(x) { x-(cfg2$w_start-1) })
#   breaks_x <- seq(cfg2$w_start,cfg2$w_end, length.out=5)
#   breaks_y <- seq(0,20,2)
#   prob2_m <- function(type, which, outer) {
#     1000 * sapply(grid, function(j) {
#       prob_m(type=type, m=cfg2$m, j=j, which=which)
#     })
#   }
#   
#   plot_data <- data.frame(
#     x = rep(grid,2) + (cfg2$w_start-1),
#     Rate = c(prob2_m(type="mort (HIV-)", which="est", outer=o),
#              prob2_m(type="mort (HIV+)", which="est", outer=o)),
#     ci_lo = c(prob2_m(type="mort (HIV-)", which="ci_lo", outer=o),
#               prob2_m(type="mort (HIV+)", which="ci_lo", outer=o)),
#     ci_up = c(prob2_m(type="mort (HIV-)", which="ci_up", outer=o),
#               prob2_m(type="mort (HIV+)", which="ci_up", outer=o)),
#     color = rep(c("HIV-","HIV+"), each=length(grid))
#   )
#   
#   # y_max <- 25
#   plot_08 <- ggplot(
#     plot_data,
#     aes(x=x, y=Rate, color=color)
#   ) +
#     geom_line() +
#     geom_ribbon(
#       aes(ymin=ci_lo, ymax=ci_up, fill=color),
#       alpha = 0.2,
#       linetype = "blank"
#     ) +
#     # coord_cartesian(ylim=c(0,y_max)) +
#     labs(
#       title = "Marginalized mortality rates by HIV status, over calendar time",
#       x = "Year",
#       color = "HIV Status",
#       fill = "HIV Status",
#       y = "Deaths per 1,000 person-years"
#     ) +
#     theme(legend.position = "bottom") +
#     scale_x_continuous(breaks=breaks_x) +
#     scale_y_continuous(breaks=breaks_y) +
#     scale_color_manual(values=c("forestgreen", "#56B4E9")) +
#     scale_fill_manual(values=c("forestgreen", "#56B4E9"))
#   
#   ggsave(
#     filename = paste0("../Figures + Tables/", cfg2$d, " p8 (mort CIs, by year,",
#                       " marginal) - model ", cfg2$m, ".pdf"),
#     plot = plot_08, device="pdf", width=8, height=5
#   )
#   
# }



#############################################.
##### VIZ: HRs as functions of age/time #####
#############################################.

if (cfg2$process_analysis) {
  
  # Functions to return spline basis function as a matrix
  A_b5 <- function(j) {
    t(matrix(c(b5(j,1),b5(j,2),b5(j,3),b5(j,4))))
  }
  A_b6 <- function(w_2) {
    t(matrix(c(b6(w_2,1), b6(w_2,2), b6(w_2,3), b6(w_2,4))))
  }
  A_b7 <- function(j) {
    t(matrix(c(b7(j,1),b7(j,2),b7(j,3),b7(j,4),b7(j,5))))
  }
  A_b8 <- function(w_2) {
    t(matrix(c(b8(w_2,1), b8(w_2,2), b8(w_2,3), b8(w_2,4), b8(w_2,5))))
  }
  A_b9 <- function(w_2) {
    t(matrix(c(b9(w_2,1), b9(w_2,2), b9(w_2,3), b9(w_2,4))))
  }
  A_b10 <- function(j) {
    t(matrix(c(b10(j,1),b10(j,2),b10(j,3),b10(j,4))))
  }
  A_b11 <- function(j) {
    t(matrix(c(b11(j,1),b11(j,2),b11(j,3),b11(j,4),b11(j,5))))
  }
  
  # Extract estimates and SEs
  plot_names <- c(
    "HR_mort_hiv_cal",
    "HR_mort_age",
    "HR_mort_cal",
    "HR_sero_cal",
    # "HR_sero_male_age",
    # "HR_sero_female_age",
    "HR_sero_age",
    "HR_init_age"
  )
  for (plot_name in plot_names) {
    
    # Set graph-specific variables
    if (plot_name=="HR_mort_hiv_cal") {
      title <- "HR of mortality, HIV+ vs. HIV- individuals"
      x_axis <- "cal time"
      # if (cfg2$m %in% c(23:25)) {
      #   params <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5")
      #   A <- A_b7
      # } else if (cfg2$m==26) {
      #   params <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5")
      #   A <- A_b11
      # }
      if (cfg2$m==30) {
        params <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4")
        A <- function(j, w_1) { t(matrix(c(1, j, w_1, j*w_1))) }
      } else if (cfg2$m==31) {
        params <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5",
                    "beta_x6")
        A <- function(j, w_1) {
          t(matrix(c(1, j, max(w_1-0.3,0), j*max(w_1-0.3,0), max(w_1-0.5,0),
                     j*max(w_1-0.5,0))))
        }
      }
    } else if (plot_name=="HR_mort_age") {
      title <- "HR of mortality (age)"
      x_axis <- "age"
      if (cfg2$m==24) {
        params <- c("g_y2", "g_y3", "g_y4", "g_y5")
        A <- A_b6
      } else if (cfg2$m %in% c(25:26)) {
        params <- c("g_y2", "g_y3", "g_y4", "g_y5")
        A <- A_b9
      } else if (cfg2$m %in% c(30:31)) {
        params <- c("g_y1", "g_y2", "g_y3", "g_y4")
        A <- A_b9
      }
    } else if (plot_name=="HR_mort_cal") {
      title <- "HR of mortality (calendar time)"
      x_axis <- "cal time"
      if (cfg2$m %in% c(24:25)) {
        params <- c("t_y1", "t_y2", "t_y3", "t_y4")
        A <- A_b5
      } else if (cfg2$m %in% c(26,30,31)) {
        params <- c("t_y1", "t_y2", "t_y3", "t_y4")
        A <- A_b10
      }
    } else if (plot_name=="HR_sero_cal") {
      title <- "HR of seroconversion (calendar time)"
      if (cfg2$m %in% c(23:25)) {
        params <- c("t_x1", "t_x2", "t_x3", "t_x4")
        A <- A_b5
      } else if (cfg2$m %in% c(26,30,31)) {
        params <- c("t_x1", "t_x2", "t_x3", "t_x4")
        A <- A_b10
      }
      x_axis <- "cal time"
    } else if (plot_name=="HR_sero_age") {
      title <- "HR of seroconversion (age)"
      x_axis <- "age"
      if (cfg2$m %in% c(30:31)) {
        params <- c("g_x1", "g_x2", "g_x3", "g_x4")
        A <- A_b9
      }
    # } else if (plot_name=="HR_sero_male_age") {
    #   title <- "HR of seroconversion among males (age)"
    #   x_axis <- "age"
    #   if (cfg2$m==23) {
    #     params <- c("g_x1", "g_x2", "g_x3", "g_x4", "g_x5")
    #     A <- A_b8
    #   } else if (cfg2$m==24) {
    #     params <- c("g_x1", "g_x2", "g_x3", "g_x4")
    #     A <- A_b6
    #   } else if (cfg2$m %in% c(25:26)) {
    #     params <- c("g_x1", "g_x2", "g_x3", "g_x4")
    #     A <- A_b9
    #   } else if (cfg2$m==30) {
    #     
    #   }
    # } else if (plot_name=="HR_sero_female_age") {
    #   title <- "HR of seroconversion among females (age)"
    #   x_axis <- "age"
    #   if (cfg2$m==23) {
    #     params <- c("g_x6", "g_x7", "g_x8", "g_x9", "g_x10")
    #     A <- A_b8
    #   } else if (cfg2$m==24) {
    #     params <- c("g_x5", "g_x6", "g_x7", "g_x8")
    #     A <- A_b6
    #   } else if (cfg2$m %in% c(25:26)) {
    #     params <- c("g_x5", "g_x6", "g_x7", "g_x8")
    #     A <- A_b9
    #   } else if (cfg2$m==30) {
    #     
    #   }
    } else if (plot_name=="HR_init_age") {
      title <- "HR of HIV+ initial status (age)"
      x_axis <- "age"
      if (cfg2$m==24) {
        params <- c("g_s2", "g_s3", "g_s4", "g_s5")
        A <- A_b6
      } else if (cfg2$m %in% c(25:26)) {
        params <- c("g_s2", "g_s3", "g_s4", "g_s5")
        A <- A_b9
      } else if (cfg2$m %in% c(30:31)) {
        params <- c("g_s1", "g_s2", "g_s3", "g_s4")
        A <- A_b9
      }
    }
    
    if (cfg2$m %in% c(24:26)) {
      indices <- which(names(cfg2$ests$opt$par) %in% params)
      beta <- matrix(cfg2$ests$opt$par[indices])
      Sigma <- cfg2$ests$hessian_inv[indices,indices]
      indices_M <- indices_F <- indices
      beta_M <- beta_F <- beta
      Sigma_M <- Sigma_F <- Sigma
    } else if (cfg2$m %in% c(30:31)) {
      indices_M <- which(names(cfg2$ests_M$opt$par) %in% params)
      beta_M <- matrix(cfg2$ests_M$opt$par[indices_M])
      Sigma_M <- cfg2$ests_M$hessian_inv[indices_M,indices_M]
      indices_F <- which(names(cfg2$ests_F$opt$par) %in% params)
      beta_F <- matrix(cfg2$ests_F$opt$par[indices_F])
      Sigma_F <- cfg2$ests_F$hessian_inv[indices_F,indices_F]
    }
    
    hr <- function(x, sex, log=F) {
      if (sex=="M") {
        est <- c(A(x) %*% beta_M)
        se <- c(sqrt(A(x) %*% Sigma_M %*% t(A(x))))
      } else if (sex=="F") {
        est <- c(A(x) %*% beta_F)
        se <- c(sqrt(A(x) %*% Sigma_F %*% t(A(x))))
      }
      if (log) {
        return(est + c(0,-1.96,1.96)*se)
      } else {
        return(exp(est + c(0,-1.96,1.96)*se))
      }
    }
    
    hr2 <- function(j, w_1, sex, log=F) {
      if (sex=="M") {
        est <- c(A(j, w_1) %*% beta_M)
        se <- c(sqrt(A(j, w_1) %*% Sigma_M %*% t(A(j, w_1))))
      } else if (sex=="F") {
        est <- c(A(j, w_1) %*% beta_F)
        se <- c(sqrt(A(j, w_1) %*% Sigma_F %*% t(A(j, w_1))))
      }
      if (log) {
        return(est + c(0,-1.96,1.96)*se)
      } else {
        return(exp(est + c(0,-1.96,1.96)*se))
      }
    }
    
    if (x_axis=="cal time") {
      x_grid <- seq(cfg2$w_start,cfg2$w_end,0.1)
      grid <- sapply(x_grid, function(x) { (x-cfg2$w_start+1)/10 })
      breaks <- seq(cfg2$w_start, cfg2$w_end, length.out=5)
    } else if (x_axis=="age") {
      x_grid <- seq(13,60,0.1)
      grid <- sapply(x_grid, function(x) { x / 100 })
      breaks <- seq(20,60, length.out=5)
    }
    
    if (plot_name=="HR_mort_hiv_cal") {
      
      log <- F
      df_plot <- data.frame(
        x = rep(x_grid,6),
        y = c(
          sapply(grid, function(j) { hr2(j=j, w_1=0.20, sex="M", log=log)[1] }),
          sapply(grid, function(j) { hr2(j=j, w_1=0.20, sex="F", log=log)[1] }),
          sapply(grid, function(j) { hr2(j=j, w_1=0.40, sex="M", log=log)[1] }),
          sapply(grid, function(j) { hr2(j=j, w_1=0.40, sex="F", log=log)[1] }),
          sapply(grid, function(j) { hr2(j=j, w_1=0.60, sex="M", log=log)[1] }),
          sapply(grid, function(j) { hr2(j=j, w_1=0.60, sex="F", log=log)[1] })
        ),
        ci_lo = c(
          sapply(grid, function(j) { hr2(j=j, w_1=0.20, sex="M", log=log)[2] }),
          sapply(grid, function(j) { hr2(j=j, w_1=0.20, sex="F", log=log)[2] }),
          sapply(grid, function(j) { hr2(j=j, w_1=0.40, sex="M", log=log)[2] }),
          sapply(grid, function(j) { hr2(j=j, w_1=0.40, sex="F", log=log)[2] }),
          sapply(grid, function(j) { hr2(j=j, w_1=0.60, sex="M", log=log)[2] }),
          sapply(grid, function(j) { hr2(j=j, w_1=0.60, sex="F", log=log)[2] })
        ),
        ci_up = c(
          sapply(grid, function(j) { hr2(j=j, w_1=0.20, sex="M", log=log)[3] }),
          sapply(grid, function(j) { hr2(j=j, w_1=0.20, sex="F", log=log)[3] }),
          sapply(grid, function(j) { hr2(j=j, w_1=0.40, sex="M", log=log)[3] }),
          sapply(grid, function(j) { hr2(j=j, w_1=0.40, sex="F", log=log)[3] }),
          sapply(grid, function(j) { hr2(j=j, w_1=0.60, sex="M", log=log)[3] }),
          sapply(grid, function(j) { hr2(j=j, w_1=0.60, sex="F", log=log)[3] })
        ),
        sex = rep(rep(c("Male", "Female"),3), each=length(x_grid)),
        age = rep(
          rep(c("Age: 20", "Age: 40", "Age: 60"), each=2),
          each=length(x_grid)
        )
      )
      
      plot <- ggplot(
        data = df_plot,
        aes(x=x, y=y)) +
        geom_line(color="forestgreen") +
        geom_ribbon(
          aes(ymin=ci_lo, ymax=ci_up),
          alpha = 0.2,
          fill = "forestgreen"
        ) +
        scale_x_continuous(
          breaks = breaks,
          minor_breaks = seq(cfg2$w_start, cfg2$w_end, 1)
        ) +
        scale_y_continuous(
          breaks = c(1:10),
          minor_breaks = NULL
          # labels = c(c(1:10), rep("", 5))
        ) +
        coord_cartesian(ylim=c(1,10)) +
        facet_grid(rows=dplyr::vars(sex), cols=dplyr::vars(age)) +
        labs(
          x = "Year",
          y = "Hazard Ratio",
          title = title
        )
      
      ggsave(
        filename = paste0("../Figures + Tables/", cfg2$d, " ", plot_name,
                          " - model ", cfg2$m, ".pdf"),
        plot=plot, device="pdf", width=9, height=6
      )
      
      # Save plot for paper
      plot_paper <- plot + labs(title=NULL)
      ggsave(
        filename = paste0("../Figures + Tables/", cfg2$d, " paper_", plot_name,
                          " - model ", cfg2$m, ".pdf"),
        plot=plot_paper, device="pdf", width=9, height=6
      )
      
    } else {
      # df_plot <- data.frame(
      #   x = x_grid,
      #   y = sapply(grid, function(x) { hr(x, log=log)[1] }),
      #   ci_lo = sapply(grid, function(x) { hr(x, log=log)[2] }),
      #   ci_up = sapply(grid, function(x) { hr(x, log=log)[3] })
      # )
      df_plot <- data.frame(
        x = rep(x_grid,4),
        y = c(
          sapply(grid, function(x) { hr(x, sex="M", log=T)[1] }),
          sapply(grid, function(x) { hr(x, sex="M", log=F)[1] }),
          sapply(grid, function(x) { hr(x, sex="F", log=T)[1] }),
          sapply(grid, function(x) { hr(x, sex="F", log=F)[1] })
        ),
        ci_lo = c(
          sapply(grid, function(x) { hr(x, sex="M", log=T)[2] }),
          sapply(grid, function(x) { hr(x, sex="M", log=F)[2] }),
          sapply(grid, function(x) { hr(x, sex="F", log=T)[2] }),
          sapply(grid, function(x) { hr(x, sex="F", log=F)[2] })
        ),
        ci_up = c(
          sapply(grid, function(x) { hr(x, sex="M", log=T)[3] }),
          sapply(grid, function(x) { hr(x, sex="M", log=F)[3] }),
          sapply(grid, function(x) { hr(x, sex="F", log=T)[3] }),
          sapply(grid, function(x) { hr(x, sex="F", log=F)[3] })
        ),
        sex = rep(c("Male", "Female"), each=2*length(x_grid)),
        log = rep(c("log(HR)", "HR", "log(HR)", "HR"), each=length(x_grid))
      )
      
      plot <- ggplot(
        data = df_plot,
        aes(x=x, y=y)) +
        geom_line(color="forestgreen") +
        geom_ribbon(
          aes(ymin=ci_lo, ymax=ci_up),
          alpha = 0.2,
          fill = "forestgreen"
        ) +
        scale_x_continuous(breaks=breaks) +
        facet_grid(rows=dplyr::vars(log), cols=dplyr::vars(sex), scales="free") +
        labs(
          x = ifelse(x_axis=="cal time", "Year", "Age"),
          # y = ifelse(log, "Log hazard ratio", "Hazard Ratio"),
          y = "Hazard Ratio",
          title = title
        )
      # assign(paste0("plot_", log), plot)
      
    }
    
    # plot <- ggpubr::ggarrange(plot_FALSE, plot_TRUE)
    ggsave(
      filename = paste0("../Figures + Tables/", cfg2$d, " ", plot_name,
                        " - model ", cfg2$m, ".pdf"),
      plot=plot, device="pdf", width=8, height=4 # 6x4
    )
    
  }
  
}



##################################################.
##### TAB: Table of probabilities for UNAIDS #####
##################################################.

if (F) {
  
  df_tab <- data.frame(
    "Year" = integer(),
    "Age" = double(),
    "Sex" = double(),
    "Rate" = double(),
    "Status" = character(),
    "ci_lo" = double(),
    "ci_up" = double()
  )
  
  for (year in c(cfg2$w_start:cfg2$w_end)) {
    for (age in seq(15,60,5)) {
    # for (age in c(15:59)) { # !!!!!
      for (sex in c(0,1)) {
        for (status in c("HIV-", "HIV+")) {
          year_sc <- year - cfg2$w_start + 1
          type <- paste0("mort (", status, ")")
          rate <- 1000 * prob(
            type=type, m=cfg2$m, j=year_sc, w_1=sex, w_2=age, which="est"
          )
          ci_lo <- 1000 * prob(
            type=type, m=cfg2$m, j=year_sc, w_1=sex, w_2=age, which="ci_lo"
          )
          ci_up <- 1000 * prob(
            type=type, m=cfg2$m, j=year_sc, w_1=sex, w_2=age, which="ci_up"
          )
          df_tab[nrow(df_tab)+1,] <- list(year,age,sex,rate,status,ci_lo,ci_up)
        }
      }
    }
  }
  
  utils::write.table(
    x = df_tab,
    file = paste0("../Figures + Tables/", cfg2$d, " temp_table_of_vals.csv"),
    sep = ",",
    row.names = FALSE
  )
  
}
