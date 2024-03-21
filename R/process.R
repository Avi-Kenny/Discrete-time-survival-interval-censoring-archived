#######################################.
##### VIZ: All params (one model) #####
#######################################.

sim <- readRDS("sim.rds")

# p_names should match those used by negloglik_miss
p_names <- c("a_x", "g_x1", "g_x2", "a_y", "g_y1", "g_y2", "beta_x", "beta_z",
            "t_x", "t_y", "a_s", "t_s", "g_s1", "g_s2")
# true_vals <- log(c(0.005,1.3,1.2,0.003,1.2,1.1,1.5,0.7,1,1,0.05,1,2,1.5)) # Monthly
true_vals <- log(c(0.05,1.3,1.2,0.03,1.2,1.1,1.5,0.7,1,1,0.05,1,2,1.5)) # Yearly

# prm_sim <- sim$levels$params[[1]]
# true_vals2 <- c(prm_sim$a_x, prm_sim$g_x[1], prm_sim$g_x[2], prm_sim$a_y,
#                prm_sim$g_y[1], prm_sim$g_y[2], prm_sim$beta_x, prm_sim$beta_z,
#                # prm_sim$t_x, prm_sim$t_y,
#                prm_sim$a_s, prm_sim$t_s,
#                prm_sim$g_s[1], prm_sim$g_s[2])

v <- paste0("lik_M_",p_names,"_est")
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

# Summary stats
summ_mean <- summ_sd <- summ_cov <- list()
for (i in c(1:length(p_names))) {
  p <- p_names[i]
  summ_mean[[i]] <- list(stat="mean",
                         name=paste0(p,"__sd_est"),
                         x=paste0("lik_M_",p,"_se"))
  summ_sd[[i]] <- list(stat="sd", name=paste0(p,"__sd_actual"),
                       x=paste0("lik_M_",p,"_est"))
  summ_cov[[i]] <- list(stat="coverage", name=paste0(p,"__cov"),
                        truth=true_vals[i],
                        estimate=paste0("lik_M_",p,"_est"),
                        se=paste0("lik_M_",p,"_se"), na.rm=T)
}
summ <- do.call(SimEngine::summarize, c(list(sim), summ_mean, summ_sd, summ_cov))
l_id <- 1
summ2 <- summ[summ$level_id==l_id]
df_results <- data.frame(
  "param" = character(),
  "sd_est" = double(),
  "sd_actual" = double(),
  "coverage" = double()
)
for (p in p_names) {
  df_results[nrow(df_results)+1,] <- c(
    p,
    round(summ[[paste0(p,"__sd_actual")]], 3),
    round(summ[[paste0(p,"__sd_est")]], 3),
    round(summ[[paste0(p,"__cov")]], 3)
  )
}
print(df_results)



#######################################################.
##### Function for plotting modeled probabilities #####
#######################################################.

#' Return modeled probability
#' @param type One of c("sero", "init", "mort (HIV-)", "mort (HIV+ART-)",
#'     "mort (HIV+ART+)")
#' @param j Calendar time (in years, starting at 1, where 1 represents the
#'     first year in the dataset)
#' @param w_1 Sex (0=female, 1=male)
#' @param w_2 Age (in years)
#' @param m An integer representing the model version number
#' @return Numeric probability
prob <- function(type, m, j, w_1, w_2) {
  if (m==7) {
    p <- list(
      a_x=-3.5607, g_x1=-0.3244, g_x2=-0.2809, a_y=-5.7446, g_y1=0.3544,
      g_y2=4.4057, beta_x=1.8096, beta_z=1.8153, t_x=-0.786, t_y=-0.7826,
      a_s=-2.87, t_s=0.6349, g_s1=-0.3768, g_s2=0.6409
    )
  } else if (m==8) {
    p <- list(
      a_x=-2.2421, g_x1=-0.5886, g_x2=-0.9116, a_y=-6.3142, g_y1=0.3968,
      g_y2=7.1492, g_y3=-3.9706, g_y4=1.9161, beta_x=1.6961, beta_z=1.4211,
      t_x=-1.7858, t_y=-0.7475, a_s=-5.148, t_s=1.9954, g_s1=-0.4421, g_s2=0.344
    )
  } else if (m==9) {
    p <- list(
      a_x=-2.273, g_x1=-1.131, g_x2=-2.875, a_y=-3.350, g_y1=0.2217,
      g_y2=0.6883, g_y3=0.7507, g_y4=0.004464, g_y5=3.886,
      beta_x=1.376, beta_z=0.8172, t_x=-1.332, t_y=-0.6786,
      a_s=-3.273, t_s=0.7709, g_s1=-0.4208, g_s2=0.9081
    )
  }
  j <- j/10
  w_2 <- w_2/100
  
  if (type=="sero") {
    
    prob <- exp2(p$a_x + p$t_x*j + p$g_x1*w_1 + p$g_x2*w_2)
    
  } else if (type=="init") {
    
    prob <- exp2(p$a_s + p$t_s*j + p$g_s1*w_1 + p$g_s2*w_2)
    
  } else {
    
    if (type=="mort (HIV-)") { x <- 0; z <- 0;}
    if (type=="mort (HIV+ART-)") { x <- 1; z <- 0;}
    if (type=="mort (HIV+ART+)") { x <- 0; z <- 1;}
    
    if (m==7) {
      prob <- exp2(
        p$a_y + p$t_y*j + p$g_y1*w_1 + p$g_y2*w_2 + p$beta_x*x + p$beta_z*z
      )
    } else if (m==8) {
      prob <- exp2(
        p$a_y + p$t_y*j + p$g_y1*w_1 + p$g_y2*w_2 + p$g_y3*w_2^2 +
          p$g_y4*w_2^3 + p$beta_x*x + p$beta_z*z
      )
    } else if (m==9) {
      prob <- exp2(
        p$a_y + p$t_y*j + p$g_y1*w_1 + p$g_y2*b(w_2,1) + p$g_y3*b(w_2,2) +
          p$g_y4*b(w_2,3) + p$g_y5*b(w_2,4) + p$beta_x*x + p$beta_z*z
      )
    }
    
  }

  return(prob)
  
}

#' Return plot of modeled probabilities
#' @param x One of c("Year", "Age"); the variable to go on the X-axis
#' @param type One of c("sero", "init", "mort")
#' @param m An integer representing the model version number
#' @return ggplot2 object
plot_mod <- function(x_axis, type, m) {
  
  if (x_axis=="Age") {
    grid <- seq(13,100,0.1)
    color <- "Year"
    plot_data <- data.frame(
      x = rep(grid,6),
      Probability = c(
        sapply(grid, function(w_2) { prob(type, m, j=0, w_1=0, w_2) }),
        sapply(grid, function(w_2) { prob(type, m, j=0, w_1=1, w_2) }),
        sapply(grid, function(w_2) { prob(type, m, j=10, w_1=0, w_2) }),
        sapply(grid, function(w_2) { prob(type, m, j=10, w_1=1, w_2) }),
        sapply(grid, function(w_2) { prob(type, m, j=20, w_1=0, w_2) }),
        sapply(grid, function(w_2) { prob(type, m, j=20, w_1=1, w_2) })
      ),
      Sex = rep(rep(c("Female", "Male"),3), each=length(grid)),
      color = rep(c("2000","2010","2020"), each=2*length(grid))
    )
  } else if (x_axis=="Year") {
    grid <- seq(0,22,0.01)
    color <- "Age"
    plot_data <- data.frame(
      x = rep(grid,6) + 1999,
      Probability = c(
        sapply(grid, function(j) { prob(type, m, j, w_1=0, w_2=20) }),
        sapply(grid, function(j) { prob(type, m, j, w_1=1, w_2=20) }),
        sapply(grid, function(j) { prob(type, m, j, w_1=0, w_2=35) }),
        sapply(grid, function(j) { prob(type, m, j, w_1=1, w_2=35) }),
        sapply(grid, function(j) { prob(type, m, j, w_1=0, w_2=50) }),
        sapply(grid, function(j) { prob(type, m, j, w_1=1, w_2=50) })
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
    labs(title=title, x=x_axis, color=color) +
    theme(plot.background = element_rect(color="black"))
  
  return(plot)
  
}



###############################################.
##### VIZ: Plotting modeled probabilities #####
###############################################.

m <- 9

# Construct spline basis
if (m==9) {
  b <- construct_basis(m=9, s=1)
}

# Seroconversion prob as a function of age
p01 <- plot_mod(x_axis="Age", type="sero", m=m)

# Seroconversion prob as a function of calendar time
p02 <- plot_mod(x_axis="Year", type="sero", m=m)

# HIV+ initial status as a function of age
p03 <- plot_mod(x_axis="Age", type="init", m=m)

# HIV+ initial status as a function of calendar time
p04 <- plot_mod(x_axis="Year", type="init", m=m)

# Mortality prob as a function of age
p05 <- plot_mod(x_axis="Age", type="mort (HIV-)", m=m)
p06 <- plot_mod(x_axis="Age", type="mort (HIV+ART-)", m=m)
# p07 <- plot_mod(x_axis="Age", type="mort (HIV+ART+)", m=m)

# Mortality prob as a function of calendar time
p08 <- plot_mod(x_axis="Year", type="mort (HIV-)", m=m)
p09 <- plot_mod(x_axis="Year", type="mort (HIV+ART-)", m=m)
# p10 <- plot_mod(x_axis="Year", type="mort (HIV+ART+)", m=m)

print(ggpubr::ggarrange(p01, p02)) # Export 10"x5"
print(ggpubr::ggarrange(p03, p04)) # Export 10"x5"
print(ggpubr::ggarrange(p05, p08, p06, p09)) # Export 10"x10"
 