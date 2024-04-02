#######################################.
##### VIZ: All params (one model) #####
#######################################.

sim <- readRDS("sim.rds")

# p_names should match those used by negloglik()
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
      a_x=-2.231, g_x1=-0.4977, g_x2=-0.9101, a_y=-6.340, g_y1=0.5996,
      g_y2=2.410, g_y3=2.814, g_y4=7.096, g_y5=6.013,
      beta_x=1.295, beta_z=1.286, t_x=-1.366, t_y=-0.7141,
      a_s=-2.098, t_s=0.3321, g_s1=-0.8771, g_s2=0.8316
    )
  } else if (m==11) {
    p <- list(
      a_x=-8.0742, g_x1=-0.7309, g_x2=4.2919, g_x3=1.0879, g_x4=6.5164,
      g_x5=-8.508, t_x=-1.3492, a_s=-3.7151, g_s1=-0.4807, g_s2=1.1155,
      t_s=1.1718, beta_x=1.1894, beta_z=0.8276, a_y=-8.657, g_y1=0.4471,
      g_y2=4.4427, g_y3=5.073, g_y4=11.3514, g_y5=5.022, t_y=-0.7018
    )
  } else if (m==12) {
    p <- list(
      a_x=-4.6674897, g_x1=-0.7444126, g_x2=3.6579576, g_x3=-1.9702465, g_x4=-0.8989252,
      g_x5=-6.6903666, t_x=-1.1693100, a_s=-3.2990485, g_s1=-0.6594199, g_s2=0.8443497,
      t_s=1.0513203, beta_x=-0.7316351, beta_z=0.4216804, a_y=-5.6547364, g_y1=0.2733265,
      g_y2=1.7590798, g_y3=2.8024851, g_y4=6.1807478, g_y5=2.6495728, t_y=-0.6332008
    )
  } else if (m==13) {
    p <- list(
      a_x=-6.425, g_x1=-0.6614, g_x2=3.873, g_x3=0.3610, g_x4=0.4009,
      g_x5=-8.318, t_x1=-1.543, t_x2=-0.4422, t_x3=0.2004, t_x4=-1.972,
      a_s=-3.281, g_s1=-0.5525, g_s2=0.9998, t_s=0.9649, beta_x=1.009,
      beta_z=0.7339, a_y=-6.016, g_y1=0.3872, g_y2=1.862, g_y3=3.037,
      g_y4=6.715, g_y5=3.417, t_y=-0.6883
    )
  }
  
  j <- j/10
  w_2 <- w_2/100
  
  if (type=="sero") {
    
    if (m<11) {
      prob <- exp2(p$a_x + p$t_x*j + p$g_x1*w_1 + p$g_x2*w_2)
    } else if (m==11) {
      prob <- exp2(
        p$a_x + p$t_x*j + p$g_x1*w_1 + p$g_x2*b1(w_2,1) + p$g_x3*b1(w_2,2) +
          p$g_x4*b1(w_2,3) + p$g_x5*b1(w_2,4)
      )
    } else if (m==12) {
      prob <- exp2(
        p$a_x + p$t_x*j + p$g_x1*w_1 + p$g_x2*b2(w_2,1) + p$g_x3*b2(w_2,2) +
          p$g_x4*b2(w_2,3) + p$g_x5*b2(w_2,4)
      )
    } else if (m==13) {
      prob <- exp2(
        p$a_x + p$t_x1*b4(j,1) + p$t_x2*b4(j,2) + p$t_x3*b4(j,3) +
          p$t_x4*b4(j,4) + p$g_x1*w_1 + p$g_x2*b2(w_2,1) + p$g_x3*b2(w_2,2) +
          p$g_x4*b2(w_2,3) + p$g_x5*b2(w_2,4)
      )
    }
    
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
    } else if (m %in% c(9,11)) {
      prob <- exp2(
        p$a_y + p$t_y*j + p$g_y1*w_1 + p$g_y2*b1(w_2,1) + p$g_y3*b1(w_2,2) +
          p$g_y4*b1(w_2,3) + p$g_y5*b1(w_2,4) + p$beta_x*x + p$beta_z*z
      )
    } else if (m %in% c(12,13)) {
      prob <- exp2(
        p$a_y + p$t_y*j + p$g_y1*w_1 + p$g_y2*b3(w_2,1) + p$g_y3*b3(w_2,2) +
          p$g_y4*b3(w_2,3) + p$g_y5*b3(w_2,4) + p$beta_x*x + p$beta_z*z
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
    grid <- seq(13,90,0.1)
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

m <- 13
b1 <- construct_basis("age (0-100), 4DF")
b2 <- construct_basis("age (13,20,30,60,90)")
b3 <- construct_basis("age (13,30,60,75,90)")
b4 <- construct_basis("year (00,05,10,15,20)")

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
 