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
      # a_x=-2.0049, g_x1=-0.6199, g_x2=-1.1939, a_y=-6.1569, g_y1=0.4559,   # Exp2 link
      # g_y2=4.4967, g_y3=0.1729, g_y4=0.6047, beta_x=1.5029, beta_z=1.0327, # Exp2 link
      # t_x=-1.0738, t_y=-0.7204, a_s=-3.0524, t_s=0.9542, g_s1=-0.679,      # Exp2 link
      # g_s2=0.5649                                                          # Exp2 link
      a_x=-1.629, g_x1=-0.643, g_x2=-2.100, a_y=-5.945, g_y1=0.382,          # ICLL link
      g_y2=4.633, g_y3=-0.742, g_y4=1.280, beta_x=1.099, beta_z=1.046,       # ICLL link
      t_x=-1.739, t_y=-0.720, a_s=-2.272, t_s=0.249, g_s1=-0.485,            # ICLL link
      g_s2=1.511                                                             # ICLL link
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
  } else if (m==14) {
    p <- list(
      a_x=-6.356, g_x1=-0.7200, g_x2=3.587, g_x3=0.8982, g_x4=1.234,
      g_x5=-7.876, t_x1=-2.196, t_x2=-0.3959, t_x3=-0.007392, t_x4=-2.254,
      a_s=-3.066, g_s1=-0.6221, g_s2=0.8340, t_s=0.8966, beta_x=0.9409,
      beta_z=0.7236, a_y=-6.504, g_y1=0.4055, g_y2=1.991, g_y3=3.096,
      g_y4=6.936, g_y5=3.469, t_y1=-0.2982, t_y2=-0.8745, t_y3=-0.5772,
      t_y4=-1.072
    )
  } else if (m==15) {
    p <- list(
      a_x=-6.67, g_x1=4.86, g_x2=1.33, g_x3=0.871, g_x4=-6.53, g_x5=3.66,
      g_x6=1.07, g_x7=2.60, g_x8=-7.59, t_x1=-1.82, t_x2=-0.728, t_x3=-0.593,
      t_x4=-2.06, a_s=-3.08, g_s1=-0.742, g_s2=0.753, t_s=0.936, beta_x=1.12,
      beta_z=1.00, a_y=-6.56, g_y1=0.431, g_y2=1.95, g_y3=3.29, g_y4=7.18,
      g_y5=3.60, t_y1=-0.319, t_y2=-1.00, t_y3=-0.934, t_y4=-1.12
    )
  } else if (m==16) {
    p <- list(
      a_x=-6.66, g_x1=5.02, g_x2=1.27, g_x3=0.508, g_x4=-6.98, g_x5=3.77,
      g_x6=0.790, g_x7=2.33, g_x8=-8.08, t_x1=-1.76, t_x2=-0.824, t_x3=-0.990,
      t_x4=-2.07, a_s=-2.43, g_s1=-0.716, g_s2=0.755, t_s1=0.670, t_s2=0.545,
      t_s3=0.509, t_s4=1.76, beta_x=1.16, beta_z=0.979, a_y=-6.52, g_y1=0.434,
      g_y2=1.94, g_y3=3.25, g_y4=7.04, g_y5=3.62, t_y1=-0.274, t_y2=-1.02,
      t_y3=-0.964, t_y4=-1.10
    )
  } else if (m==17) {
    p <- list(
      a_x=-6.9521, g_x1=4.2808, g_x2=-0.1476, g_x3=1.4616, g_x4=-5.7527,
      g_x5=2.4626, g_x6=0.1639, g_x7=3.6117, g_x8=-7.3954, t_x1=-1.093,
      t_x2=-0.888, t_x3=-1.3534, t_x4=-1.4756, a_s=-2.8117, g_s1=-0.3086,
      g_s2=0.8134, g_s3=-0.1185, g_s4=3.2074, g_s5=-1.4434, beta_x=1.369,
      beta_z=1.1501, a_y=-6.5127, g_y1=0.4162, g_y2=1.7359, g_y3=3.1651,
      g_y4=6.5974, g_y5=4.0223, t_y1=-0.1488, t_y2=-0.8869, t_y3=-0.6838,
      t_y4=-1.0119
    )
  } else if (m==18) {
    p <- list(
      a_x=-7.855, g_x1=4.791, g_x2=-0.935, g_x3=0.844, g_x4=-8.116, g_x5=2.804,
      g_x6=-2.155, g_x7=5.066, g_x8=-5.616, t_x1=-1.084, t_x2=-0.764,
      t_x3=-1.505, t_x4=-2.598, a_s=-3.043, g_s1=-0.281, g_s2=1.184, g_s3=0.219,
      g_s4=3.380, g_s5=-3.239, beta_x1=2.636, beta_x2=-1.051, beta_z1=3.983,
      beta_z2=-1.870, a_y=-7.221, g_y1=0.411, g_y2=1.692, g_y3=3.440,
      g_y4=6.288, g_y5=3.826, t_y1=0.413, t_y2=0.030, t_y3=0.530, t_y4=-0.235
    )
  }
  
  j <- j/10
  w_2 <- w_2/100
  
  if (type=="sero") {
    
    if (m<11) {
      prob <- icll(p$a_x + p$t_x*j + p$g_x1*w_1 + p$g_x2*w_2)
    } else if (m==11) {
      prob <- icll(
        p$a_x + p$t_x*j + p$g_x1*w_1 + p$g_x2*b1(w_2,1) + p$g_x3*b1(w_2,2) +
          p$g_x4*b1(w_2,3) + p$g_x5*b1(w_2,4)
      )
    } else if (m==12) {
      prob <- icll(
        p$a_x + p$t_x*j + p$g_x1*w_1 + p$g_x2*b2(w_2,1) + p$g_x3*b2(w_2,2) +
          p$g_x4*b2(w_2,3) + p$g_x5*b2(w_2,4)
      )
    } else if (m %in% c(13,14)) {
      prob <- icll(
        p$a_x + p$t_x1*b4(j,1) + p$t_x2*b4(j,2) + p$t_x3*b4(j,3) +
          p$t_x4*b4(j,4) + p$g_x1*w_1 + p$g_x2*b2(w_2,1) + p$g_x3*b2(w_2,2) +
          p$g_x4*b2(w_2,3) + p$g_x5*b2(w_2,4)
      )
    } else if (m %in% c(15:18)) {
      prob <- icll(
        p$a_x + p$t_x1*b4(j,1) + p$t_x2*b4(j,2) + p$t_x3*b4(j,3) +
          p$t_x4*b4(j,4) + w_1*(
            p$g_x1*b2(w_2,1) + p$g_x2*b2(w_2,2) +
              p$g_x3*b2(w_2,3) + p$g_x4*b2(w_2,4)
          ) + (1-w_1)*(
            p$g_x5*b2(w_2,1) + p$g_x6*b2(w_2,2) +
              p$g_x7*b2(w_2,3) + p$g_x8*b2(w_2,4)
          )
      )
    }
    
  } else if (type=="init") {
    
    if (m<16) {
      prob <- icll(p$a_s + p$t_s*j + p$g_s1*w_1 + p$g_s2*w_2)
    } else if (m==16) {
      prob <- icll(
        p$a_s + p$t_s1*b4(j,1) + p$t_s2*b4(j,2) + p$t_s3*b4(j,3) +
          p$t_s4*b4(j,4) + p$g_s1*w_1 + p$g_s2*w_2
      )
    } else if (m %in% c(17,18)) {
      prob <- icll(
        p$a_s + p$g_s1*w_1 + p$g_s2*b3(w_2,1) + p$g_s3*b3(w_2,2) + 
          p$g_s4*b3(w_2,3) + p$g_s5*b3(w_2,4)
      )
    }
    
  } else {
    
    if (type=="mort (HIV-)") { x <- 0; z <- 0;}
    if (type=="mort (HIV+ART-)") { x <- 1; z <- 0;}
    if (type=="mort (HIV+ART+)") { x <- 1; z <- 1;}
    
    if (m==7) {
      prob <- icll(
        p$a_y + p$t_y*j + p$g_y1*w_1 + p$g_y2*w_2 + p$beta_x*x + p$beta_z*z
      )
    } else if (m==8) {
      prob <- icll(
        p$a_y + p$t_y*j + p$g_y1*w_1 + p$g_y2*w_2 + p$g_y3*w_2^2 +
          p$g_y4*w_2^3 + p$beta_x*x*(1-z) + p$beta_z*x*z
      )
    } else if (m %in% c(9,11)) {
      prob <- icll(
        p$a_y + p$t_y*j + p$g_y1*w_1 + p$g_y2*b1(w_2,1) + p$g_y3*b1(w_2,2) +
          p$g_y4*b1(w_2,3) + p$g_y5*b1(w_2,4) + p$beta_x*x*(1-z) + p$beta_z*z*z
      )
    } else if (m %in% c(12,13)) {
      prob <- icll(
        p$a_y + p$t_y*j + p$g_y1*w_1 + p$g_y2*b3(w_2,1) + p$g_y3*b3(w_2,2) +
          p$g_y4*b3(w_2,3) + p$g_y5*b3(w_2,4) + p$beta_x*x*(1-z) + p$beta_z*x*z
      )
    } else if (m %in% c(14:17)) {
      prob <- icll(
        p$a_y + p$t_y1*b4(j,1) + p$t_y2*b4(j,2) + p$t_y3*b4(j,3) +
          p$t_y4*b4(j,4) + p$g_y1*w_1 + p$g_y2*b3(w_2,1) + p$g_y3*b3(w_2,2) +
          p$g_y4*b3(w_2,3) + p$g_y5*b3(w_2,4) + p$beta_x*x*(1-z) + p$beta_z*x*z
      )
    } else if (m==18) {
      prob <- icll(
        p$a_y + p$t_y1*b4(j,1) + p$t_y2*b4(j,2) + p$t_y3*b4(j,3) +
          p$t_y4*b4(j,4) + p$g_y1*w_1 + p$g_y2*b3(w_2,1) + p$g_y3*b3(w_2,2) +
          p$g_y4*b3(w_2,3) + p$g_y5*b3(w_2,4) +
          (p$beta_x1+p$beta_x2*j)*x*(1-z) + (p$beta_z1+p$beta_z2*j)*x*z
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

m <- 18
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
p07 <- plot_mod(x_axis="Age", type="mort (HIV+ART+)", m=m)

# Mortality prob as a function of calendar time
p08 <- plot_mod(x_axis="Year", type="mort (HIV-)", m=m)
p09 <- plot_mod(x_axis="Year", type="mort (HIV+ART-)", m=m)
p10 <- plot_mod(x_axis="Year", type="mort (HIV+ART+)", m=m)

print(ggpubr::ggarrange(p01, p02)) # Export 10"x5"
print(ggpubr::ggarrange(p03, p04)) # Export 10"x5"
print(ggpubr::ggarrange(p05, p08, p06, p09)) # Export 10"x10"
print(ggpubr::ggarrange(p05, p06, p07, p08, p09, p10, ncol=3, nrow=2)) # Export 10"x15"



# !!!!! 3D plot
if (F) {
  library(plotly)
  len <- 50
  age_vec <- seq(13,90, length.out=len)
  year_vec <- seq(0,22, length.out=len)
  prob_mtx <- matrix(NA, nrow=len, ncol=len)
  for (row in c(1:len)) {
    for (col in c(1:len)) {
      age <- age_vec[row]
      year <- year_vec[col]
      prob_mtx[row,col] <- prob(type="sero", m=15, j=year, w_1=0, w_2=age)
    }
  }
  year_vec_display <- year_vec + 1999
  fig <- plot_ly(x=age_vec, y=year_vec_display, z=prob_mtx)
  fig <- fig %>% add_surface()
  fig
}
