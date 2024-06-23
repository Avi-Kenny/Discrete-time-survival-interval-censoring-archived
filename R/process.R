#######################################.
##### VIZ: All params (one model) #####
#######################################.

if (F) {
  
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
  
}




#######################################################.
##### Function for plotting modeled probabilities #####
#######################################################.

#' Return modeled probability
#' @param type One of c("sero", "init", "mort (HIV-)", "mort (HIV+)",
#'     "mort (HIV+ART-)", "mort (HIV+ART+)")
#' @param j Calendar time (in years, starting at 1, where 1 represents the
#'     first year in the dataset)
#' @param w_1 Sex (0=female, 1=male)
#' @param w_2 Age (in years)
#' @param m An integer representing the model version number
#' @return Numeric probability
prob <- function(type, m, j, w_1, w_2) {
  
  if (m==7) {
    p <- list(a_x=-3.5607, g_x1=-0.3244, g_x2=-0.2809, a_y=-5.7446, g_y1=0.3544, g_y2=4.4057, beta_x=1.8096, beta_z=1.8153, t_x=-0.786, t_y=-0.7826, a_s=-2.87, t_s=0.6349, g_s1=-0.3768, g_s2=0.6409)
  } else if (m==8) {
    p <- list(a_x=-1.629, g_x1=-0.643, g_x2=-2.100, a_y=-5.945, g_y1=0.382, g_y2=4.633, g_y3=-0.742, g_y4=1.280, beta_x=1.099, beta_z=1.046, t_x=-1.739, t_y=-0.720, a_s=-2.272, t_s=0.249, g_s1=-0.485, g_s2=1.511)
  } else if (m==9) {
    p <- list(a_x=-2.231, g_x1=-0.4977, g_x2=-0.9101, a_y=-6.340, g_y1=0.5996, g_y2=2.410, g_y3=2.814, g_y4=7.096, g_y5=6.013, beta_x=1.295, beta_z=1.286, t_x=-1.366, t_y=-0.7141, a_s=-2.098, t_s=0.3321, g_s1=-0.8771, g_s2=0.8316)
  } else if (m==11) {
    p <- list(a_x=-8.0742, g_x1=-0.7309, g_x2=4.2919, g_x3=1.0879, g_x4=6.5164, g_x5=-8.508, t_x=-1.3492, a_s=-3.7151, g_s1=-0.4807, g_s2=1.1155, t_s=1.1718, beta_x=1.1894, beta_z=0.8276, a_y=-8.657, g_y1=0.4471, g_y2=4.4427, g_y3=5.073, g_y4=11.3514, g_y5=5.022, t_y=-0.7018)
  } else if (m==12) {
    p <- list(a_x=-4.667, g_x1=-0.744, g_x2=3.657, g_x3=-1.970, g_x4=-0.898, g_x5=-6.690, t_x=-1.169, a_s=-3.299, g_s1=-0.659, g_s2=0.844, t_s=1.051, beta_x=-0.731, beta_z=0.421, a_y=-5.654, g_y1=0.273, g_y2=1.759, g_y3=2.802, g_y4=6.180, g_y5=2.649, t_y=-0.633)
  } else if (m==13) {
    p <- list(a_x=-6.425, g_x1=-0.6614, g_x2=3.873, g_x3=0.3610, g_x4=0.4009, g_x5=-8.318, t_x1=-1.543, t_x2=-0.4422, t_x3=0.2004, t_x4=-1.972, a_s=-3.281, g_s1=-0.5525, g_s2=0.9998, t_s=0.9649, beta_x=1.009, beta_z=0.7339, a_y=-6.016, g_y1=0.3872, g_y2=1.862, g_y3=3.037, g_y4=6.715, g_y5=3.417, t_y=-0.6883)
  } else if (m==14) {
    p <- list(a_x=-6.356, g_x1=-0.7200, g_x2=3.587, g_x3=0.8982, g_x4=1.234, g_x5=-7.876, t_x1=-2.196, t_x2=-0.3959, t_x3=-0.007392, t_x4=-2.254, a_s=-3.066, g_s1=-0.6221, g_s2=0.8340, t_s=0.8966, beta_x=0.9409, beta_z=0.7236, a_y=-6.504, g_y1=0.4055, g_y2=1.991, g_y3=3.096, g_y4=6.936, g_y5=3.469, t_y1=-0.2982, t_y2=-0.8745, t_y3=-0.5772, t_y4=-1.072)
  } else if (m==15) {
    p <- list(a_x=-6.67, g_x1=4.86, g_x2=1.33, g_x3=0.871, g_x4=-6.53, g_x5=3.66, g_x6=1.07, g_x7=2.60, g_x8=-7.59, t_x1=-1.82, t_x2=-0.728, t_x3=-0.593, t_x4=-2.06, a_s=-3.08, g_s1=-0.742, g_s2=0.753, t_s=0.936, beta_x=1.12, beta_z=1.00, a_y=-6.56, g_y1=0.431, g_y2=1.95, g_y3=3.29, g_y4=7.18, g_y5=3.60, t_y1=-0.319, t_y2=-1.00, t_y3=-0.934, t_y4=-1.12)
  } else if (m==16) {
    p <- list(a_x=-6.66, g_x1=5.02, g_x2=1.27, g_x3=0.508, g_x4=-6.98, g_x5=3.77, g_x6=0.790, g_x7=2.33, g_x8=-8.08, t_x1=-1.76, t_x2=-0.824, t_x3=-0.990, t_x4=-2.07, a_s=-2.43, g_s1=-0.716, g_s2=0.755, t_s1=0.670, t_s2=0.545, t_s3=0.509, t_s4=1.76, beta_x=1.16, beta_z=0.979, a_y=-6.52, g_y1=0.434, g_y2=1.94, g_y3=3.25, g_y4=7.04, g_y5=3.62, t_y1=-0.274, t_y2=-1.02, t_y3=-0.964, t_y4=-1.10)
  } else if (m==17) {
    p <- list(a_x=-6.9521, g_x1=4.2808, g_x2=-0.1476, g_x3=1.4616, g_x4=-5.7527, g_x5=2.4626, g_x6=0.1639, g_x7=3.6117, g_x8=-7.3954, t_x1=-1.093, t_x2=-0.888, t_x3=-1.3534, t_x4=-1.4756, a_s=-2.8117, g_s1=-0.3086, g_s2=0.8134, g_s3=-0.1185, g_s4=3.2074, g_s5=-1.4434, beta_x=1.369, beta_z=1.1501, a_y=-6.5127, g_y1=0.4162, g_y2=1.7359, g_y3=3.1651, g_y4=6.5974, g_y5=4.0223, t_y1=-0.1488, t_y2=-0.8869, t_y3=-0.6838, t_y4=-1.0119)
  } else if (m==18) {
    p <- list(a_x=-8.725, g_x1=3.414, g_x2=-0.561, g_x3=0.040, g_x4=-6.618, g_x5=-0.058, g_x6=0.498, g_x7=6.534, g_x8=-4.481, t_x1=-1.146, t_x2=-1.745, t_x3=-0.895, t_x4=-2.134, a_s=-2.914, g_s1=-0.435, g_s2=1.418, g_s3=-0.779, g_s4=3.834, g_s5=-2.681, beta_x1=1.813, beta_x2=-1.234, beta_z1=3.826,beta_z2=-3.966, a_y=-7.504, g_y1=0.598, g_y2=2.118, g_y3=3.794, g_y4=6.900, g_y5=4.007, t_y1=0.303, t_y2=-0.594, t_y3=-1.262, t_y4=-1.038)
  } else if (m==19) {
    p <- list(a_x=-8.567, g_x1=2.507, g_x2=-0.537, g_x3=-1.743, g_x4=-8.788, g_x5=-1.952, g_x6=0.833, g_x7=6.446, g_x8=-4.131, t_x1=0.772, t_x2=-0.797, t_x3=0.015, t_x4=-1.276, a_s=-2.914, g_s1=-0.455, g_s2=1.383, g_s3=-0.784, g_s4=3.879, g_s5=-2.659, beta_x1=1.818, beta_x2=-1.261, beta_z1=3.829, beta_z2=-4.027, a_y=-7.330, g_y1=0.568, g_y2=1.865, g_y3=3.614, g_y4=6.156, g_y5=3.917, t_y1=-0.142, t_y2=0.151, t_y3=0.308, t_y4=-0.183)
  } else if (m==20) {
    p <- list(a_x=-8.020, g_x1=2.288, g_x2=-1.682, g_x3=2.458, g_x4=-2.482, g_x5=-1.755, g_x6=1.293, g_x7=7.403, g_x8=-3.002, t_x1=-0.367, t_x2=-1.639, t_x3=-0.807, t_x4=-1.199, a_s=-3.130, g_s1=-0.435, g_s2=1.359, g_s3=-0.208, g_s4=3.285, g_s5=-4.575, beta_x1=0.024, beta_x2=-0.338, beta_x3=2.098, beta_x4=-1.204, beta_z1=0.549, beta_z2=-0.332, beta_z3=2.451, beta_z4=-1.659, a_y=-6.999, g_y1=0.582, g_y2=2.339, g_y3=3.775, g_y4=7.476, g_y5=3.943, t_y1=-0.835, t_y2=0.060, t_y3=-1.958, t_y4=-0.761)
  } else if (m==21) {
    p <- list(a_x=-8.106, g_x1=2.437, g_x2=-1.533, g_x3=2.540, g_x4=-2.320, g_x5=-1.679, g_x6=1.008, g_x7=7.416, g_x8=-2.923, t_x1=-0.284, t_x2=-1.655, t_x3=-0.782, t_x4=-1.090, a_s=-3.130, g_s1=-0.443, g_s2=1.369, g_s3=-0.218, g_s4=3.286, g_s5=-4.564, beta_x1=0.119, beta_x2=-0.316, beta_x3=2.065, beta_x4=-1.207, a_y=-6.998, g_y1=0.559, g_y2=2.351, g_y3=3.749, g_y4=7.466, g_y5=3.935, t_y1=-0.824, t_y2=0.082, t_y3=-1.909, t_y4=-0.774) # Before removing 75+ year olds
    p <- list(a_x=-8.255, g_x1=3.370, g_x2=-3.520, g_x3=2.301, g_x4=-0.685, g_x5=-1.084, g_x6=0.012, g_x7=7.442, g_x8=-1.853, t_x1=0.325, t_x2=-2.542, t_x3=-0.288, t_x4=-0.498, a_s=-3.068, g_s1=-0.466, g_s2=1.208, g_s3=0.093, g_s4=3.367, g_s5=-4.529, beta_x1=0.369, beta_x2=-0.390, beta_x3=1.893, beta_x4=-1.047, a_y=-6.839, g_y1=0.493, g_y2=2.066, g_y3=4.414, g_y4=7.743, g_y5=3.549, t_y1=-1.175, t_y2=-0.015, t_y3=-2.257, t_y4=-0.998) # After removing 75+ year olds
  } else if (m==22) {
    stop()
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
    } else if (m %in% c(19:21)) {
      prob <- icll(
        p$a_x + p$t_x1*b5(j,1) + p$t_x2*b5(j,2) + p$t_x3*b5(j,3) +
          p$t_x4*b5(j,4) + w_1*(
            p$g_x1*b2(w_2,1) + p$g_x2*b2(w_2,2) +
              p$g_x3*b2(w_2,3) + p$g_x4*b2(w_2,4)
          ) + (1-w_1)*(
            p$g_x5*b2(w_2,1) + p$g_x6*b2(w_2,2) +
              p$g_x7*b2(w_2,3) + p$g_x8*b2(w_2,4)
          )
      )
    } else if (m==22) {
      prob <- icll(
        p$a_x + p$t_x1*b5(j,1) + p$t_x2*b5(j,2) + p$t_x3*b5(j,3) +
          p$t_x4*b5(j,4) + w_1*(
            p$g_x1*b6(w_2,1) + p$g_x2*b6(w_2,2) +
              p$g_x3*b6(w_2,3) + p$g_x4*b6(w_2,4)
          ) + (1-w_1)*(
            p$g_x5*b6(w_2,1) + p$g_x6*b6(w_2,2) +
              p$g_x7*b6(w_2,3) + p$g_x8*b6(w_2,4)
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
    } else if (m %in% c(17:21)) {
      prob <- icll(
        p$a_s + p$g_s1*w_1 + p$g_s2*b3(w_2,1) + p$g_s3*b3(w_2,2) + 
          p$g_s4*b3(w_2,3) + p$g_s5*b3(w_2,4)
      )
    } else if (m==22) {
      prob <- icll(
        p$a_s + p$g_s1*w_1 + p$g_s2*b6(w_2,1) + p$g_s3*b6(w_2,2) + 
          p$g_s4*b6(w_2,3) + p$g_s5*b6(w_2,4)
      )
    }
    
  } else {
    
    if (type=="mort (HIV-)") { x <- 0; z <- 0;}
    if (type=="mort (HIV+)") { x <- 1; z <- NA;}
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
    } else if (m==19) {
      prob <- icll(
        p$a_y + p$t_y1*b5(j,1) + p$t_y2*b5(j,2) + p$t_y3*b5(j,3) +
          p$t_y4*b5(j,4) + p$g_y1*w_1 + p$g_y2*b3(w_2,1) + p$g_y3*b3(w_2,2) +
          p$g_y4*b3(w_2,3) + p$g_y5*b3(w_2,4) +
          (p$beta_x1+p$beta_x2*j)*x*(1-z) + (p$beta_z1+p$beta_z2*j)*x*z
      )
    } else if (m==20) {
      prob <- icll(
        p$a_y + p$t_y1*b5(j,1) + p$t_y2*b5(j,2) + p$t_y3*b5(j,3) +
          p$t_y4*b5(j,4) + p$g_y1*w_1 + p$g_y2*b3(w_2,1) + p$g_y3*b3(w_2,2) +
          p$g_y4*b3(w_2,3) + p$g_y5*b3(w_2,4) +
          x*(1-z)*(
            p$beta_x1*b5(j,1) + p$beta_x2*b5(j,2) + p$beta_x3*b5(j,3) +
              p$beta_x4*b5(j,4)
          ) + x*z*(
            p$beta_z1*b5(j,1) + p$beta_z2*b5(j,2) + p$beta_z3*b5(j,3) +
              p$beta_z4*b5(j,4)
          )
      )
    } else if (m==21) {
      prob <- icll(
        p$a_y + p$t_y1*b5(j,1) + p$t_y2*b5(j,2) + p$t_y3*b5(j,3) +
          p$t_y4*b5(j,4) + p$g_y1*w_1 + p$g_y2*b3(w_2,1) + p$g_y3*b3(w_2,2) +
          p$g_y4*b3(w_2,3) + p$g_y5*b3(w_2,4) +
          x*(p$beta_x1*b5(j,1) + p$beta_x2*b5(j,2) + p$beta_x3*b5(j,3) +
               p$beta_x4*b5(j,4))
      )
    } else if (m==22) {
      prob <- icll(
        p$a_y + p$t_y1*b5(j,1) + p$t_y2*b5(j,2) + p$t_y3*b5(j,3) +
          p$t_y4*b5(j,4) + p$g_y1*w_1 + p$g_y2*b6(w_2,1) + p$g_y3*b6(w_2,2) +
          p$g_y4*b6(w_2,3) + p$g_y5*b6(w_2,4) +
          x*(p$beta_x1*b5(j,1) + p$beta_x2*b5(j,2) + p$beta_x3*b5(j,3) +
               p$beta_x4*b5(j,4))
      )
    }
    
  }

  return(prob)
  
}

#' Return plot of modeled probabilities
#' @param x One of c("Year", "Age"); the variable to go on the X-axis
#' @param type One of c("sero", "init", "mort (HIV-)", "mort (HIV+)",
#'     "mort (HIV+ART-)", "mort (HIV+ART+)")
#' @param m An integer representing the model version number
#' @param y_max Maximum Y value for the plot
#' @return ggplot2 object
plot_mod <- function(x_axis, type, m, w_start, y_max) {
  
  if (w_start==2000) {
    j_vals <- c(1,11,21) # This was formerly incorrectly set to c(0,10,20)
    j_labs <- c("2000","2010","2020")
    grid <- seq(2000,2023,0.01) %>% (function(x) { x-(w_start-1) })
  } else if (w_start==2010) {
    j_vals <- c(1,6,11)
    j_labs <- c("2010","2015","2020")
    grid <- seq(2010,2023,0.01) %>% (function(x) { x-(w_start-1) })
  }
  
  if (x_axis=="Age") {
    grid <- seq(13,75,0.1)
    color <- "Year"
    plot_data <- data.frame(
      x = rep(grid,6),
      Probability = c(
        sapply(grid, function(w_2) { prob(type, m, j=j_vals[1], w_1=0, w_2) }),
        sapply(grid, function(w_2) { prob(type, m, j=j_vals[1], w_1=1, w_2) }),
        sapply(grid, function(w_2) { prob(type, m, j=j_vals[2], w_1=0, w_2) }),
        sapply(grid, function(w_2) { prob(type, m, j=j_vals[2], w_1=1, w_2) }),
        sapply(grid, function(w_2) { prob(type, m, j=j_vals[3], w_1=0, w_2) }),
        sapply(grid, function(w_2) { prob(type, m, j=j_vals[3], w_1=1, w_2) })
      ),
      Sex = rep(rep(c("Female", "Male"),3), each=length(grid)),
      color = rep(j_labs, each=2*length(grid))
    )
  } else if (x_axis=="Year") {
    color <- "Age"
    plot_data <- data.frame(
      x = rep(grid,6) + (w_start-1),
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
    ylim(0, y_max) +
    labs(title=title, x=x_axis, color=color) +
    theme(plot.background = element_rect(color="black"))
  
  return(plot)
  
}



###############################################.
##### VIZ: Plotting modeled probabilities #####
###############################################.

m <- 21
w_start <- 2010
hivart <- "HIV" # One of c("HIV", "HIV+ART")
y_max <- c(0.06, 0.8, 0.2, 0.05)
b1 <- construct_basis("age (0-100), 4DF")
b2 <- construct_basis("age (13,20,30,60,90)")
b3 <- construct_basis("age (13,30,60,75,90)")
b4 <- construct_basis("year (00,05,10,15,20)", window_start=w_start)
b5 <- construct_basis("year (10,13,17,20,23)", window_start=w_start)
b6 <- construct_basis("age (15,30,45,60,75)", window_start=w_start)

# Seroconversion prob as a function of age
p01 <- plot_mod(x_axis="Age", type="sero", m=m, w_start=w_start, y_max=y_max[1])

# Seroconversion prob as a function of calendar time
p02 <- plot_mod(x_axis="Year", type="sero", m=m, w_start=w_start, y_max=y_max[1])

# HIV+ initial status as a function of age
p03 <- plot_mod(x_axis="Age", type="init", m=m, w_start=w_start, y_max=y_max[2])

# HIV+ initial status as a function of calendar time
p04 <- plot_mod(x_axis="Year", type="init", m=m, w_start=w_start, y_max=y_max[2])

# HIV only
if (hivart=="HIV") {
  
  # Mortality prob as a function of age
  p05 <- plot_mod(x_axis="Age", type="mort (HIV-)", m=m, w_start=w_start, y_max=y_max[3])
  p06 <- plot_mod(x_axis="Age", type="mort (HIV+)", m=m, w_start=w_start, y_max=y_max[3])

  # Mortality prob as a function of calendar time
  p08 <- plot_mod(x_axis="Year", type="mort (HIV-)", m=m, w_start=w_start, y_max=y_max[4])
  p09 <- plot_mod(x_axis="Year", type="mort (HIV+)", m=m, w_start=w_start, y_max=y_max[4])

} else if (hivart=="HIV+ART") {
  
  # Mortality prob as a function of age
  p05 <- plot_mod(x_axis="Age", type="mort (HIV-)", m=m, w_start=w_start, y_max=y_max[3])
  p06 <- plot_mod(x_axis="Age", type="mort (HIV+ART-)", m=m, w_start=w_start, y_max=y_max[3])
  p07 <- plot_mod(x_axis="Age", type="mort (HIV+ART+)", m=m, w_start=w_start, y_max=y_max[3])
  
  # Mortality prob as a function of calendar time
  p08 <- plot_mod(x_axis="Year", type="mort (HIV-)", m=m, w_start=w_start, y_max=y_max[4])
  p09 <- plot_mod(x_axis="Year", type="mort (HIV+ART-)", m=m, w_start=w_start, y_max=y_max[4])
  p10 <- plot_mod(x_axis="Year", type="mort (HIV+ART+)", m=m, w_start=w_start, y_max=y_max[4])
  
}

print(ggpubr::ggarrange(p01, p02)) # Export 10"x5"
print(ggpubr::ggarrange(p03, p04)) # Export 10"x5"
if (hivart=="HIV") {
  print(ggpubr::ggarrange(p05, p08, p06, p09)) # Export 10"x10"
} else if (hivart=="HIV+ART") {
  print(ggpubr::ggarrange(p05, p06, p07, p08, p09, p10, ncol=3, nrow=2)) # Export 15"x10"
}



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
