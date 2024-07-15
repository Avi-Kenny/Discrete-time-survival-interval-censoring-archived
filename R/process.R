#################.
##### Setup #####
#################.

cfg2 <- list(
  process_sims = F,
  process_analysis = T,
  # m = 24,
  m = 25,
  w_start = 2010,
  ests = readRDS("objs/ests_25_full_20240715.rds")
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

# Get current date
cfg2$d <- format(Sys.time(), "%Y-%m-%d")



#####################################################.
##### VIZ (simulations): all params (one model) #####
#####################################################.

if (cfg2$process_sims) {
  
  sim <- readRDS(paste0("SimEngine.out/sim_20231017 (v8, added time trend to b",
                        "aseline model).rds"))
  
  # p_names should match those used by negloglik()
  p_names <- c("a_x", "g_x1", "g_x2", "a_y", "g_y1", "g_y2", "beta_x", "beta_z",
               "t_x", "t_y", "a_s", "t_s", "g_s1", "g_s2")
  true_vals <- log(c(0.05,1.3,1.2,0.03,1.2,1.1,1.5,0.7,1,1,0.05,1,2,1.5)) # Yearly
  
  # prm_sim <- sim$levels$params[[1]]
  # true_vals2 <- c(prm_sim$a_x, prm_sim$g_x[1], prm_sim$g_x[2], prm_sim$a_y,
  #                prm_sim$g_y[1], prm_sim$g_y[2], prm_sim$beta_x, prm_sim$beta_z,
  #                # prm_sim$t_x, prm_sim$t_y,
  #                prm_sim$a_s, prm_sim$t_s,
  #                prm_sim$g_s[1], prm_sim$g_s[2])
  
  v <- paste0("lik_M_",p_names,"_est")
  r <- dplyr::filter(sim$results, params=="70pct testing")
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
    labs(y=NULL, title="70% testing") +
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
  utils::write.table(
    x = df_results,
    file = "../Figures + Tables/sims_sd_and_coverage.csv",
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
#' @param w_1 Sex (0=female, 1=male)
#' @param w_2 Age (in years)
#' @param m An integer representing the model version number
#' @param which One of c("est", "ci_lo", "ci_up")
#' @return Numeric probability
prob <- function(type, m, j, w_1, w_2, which="est") {
  
  j <- j/10
  w_2 <- w_2/100
  
  if (m %in% c(0:23)) {
    
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
      p <- list(a_x=-8.255, g_x1=3.370, g_x2=-3.520, g_x3=2.301, g_x4=-0.685, g_x5=-1.084, g_x6=0.012, g_x7=7.442, g_x8=-1.853, t_x1=0.325, t_x2=-2.542, t_x3=-0.288, t_x4=-0.498, a_s=-3.068, g_s1=-0.466, g_s2=1.208, g_s3=0.093, g_s4=3.367, g_s5=-4.529, beta_x1=0.369, beta_x2=-0.390, beta_x3=1.893, beta_x4=-1.047, a_y=-6.839, g_y1=0.493, g_y2=2.066, g_y3=4.414, g_y4=7.743, g_y5=3.549, t_y1=-1.175, t_y2=-0.015, t_y3=-2.257, t_y4=-0.998) # After removing 75+ year olds
    } else if (m==22) {
      p <- list(a_x=-5.974, g_x1=-0.153, g_x2=-0.892, g_x3=1.491, g_x4=1.431, g_x5=-6.335, g_x6=-0.520, g_x7=4.633, g_x8=-1.664, t_x1=-0.122, t_x2=-2.452, t_x3=-0.125, t_x4=-0.539, a_s=-2.965, g_s1=-0.516, g_s2=3.010, g_s3=0.457, g_s4=3.360, g_s5=-2.050, beta_x1=0.854, beta_x2=-0.511, beta_x3=1.923, beta_x4=-0.401, a_y=-7.149, g_y1=0.554, g_y2=2.471, g_y3=2.540, g_y4=6.911, g_y5=3.670, t_y1=-1.368, t_y2=0.105, t_y3=-2.186, t_y4=-1.376)
    } else if (m==23) {
      p <- list(g_x1=2.320, g_x2=-1.314, g_x3=-3.709, g_x4=-8.004, g_x5=7.655, g_x6=2.610, g_x7=-8.068, g_x8=0.399, g_x9=-3.308, g_x10=-2.569, t_x1=-1.323, t_x2=-5.764, t_x3=-11.607, t_x4=-0.073, a_s=-3.095, g_s1=-0.514, g_s2=3.104, g_s3=0.499, g_s4=4.633, g_s5=-0.516, beta_x1=1.054, beta_x2=0.747, beta_x3=-0.315, beta_x4=2.192, beta_x5=-0.820, a_y=-8.450, g_y1=0.527, g_y2=2.843, g_y3=3.086, g_y4=8.093, g_y5=3.746, t_y1=-0.530, t_y2=0.400, t_y3=-0.635, t_y4=-0.973)
    }    
    
  }

  if (type=="sero") {
    
    if (m %in% c(0:23)) {
      
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
      } else if (m==23) {
        prob <- icll(
          p$t_x1*b5(j,1) + p$t_x2*b5(j,2) + p$t_x3*b5(j,3) + p$t_x4*b5(j,4) +
            w_1*(
              p$g_x1*b8(w_2,1) + p$g_x2*b8(w_2,2) + p$g_x3*b8(w_2,3) +
                p$g_x4*b8(w_2,4) + p$g_x5*b8(w_2,5)
            ) + (1-w_1)*(
              p$g_x6*b8(w_2,1) + p$g_x7*b8(w_2,2) + p$g_x8*b8(w_2,3) +
                p$g_x9*b8(w_2,4) + p$g_x10*b8(w_2,5)
            )
        )
      }
      
    }
    if (m==24) {
      
      A <- function(j,w_2) { t(matrix(c(
        1, b5(j,1), b5(j,2), b5(j,3), b5(j,4), w_1*b6(w_2,1), w_1*b6(w_2,2),
        w_1*b6(w_2,3), w_1*b6(w_2,4), (1-w_1)*b6(w_2,1), (1-w_1)*b6(w_2,2),
        (1-w_1)*b6(w_2,3), (1-w_1)*b6(w_2,4)
      ))) }
      p2 <- c("a_x", "t_x1", "t_x2", "t_x3", "t_x4", "g_x1", "g_x2", "g_x3",
              "g_x4", "g_x5", "g_x6", "g_x7", "g_x8")
      
    } else if (m==25) {
      
      A <- function(j,w_2) { t(matrix(c(
        1, b5(j,1), b5(j,2), b5(j,3), b5(j,4), w_1*b9(w_2,1), w_1*b9(w_2,2),
        w_1*b9(w_2,3), w_1*b9(w_2,4), (1-w_1)*b9(w_2,1), (1-w_1)*b9(w_2,2),
        (1-w_1)*b9(w_2,3), (1-w_1)*b9(w_2,4)
      ))) }
      p2 <- c("a_x", "t_x1", "t_x2", "t_x3", "t_x4", "g_x1", "g_x2", "g_x3",
              "g_x4", "g_x5", "g_x6", "g_x7", "g_x8")
      
    }
    
  } else if (type=="init") {
    
    if (m %in% c(0:23)) {
      
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
      } else if (m %in% c(22:23)) {
        prob <- icll(
          p$a_s + p$g_s1*w_1 + p$g_s2*b6(w_2,1) + p$g_s3*b6(w_2,2) + 
            p$g_s4*b6(w_2,3) + p$g_s5*b6(w_2,4)
        )
      }      
      
    }
    if (m==24) {
      A <- function(j,w_2) { t(matrix(c(
        1, w_1, b6(w_2,1), b6(w_2,2), b6(w_2,3), b6(w_2,4)
      ))) }
      p2 <- c("a_s", "g_s1", "g_s2", "g_s3", "g_s4", "g_s5")
    } else if (m==25) {
      A <- function(j,w_2) { t(matrix(c(
        1, w_1, b9(w_2,1), b9(w_2,2), b9(w_2,3), b9(w_2,4)
      ))) }
      p2 <- c("a_s", "g_s1", "g_s2", "g_s3", "g_s4", "g_s5")
    }
    
  } else {
    
    if (type=="mort (HIV-)") { x <- 0 }
    if (type=="mort (HIV+)") { x <- 1 }
    # if (type=="mort (HIV-)") { x <- 0; z <- 0;}
    # if (type=="mort (HIV+)") { x <- 1; z <- NA;}
    # if (type=="mort (HIV+ART-)") { x <- 1; z <- 0;}
    # if (type=="mort (HIV+ART+)") { x <- 1; z <- 1;}
    
    if (m %in% c(0:23)) {
      
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
      } else if (m==23) {
        A <- function(j,w_2) { t(matrix(c(
          x*b7(j,1), x*b7(j,2), x*b7(j,3), x*b7(j,4), x*b7(j,5), 1, w_1,
          b6(w_2,1), b6(w_2,2), b6(w_2,3), b6(w_2,4), b5(j,1), b5(j,2), b5(j,3),
          b5(j,4)
        ))) }
        p2 <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5", "a_y",
                "g_y1", "g_y2", "g_y3", "g_y4", "g_y5", "t_y1", "t_y2", "t_y3",
                "t_y4")
      }
      
    }
    if (m==24) {
      A <- function(j,w_2) { t(matrix(c(
        x*b7(j,1), x*b7(j,2), x*b7(j,3), x*b7(j,4), x*b7(j,5), 1, w_1,
        b6(w_2,1), b6(w_2,2), b6(w_2,3), b6(w_2,4), b5(j,1), b5(j,2), b5(j,3),
        b5(j,4)
      ))) }
      p2 <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5", "a_y",
              "g_y1", "g_y2", "g_y3", "g_y4", "g_y5", "t_y1", "t_y2", "t_y3",
              "t_y4")
    } else if (m==25) {
      A <- function(j,w_2) { t(matrix(c(
        x*b7(j,1), x*b7(j,2), x*b7(j,3), x*b7(j,4), x*b7(j,5), 1, w_1,
        b9(w_2,1), b9(w_2,2), b9(w_2,3), b9(w_2,4), b5(j,1), b5(j,2), b5(j,3),
        b5(j,4)
      ))) }
      p2 <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5", "a_y",
              "g_y1", "g_y2", "g_y3", "g_y4", "g_y5", "t_y1", "t_y2", "t_y3",
              "t_y4")
    }
    
  }
  
  indices <- as.numeric(sapply(p2, function(p) {
    which(names(cfg2$ests$opt$par)==p)
  }))
  beta <- matrix(cfg2$ests$opt$par[indices])
  Sigma <- cfg2$ests$hessian_inv[indices,indices]
  if (which=="est") {
    fac <- 0
  } else if (which=="ci_lo") {
    fac <- -1.96
  } else if (which=="ci_up") {
    fac <- 1.96
  }
  est <- c(A(j,w_2) %*% beta)
  se <- c(sqrt(A(j,w_2) %*% Sigma %*% t(A(j,w_2))))
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
#' @param y_max Maximum Y value for the plot
#' @return ggplot2 object
plot_mod <- function(x_axis, type, m, w_start, y_max) {
  
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
    grid <- seq(w_start,2023,0.01) %>% (function(x) { x-(w_start-1) })
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
    coord_cartesian(ylim=c(0,y_max)) +
    labs(title=title, x=x_axis, color=color) +
    theme(plot.background = element_rect(color="black"))
  
  return(plot)
  
}



#' Return plot of modeled mortality probabilities (HIV- vs. HIV+) with CIs
#' @param x_axis One of c("Year", "Age"); the variable to go on the X-axis
#' @param m An integer representing the model version number
#' @param w_start An integer representing the window start calendar year
#' @param y_max Maximum Y value for the plot
#' @param title Boolean; if F, title is suppressed
#' @return ggplot2 object
plot_mort3 <- function(x_axis, m, w_start, y_max=NA, title=T) {
  
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
        prob(type=type, m=m, j=outer, w_1=sex, w_2=inner, which=which)
      })
    }
    
  } else if (x_axis=="Year") {
    
    grid <- seq(w_start,2023,0.01) %>% (function(x) { x-(w_start-1) })
    breaks <- seq(w_start,2022, length.out=5)
    color <- "Age"
    outer <- c(20,35,50)
    prob2 <- function(type, which, outer) {
      1000 * sapply(grid, function(inner) {
        prob(type=type, m=m, j=inner, w_1=sex, w_2=outer, which=which)
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
      prob(type="sero", m=m, j=outer, w_1=sex, w_2=inner, which=which)
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
  p01 <- plot_mod(x_axis="Age", type="sero", m=cfg2$m, w_start=cfg2$w_start, y_max=y_max[1])
  
  # Seroconversion prob as a function of calendar time
  p02 <- plot_mod(x_axis="Year", type="sero", m=cfg2$m, w_start=cfg2$w_start, y_max=y_max[1])
  
  # HIV+ initial status as a function of age
  p03 <- plot_mod(x_axis="Age", type="init", m=cfg2$m, w_start=cfg2$w_start, y_max=y_max[2])
  
  # HIV+ initial status as a function of calendar time
  p04 <- plot_mod(x_axis="Year", type="init", m=cfg2$m, w_start=cfg2$w_start, y_max=y_max[2])
  
  # Mortality prob as a function of age
  p05 <- plot_mod(x_axis="Age", type="mort (HIV-)", m=cfg2$m, w_start=cfg2$w_start, y_max=y_max[3])
  p06 <- plot_mod(x_axis="Age", type="mort (HIV+)", m=cfg2$m, w_start=cfg2$w_start, y_max=y_max[3])
  
  # Mortality prob as a function of calendar time
  p08 <- plot_mod(x_axis="Year", type="mort (HIV-)", m=cfg2$m, w_start=cfg2$w_start, y_max=y_max[4])
  p09 <- plot_mod(x_axis="Year", type="mort (HIV+)", m=cfg2$m, w_start=cfg2$w_start, y_max=y_max[4])
  
  # Old ART plots
  if (F) {
    # Mortality prob as a function of age
    p05 <- plot_mod(x_axis="Age", type="mort (HIV-)", m=cfg2$m, w_start=cfg2$w_start, y_max=y_max[3])
    p06 <- plot_mod(x_axis="Age", type="mort (HIV+ART-)", m=cfg2$m, w_start=cfg2$w_start, y_max=y_max[3])
    p07 <- plot_mod(x_axis="Age", type="mort (HIV+ART+)", m=cfg2$m, w_start=cfg2$w_start, y_max=y_max[3])
    
    # Mortality prob as a function of calendar time
    p08 <- plot_mod(x_axis="Year", type="mort (HIV-)", m=cfg2$m, w_start=cfg2$w_start, y_max=y_max[4])
    p09 <- plot_mod(x_axis="Year", type="mort (HIV+ART-)", m=cfg2$m, w_start=cfg2$w_start, y_max=y_max[4])
    p10 <- plot_mod(x_axis="Year", type="mort (HIV+ART+)", m=cfg2$m, w_start=cfg2$w_start, y_max=y_max[4])
    
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
                        y_max=65, title=F)
  ggsave(
    filename = paste0("../Figures + Tables/", cfg2$d, " p5 (mort CIs, by year)",
                      " - model ", cfg2$m, ".pdf"),
    plot = plot_05, device="pdf", width=8, height=5
  )
  
  # Mortality by age, with CIs
  plot_06 <- plot_mort3(x_axis="Age", m=cfg2$m, w_start=cfg2$w_start,
                        y_max=120, title=F)
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

if (cfg2$process_analysis) {
  
  dat <- readRDS("../Data/dat.rds")

  w_1_mean <- mean(dat$w_1)
  b9_1_mean <- mean(sapply(dat$w_2, function(w_2) { b9(w_2,1) }))
  b9_2_mean <- mean(sapply(dat$w_2, function(w_2) { b9(w_2,2) }))
  b9_3_mean <- mean(sapply(dat$w_2, function(w_2) { b9(w_2,3) }))
  b9_4_mean <- mean(sapply(dat$w_2, function(w_2) { b9(w_2,4) }))
  
  prob_m <- function(type, m, j, which="est") {
    
    j <- j/10
    
    if (type=="mort (HIV-)") { x <- 0; z <- NA; }
    if (type=="mort (HIV+)") { x <- 1; z <- NA; }
    
    if (m==25) {
      A <- function(j) { t(matrix(c(
        x*b7(j,1), x*b7(j,2), x*b7(j,3), x*b7(j,4), x*b7(j,5), 1, w_1_mean,
        b9_1_mean, b9_2_mean, b9_3_mean, b9_4_mean, b5(j,1), b5(j,2), b5(j,3),
        b5(j,4)
      ))) }
      p2 <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5", "a_y",
              "g_y1", "g_y2", "g_y3", "g_y4", "g_y5", "t_y1", "t_y2", "t_y3",
              "t_y4")
    }
    
    indices <- as.numeric(sapply(p2, function(p) {
      which(names(cfg2$ests$opt$par)==p)
    }))
    beta <- matrix(cfg2$ests$opt$par[indices])
    Sigma <- cfg2$ests$hessian_inv[indices,indices]
    if (which=="est") {
      fac <- 0
    } else if (which=="ci_lo") {
      fac <- -1.96
    } else if (which=="ci_up") {
      fac <- 1.96
    }
    est <- c(A(j) %*% beta)
    se <- c(sqrt(A(j) %*% Sigma %*% t(A(j))))
    prob <- icll(est+fac*se)
    
    return(prob)
    
  }
  
  grid <- seq(cfg2$w_start,2023,0.01) %>% (function(x) { x-(cfg2$w_start-1) })
  breaks_x <- seq(cfg2$w_start,2022, length.out=5)
  breaks_y <- seq(0,20,2)
  prob2_m <- function(type, which, outer) {
    1000 * sapply(grid, function(j) {
      prob_m(type=type, m=cfg2$m, j=j, which=which)
    })
  }
  
  plot_data <- data.frame(
    x = rep(grid,2) + (cfg2$w_start-1),
    Rate = c(prob2_m(type="mort (HIV-)", which="est", outer=o),
             prob2_m(type="mort (HIV+)", which="est", outer=o)),
    ci_lo = c(prob2_m(type="mort (HIV-)", which="ci_lo", outer=o),
              prob2_m(type="mort (HIV+)", which="ci_lo", outer=o)),
    ci_up = c(prob2_m(type="mort (HIV-)", which="ci_up", outer=o),
              prob2_m(type="mort (HIV+)", which="ci_up", outer=o)),
    color = rep(c("HIV-","HIV+"), each=length(grid))
  )
  
  # y_max <- 25
  plot_08 <- ggplot(
    plot_data,
    aes(x=x, y=Rate, color=color)
  ) +
    geom_line() +
    geom_ribbon(
      aes(ymin=ci_lo, ymax=ci_up, fill=color),
      alpha = 0.2,
      linetype = "blank"
    ) +
    # coord_cartesian(ylim=c(0,y_max)) +
    labs(
      title = "Marginalized mortality rates by HIV status, over calendar time",
      x = "Year",
      color = "HIV Status",
      fill = "HIV Status",
      y = "Deaths per 1,000 person-years"
    ) +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks=breaks_x) +
    scale_y_continuous(breaks=breaks_y) +
    scale_color_manual(values=c("forestgreen", "#56B4E9")) +
    scale_fill_manual(values=c("forestgreen", "#56B4E9"))
  
  ggsave(
    filename = paste0("../Figures + Tables/", cfg2$d, " p8 (mort CIs, by year,",
                      " marginal) - model ", cfg2$m, ".pdf"),
    plot = plot_08, device="pdf", width=8, height=5
  )
  
}



##################################################.
##### VIZ: HR as a function of calendar time #####
##################################################.

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
  
  # Extract estimates and SEs
  plot_names <- c("HR_mort_hiv_cal",
                  "HR_mort_age",
                  "HR_mort_cal",
                  "HR_sero_cal",
                  "HR_sero_male_age",
                  "HR_sero_female_age",
                  "HR_init_age")
  for (plot_name in plot_names) {
    
    # Set graph-specific variables
    if (plot_name=="HR_mort_hiv_cal") {
      title <- "HR of mortality, HIV+ vs. HIV- individuals (calendar time)"
      x_axis <- "cal time"
      if (cfg2$m %in% c(23:25)) {
        params <- c("beta_x1", "beta_x2", "beta_x3", "beta_x4", "beta_x5")
        A <- A_b7
      }
    } else if (plot_name=="HR_mort_age") {
      title <- "HR of mortality (age)"
      x_axis <- "age"
      if (cfg2$m==24) {
        params <- c("g_y2", "g_y3", "g_y4", "g_y5")
        A <- A_b6
      } else if (cfg2$m==25) {
        params <- c("g_y2", "g_y3", "g_y4", "g_y5")
        A <- A_b9
      }
    } else if (plot_name=="HR_mort_cal") {
      title <- "HR of mortality (calendar time)"
      x_axis <- "cal time"
      if (cfg2$m %in% c(24:25)) {
        params <- c("t_y1", "t_y2", "t_y3", "t_y4")
        A <- A_b5
      }
    } else if (plot_name=="HR_sero_cal") {
      title <- "HR of seroconversion (calendar time)"
      if (cfg2$m %in% c(23:25)) {
        params <- c("t_x1", "t_x2", "t_x3", "t_x4")
        A <- A_b5
      }
      x_axis <- "cal time"
    } else if (plot_name=="HR_sero_male_age") {
      title <- "HR of seroconversion among males (age)"
      x_axis <- "age"
      if (cfg2$m==23) {
        params <- c("g_x1", "g_x2", "g_x3", "g_x4", "g_x5")
        A <- A_b8
      } else if (cfg2$m==24) {
        params <- c("g_x1", "g_x2", "g_x3", "g_x4")
        A <- A_b6
      } else if (cfg2$m==25) {
        params <- c("g_x1", "g_x2", "g_x3", "g_x4")
        A <- A_b9
      }
    } else if (plot_name=="HR_sero_female_age") {
      title <- "HR of seroconversion among females (age)"
      x_axis <- "age"
      if (cfg2$m==23) {
        params <- c("g_x6", "g_x7", "g_x8", "g_x9", "g_x10")
        A <- A_b8
      } else if (cfg2$m==24) {
        params <- c("g_x5", "g_x6", "g_x7", "g_x8")
        A <- A_b6
      } else if (cfg2$m==25) {
        params <- c("g_x5", "g_x6", "g_x7", "g_x8")
        A <- A_b9
      }
    } else if (plot_name=="HR_init_age") {
      title <- "HR of HIV+ initial status (age)"
      x_axis <- "age"
      if (cfg2$m==24) {
        params <- c("g_s2", "g_s3", "g_s4", "g_s5")
        A <- A_b6
      } else if (cfg2$m==25) {
        params <- c("g_s2", "g_s3", "g_s4", "g_s5")
        A <- A_b9
      }
    }
    
    indices <- which(names(cfg2$ests$opt$par) %in% params)
    beta <- matrix(cfg2$ests$opt$par[indices])
    Sigma <- cfg2$ests$hessian_inv[indices,indices]
    
    hr <- function(x, log=F) {
      est <- c(A(x) %*% beta)
      se <- c(sqrt(A(x) %*% Sigma %*% t(A(x))))
      if (log) {
        return(est + c(0,-1.96,1.96)*se)
      } else {
        return(exp(est + c(0,-1.96,1.96)*se))
      }
    }
    
    if (x_axis=="cal time") {
      x_grid <- seq(cfg2$w_start,2023,0.1)
      grid <- sapply(x_grid, function(x) { (x-cfg2$w_start+1)/10 })
    } else if (x_axis=="age") {
      x_grid <- seq(13,60,0.1)
      grid <- sapply(x_grid, function(x) { x / 100 })
    }
    
    for (log in c(FALSE, TRUE)) {
      
      df_plot <- data.frame(
        x = x_grid,
        y = sapply(grid, function(x) { hr(x, log=log)[1] }),
        ci_lo = sapply(grid, function(x) { hr(x, log=log)[2] }),
        ci_up = sapply(grid, function(x) { hr(x, log=log)[3] })
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
        labs(
          x = ifelse(x_axis=="cal time", "Year", "Age"),
          y = ifelse(log, "Log hazard ratio", "Hazard Ratio"),
          title = title
        )
      assign(paste0("plot_", log), plot)
      
      # Save plot for paper
      if (plot_name=="HR_mort_hiv_cal" && !log) {
        plot_paper <- plot + labs(title=NULL)
        ggsave(
          filename = paste0(
            "../Figures + Tables/", cfg2$d, " paper_", plot_name, ".pdf"
          ), plot=plot_paper, device="pdf", width=8, height=4 # 6x4
        )
      }
      
    }
    
    plot <- ggpubr::ggarrange(plot_FALSE, plot_TRUE)
    ggsave(
      filename = paste0("../Figures + Tables/", cfg2$d, " ", plot_name, ".pdf"),
      plot=plot, device="pdf", width=8, height=4 # 6x4
    )
    
  }
  
}
