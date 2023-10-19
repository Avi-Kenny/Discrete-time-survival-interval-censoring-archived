
##################################################################.
##### Generate dataset, run summary stats, and fit Cox model #####
##################################################################.

# Runtimes (seconds)
# Note: need to update this
# 2 (n=500), 5 (n=1000), 18 (n=2000), 65 (n=4000)

# Set up vectors
# T_minus <- T_plus <- c()
c1 <- c2 <- c3 <- c4 <- c()
num_pos <- num_art <- num_ev <- num_ev_pos <- num_ev_pos2 <- num_ev_pos3 <- c()
g_y_1 <- g_y_2 <- beta_x <- c()

# Generate datasets and calculate summary stats
n_reps <- 5
for (i in c(1:n_reps)) {
  
  # Generate data
  dat <- generate_data(
    n = 500,
    max_time = 70,
    params = list(
      a_x=log(0.005), a_y=log(0.003), a_v=log(0.7), a_z=log(0.01),
      g_x=c(log(1.3),log(1.2)), g_y=c(log(1.2),log(1.1)),
      g_v=c(log(1.2),log(1.1)), g_z=c(log(1.2),log(1.1)),
      beta_x=log(1.5), beta_z=log(0.7)
    )
  )
  
  # Summary stats
  num_pos[i] <- sum(dplyr::summarize(group_by(dat,id), x=max(x))$x) # Number seroconverted
  num_art[i] <- sum(dplyr::summarize(group_by(dat,id), z=max(z))$z) # Number ART+
  num_ev[i] <- sum(dat$y) # Number of events
  num_ev_pos[i] <- sum(dat$y*dat$x) # Number of events among HIV+
  num_ev_pos2[i] <- sum(dat$y*dat$x*(1-dat$z)) # Number of events among HIV+ART-
  num_ev_pos3[i] <- sum(dat$y*dat$x*dat$z) # Number of events among HIV+ART+
  
  # T_minus[i] <- attr(dat, "T_minus")
  # T_plus[i] <- attr(dat, "T_plus")
  c1[i] <- sum(attr(dat, "case")==1)
  c2[i] <- sum(attr(dat, "case")==2)
  c3[i] <- sum(attr(dat, "case")==3)
  c4[i] <- sum(attr(dat, "case")==4)
  
  # Run Cox PH model
  model <- coxph(
    Surv(t_start, t_end, y) ~ w_1 + w_2 + x + cluster(id),
    data = dat
  )
  g_y_1[i] <- model$coefficients[1]
  g_y_2[i] <- model$coefficients[2]
  beta_x[i] <- model$coefficients[3]
  
}

# Visualize summary stats
# Export PDF 5" x 10"
stats <- c("# case 1", "# case 2", "# case 3", "# case 4", "# seroconverted",
           "# ART+", "# of events", "# of events among HIV+",
           "# of events among HIV+ART-", "# of events among HIV+ART+",
           "exp(gamma_y[1])", "exp(gamma_y[2])", "exp(beta_x)")
df_plot <- data.frame(
  x = c(c1,c2,c3,c4,num_pos,num_art,num_ev,num_ev_pos,num_ev_pos2,num_ev_pos3,
        exp(g_y_1),exp(g_y_2),exp(beta_x)),
  y = rep(0,n_reps),length(stats),
  which = rep(factor(stats, levels=stats), each=n_reps)
)
ggplot(df_plot, aes(x=x, y=y, color=which)) +
  geom_jitter(width=0, height=1, alpha=0.5, size=3) +
  facet_wrap(~which, scales="free_x", ncol=1, strip.position="left") +
  labs(y=NULL) +
  ylim(-2,2) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text.y.left = element_text(angle=0),
        legend.position="none")
