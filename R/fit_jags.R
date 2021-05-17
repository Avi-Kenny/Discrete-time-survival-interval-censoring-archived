#' Fit JAGS model
#'
#' @param dat !!!!!
#' @param mcmc A list of the form list(n.adapt=1, n.iter=1, n.burn=1, thin=1,
#'     n.chains=1)
#' @return !!!!!

fit_jags <- function(dat, mcmc) {
  
  # Declare JAGS model (with U latent variable as random effect)
  jags_code_u <- quote("
    model {
      
      for (i in 1:I) {
        
        u[i] ~ dnorm(0, 1/(sigma^2))
        
        for (j in 1:J[i]) {
          
          x[i,j+1] ~ dbern(ifelse(x[i,j]==1, 1, ilogit(
            alpha0 + alpha1*sex[i] + alpha2*(b_age[i]+j-1) + alpha3*u[i]
          )))
          
          v[i,j] ~ dbern(ilogit(
            beta0 + beta1*sex[i] + beta2*(b_age[i]+j-1) + beta3*u[i]
          ))
          
          z[i,j+1] ~ dbern(ifelse(z[i,j]==1, 1,
            ifelse(x[i,j]==0 || v[i,j]==0, 0, ilogit(
              eta0 + eta1*sex[i] + eta2*(b_age[i]+j-1) + eta3*u[i]
            ))
          ))
          
          y[i,j] ~ dbern(min(0.99999, ilogit(
              gamma0 + gamma1*sex[i] + gamma2*(b_age[i]+j-1) + gamma3*u[i]
            ) * exp(
              log(psi1)*x[i,j+1]*(1-z[i,j+1]) + log(psi2)*x[i,j+1]*z[i,j+1]
          )))
          
        }
        
      }
      
      gamma3 ~ dnorm(0, 1.0E-4)
      gamma2 ~ dnorm(0, 1.0E-4)
      gamma1 ~ dnorm(0, 1.0E-4)
      gamma0 ~ dnorm(0, 1.0E-4)
      eta3 ~ dnorm(0, 1.0E-4)
      eta2 ~ dnorm(0, 1.0E-4)
      eta1 ~ dnorm(0, 1.0E-4)
      eta0 ~ dnorm(0, 1.0E-4)
      beta3 ~ dnorm(0, 1.0E-4)
      beta2 ~ dnorm(0, 1.0E-4)
      beta1 ~ dnorm(0, 1.0E-4)
      beta0 ~ dnorm(0, 1.0E-4)
      alpha3 ~ dnorm(0, 1.0E-4)
      alpha2 ~ dnorm(0, 1.0E-4)
      alpha1 ~ dnorm(0, 1.0E-4)
      alpha0 ~ dnorm(0, 1.0E-4)
      psi2 ~ dunif(0, 10)
      psi1 ~ dunif(0, 10)
      sigma <- 1/sqrt(sigma_prec)
      sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
      
    }
  ")
  
  # !!!!! Set up data
  {
    I <- 20
    dat <- list(
      I = I,
      J = c(),
      sex = sample(c(0,1), size=I, replace=TRUE),
      b_age = sample(c(30:50), size=I, replace=TRUE),
      v = matrix(NA, nrow=I, ncol=15),
      x = cbind(rep(0,I), matrix(NA, nrow=I, ncol=15)),
      y = matrix(NA, nrow=I, ncol=15),
      z = cbind(rep(0,I), matrix(0, nrow=I, ncol=15))
    )
    alpha0 <- -5;  alpha1 <- 0.3;  alpha2 <- 0.1;  alpha3 <- 0.3
    gamma0 <- -6;  gamma1 <- 0.1;  gamma2 <- 0.05;  gamma3 <- 0.1
    psi1 <- 0.4
    psi2 <- 0.2
    for (i in 1:I) {
      x_last <- 0
      for (j in 1:15) {
        
        p_sero <- expit(
          alpha0 + alpha1*dat$sex[i] + alpha2*(dat$b_age[i]+j-1) + alpha3*dat$u[i]
        )
        dat$x[i,j] <- rbinom(1, 1, ifelse(x_last==1, 1, p_sero))
        x_last <- dat$x[i,j]
        
        p_outcome <- min(1, expit(
          gamma0 + gamma1*dat$sex[i] + gamma2*(dat$b_age[i]+j-1) + gamma3*dat$u[i]
        ) * exp(
          log(psi1)*dat$x[i,j]*(1-dat$z[i,j]) + log(psi2)*dat$x[i,j]*dat$z[i,j]
        ))
        dat$y[i,j] <- rbinom(1, 1, p_outcome)
        
        if (dat$y[i,j]==1 || j==15) {
          dat$J <- c(dat$J, j)
          break
        }
        
      }
    }
  }
  
  # Run model
  jm <- jags.model(
    file = textConnection(jags_code),
    data = dat,
    n.chains = mcmc$n.chains,
    n.adapt = mcmc$n.adapt
  )
  # mcmc <- list(n.chains=2, n.adapt=1000, n.burn=1000, n.iter=1000, thin=1)
  update(jm, n.iter = mcmc$n.burn)
  output <- coda.samples(
    model = jm,
    variable.names = c("alpha0", "alpha1", "alpha2", "alpha3"),
    n.iter = mcmc$n.iter,
    thin = mcmc$thin
  )
  # summary(output)
  
  # # MCMC diagnostics
  # if (FALSE) {
  #   var <- "beta_s_6"
  #   c1 <- output[[1]][1:mcmc$n.iter,var]
  #   c2 <- output[[2]][1:mcmc$n.iter,var]
  #   ggplot(
  #     data.frame(
  #       x = rep(c(1:mcmc$n.iter),2),
  #       y = c(c1,c2),
  #       chain = rep(c(1:2), each=mcmc$n.iter)
  #     ),
  #     aes(x=x, y=y, color=factor(chain))) +
  #     geom_line() +
  #     labs(title=var)
  # }
  
  # # Extract beta_s means
  # beta_s_hat <- c()
  # for (i in 1:6) {
  #   beta_s_hat[i] <- mean(
  #     unlist(lapply(output, function(l) {
  #       l[1:n_samp,paste0("beta_s_",i)]
  #     })),
  #     na.rm = TRUE
  #   )
  # }
  
  # # Construct covariance matrix of s terms
  # sigma_s_hat <- matrix(NA, nrow=6, ncol=6)
  # for (i in 1:6) {
  #   for (j in 1:6) {
  #     sigma_s_hat[i,j] <- cov(
  #       unlist(lapply(output, function(l) {l[1:n_samp,paste0("beta_s_",i)]})),
  #       unlist(lapply(output, function(l) {l[1:n_samp,paste0("beta_s_",j)]})),
  #       use = "complete.obs"
  #     )
  #   }
  # }
  
  return (999)

}


