#' Fit JAGS model
#'
#' @param dat A list returned by transform_jags()
#' @param mcmc A list of the form list(n.adapt=1, n.burn=1, n.iter=1, thin=1,
#'     n.chains=1)
#' @return Entire JAGS return output object

fit_jags <- function(dat, mcmc) {
  
  # Declare JAGS model (with U latent variable as random effect)
  jags_code_u <- quote("
    model {
      
      for (i in 1:I) {
        
        u[i] ~ dnorm(0, 1/(sigma^2))
        
        for (j in 2:(J[i]+1)) {
          
          x[i,j] ~ dbern(ifelse(x[i,j-1]==1, 0.99999, ilogit(
            alpha0 + alpha1*sex[i] + alpha2*(b_age[i]+(j-1)/12) + alpha3*u[i]
          )))
          
          v[i,j] ~ dbern(ilogit(
            beta0 + beta1*sex[i] + beta2*(b_age[i]+(j-1)/12) + beta3*u[i]
          ))
          
          z[i,j] ~ dbern(ifelse(z[i,j-1]==1, 0.99999,
            ifelse(x[i,j]==0 || v[i,j]==0, 0.00001, ilogit(
              eta0 + eta1*sex[i] + eta2*(b_age[i]+(j-1)/12) + eta3*u[i]
            ))
          ))
          
          y[i,j] ~ dbern(min(0.99999, ilogit(
              gamma0 + gamma1*sex[i] + gamma2*(b_age[i]+(j-1)/12) + gamma3*u[i]
            ) * exp(
              log(psi1)*x[i,j]*(1-z[i,j]) + log(psi2)*x[i,j]*z[i,j]
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
  
  # Run model
  jm <- jags.model(
    file = textConnection(jags_code_u),
    data = dat,
    n.chains = mcmc$n.chains,
    n.adapt = mcmc$n.adapt
  )
  (getS3method("update","jags"))(jm, n.iter = mcmc$n.burn)
  output <- coda.samples(
    model = jm,
    variable.names = c("alpha0", "alpha1", "alpha2", "alpha3",
                       "gamma0", "gamma1", "gamma2", "gamma3",
                       "psi1", "psi2"),
    n.iter = mcmc$n.iter,
    thin = mcmc$thin
  )
  
  return (output)
  
  # !!!!! Sort the code below
  
  # summary(output)
  # (getS3method("summary","mcmc.list"))(output)
  
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
  
}
