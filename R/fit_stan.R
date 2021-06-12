#' Fit STAN model
#'
#' @param dat A list returned by transform_mcmc()
#' @param mcmc A list of the form list(n.adapt=1, n.burn=1, n.iter=1, thin=1,
#'     n.chains=1)
#' @return Entire STAN return output object
#' @notes Much of this code mirrors code in generate_data_events.R; ensure the
#'     two are in sync with one another

fit_stan <- function(dat, mcmc) {
  
  if (cfg$local) {
    options(mc.cores = parallel::detectCores()-1)
  }
  
  rstan_options(auto_write=TRUE)
  stan_code <- quote("
    
    data {
      int I;
      int max_T_i;
      int T_i[I];
      int v[I,max_T_i+1];
      int x[I,max_T_i+1];
      int y[I,max_T_i+1];
      int z[I,max_T_i+1];
      int sex[I];
      real b_age[I];
    }
    
    parameters {
      real u[I];
      real alpha0;
      real alpha1;
      real alpha2;
      real alpha3;
      real beta0;
      real beta1;
      real beta2;
      real beta3;
      real gamma0;
      real gamma1;
      real gamma2;
      real gamma3;
      real eta0;
      real eta1;
      real eta2;
      real eta3;
      real<lower=0> psi1;
      real<lower=0> psi2;
      real<lower=0> sigma;
    }
    
    model {
      
      for (i in 1:I) {
        
        u ~ normal(0, sigma);
        
        for (j in 2:(T_i[i]+1)) {
          
          x[i,j] ~ bernoulli(x[i,j-1]==1 ? 0.99999 : inv_logit(
            alpha0 + alpha1*sex[i] + alpha2*(b_age[i]+inv(12)*(j-1)) + alpha3*u[i]
          ));
          
          v[i,j] ~ bernoulli(inv_logit(
            beta0 + beta1*sex[i] + beta2*(b_age[i]+inv(12)*(j-1)) + beta3*u[i]
          ));
          
          z[i,j] ~ bernoulli(z[i,j-1]==1 ? 0.99999 :
            (x[i,j]==0 || v[i,j]==0 ? 0.00001 : inv_logit(
              eta0 + eta1*sex[i] + eta2*(b_age[i]+inv(12)*(j-1)) + eta3*u[i]
            ))
          );
          
          y[i,j] ~ bernoulli(fmin(0.99999, inv_logit(
              gamma0 + gamma1*sex[i] + gamma2*(b_age[i]+inv(12)*(j-1)) + gamma3*u[i]
            ) * exp(
              log(psi1)*x[i,j]*(1-z[i,j]) + log(psi2)*x[i,j]*z[i,j]
          )));
          
        }
        
        s = sum(v[i]>0) ? 1 : 0
        test_first = s * min( (1:T_i[i]) + T_i[i]*(1-v[i]) )
        test_last = s * max( (1:T_i[i])*v[i] )
        case = s ? 1 : x[i,test_last+1]==0 ? 2 : x[i,test_first+1]==1 ? 4 : 3
        
        !!!!! compute delta and deltax

      }
  
  }")
  
  fit <- stan(
    model_code = stan_code,
    data = dat,
    chains = mcmc$n.chains,
    iter = mcmc$n.burn + mcmc$n.iter,
    warmup = mcmc$n.burn,
    thin = mcmc$thin
  )
  
  if (F) {
    var <- "psi1"
    c1 <- extract(fit)[[var]][1:mcmc$n.iter]
    c2 <- extract(fit)[[var]][(1:mcmc$n.iter)+1000]
    c3 <- extract(fit)[[var]][(1:mcmc$n.iter)+2000]
    ggplot(
      data.frame(
        x = rep(c(1:mcmc$n.iter),3),
        y = c(c1,c2,c3),
        chain = rep(c(1:3), each=mcmc$n.iter)
      ),
      aes(x=x, y=y, color=factor(chain))
    ) +
      geom_line() +
      labs(title=var)
    # stan_diag(fit)
    # stan_mcse(fit)
  }

  # # !!!!! Priors
  # psi1 ~ uniform(0, 10);
  # psi2 ~ uniform(0, 10);
  # sigma <- 1/sqrt(sigma_prec) # JAGS format
  # sigma_prec ~ dgamma(1.0E-3, 1.0E-3) # JAGS format
  
  
  
  
  
  
  
  
  # !!!!!!!!!!!!!!!!!!!!!!!!!
  
  # Run model
  return (output)
  
  # !!!!! Sort the code below
  
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
  
}
