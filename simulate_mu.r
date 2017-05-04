setwd("~/CoursENSAE/MonteCarlo/MCMC-Stochastic-Volatility/")

source("simulate_xy.R")

T=100
sigma=0.1
mu=0
phi=0.95
x=simulate_xy(x,mu,sigma)[[1]]

n=10000
burn=2500
sd_mu_kernel=sigma
sd_mu_prior=sigma

simulate_mu <- function(x,phi,sigma){

  mu_trace=rep(0,n-burn)
  mu=0

  for (i in 1:n){
    # noyau gaussien N(0,1)
    mu_proposal=rnorm(mean=mu,sd=sd_mu_kernel,n=1)
  
    # calcul P(x| mu,sigma,phi)
    p_x=log_px(x,mu,sigma,phi)
    # Calcul de P(x|mu,sigma,phi)
    p_x_proposal=log_px(x,mu_proposal,sigma,phi)
      
    #on a suppose que mu suivait une prior N(0,sigma)
    r=(exp(p_x_proposal)*dnorm(mean = 0,sd =sd_mu_prior, x = mu_proposal))/
      (exp(p_x)*dnorm(mean = 0,sd = sd_mu_prior, x = mu))
  
    alpha=runif(n = 1,min = 0,max = 1)
  
    if (r>alpha){
      if (i>burn){
        mu_trace[i-burn]=mu_proposal}
      mu=mu_proposal
    }else{
      if (i>burn){
      mu_trace[i-burn]=mu
      }
    }
  }
  
  plot(mu_trace)
  hist(mu_trace)
  
  sample_mu=sample(mu_trace,1)
  return(sample_mu)
}


simulate_mu(x,phi,sigma)

