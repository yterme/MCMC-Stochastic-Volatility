setwd("~/CoursENSAE/MonteCarlo/")

source("simulate_xy.R")

library(statmod)
T=100
sigma=0.1
mu=0
phi=0.95
x=simulate_xy(x,mu,sigma)[[1]]

n=10000
burn=2500
mu_trace=rep(0,n-burn)

mean_sigma_kernel=0.1
mean_sigma_prior= 0.1
sd_sigma_kernel=sigma
shape_sigma_prior=sigma



simulate_sigma <- function(x,mu,phi){
  sigma_trace=rep(0,n-burn)
  sigma=0.1
  for (i in 1:n){
    # noyau gaussien N(0,1)
    sigma_proposal=rnorm(mean=sigma,sd=sd_sigma_kernel,n=1)
    
    # calcul P(x| mu,sigma,phi)
    p_x=log_px(x,mu,sigma,phi)
    # Calcul de P(x|mu,sigma,phi)
    p_x_proposal=log_px(x,mu,sigma_proposal,phi)
    
    #on a suppose que mu suivait une prior N(0,sigma)
    r=(exp(p_x_proposal)*dinvgauss(mean = mean_sigma_prior,shape = shape_sigma_prior, x = mu_proposal))/
      (exp(p_x)*dinvgauss(mean = mean_sigma_prior, shape = shape_sigma_prior, x = mu))
    
    alpha=runif(n = 1,min = 0,max = 1)
    
    if (r>alpha){
      if (i>burn){
        sigma_trace[i-burn]=sigma_proposal}
      sigma=sigma_proposal
    }else{
      if (i>burn){
        sigma_trace[i-burn]=sigma
      }
    }
  }
  
  plot(sigma_trace)
  hist(sigma_trace)
  
  sample_sigma=sample(sigma_trace,1)
  return(sample_sigma)
}


simulate_sigma(x,mu,phi)

