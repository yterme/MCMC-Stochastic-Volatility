setwd("~/CoursENSAE/MonteCarlo/MCMC-Stochastic-Volatility/")
source("simulate_xy.R")

library(truncnorm)
T=100
i=1
sigma=0.1
phi=0.95
mu=0
x=simulate_xy(mu,sigma,phi)[[1]]

n=10000
burn=2500
sd_phi_kernel=1
sd_phi_prior=1

#simuler phi

simulate_phi <- function(x,mu,sigma){
  #Initialisation de la marche aléatoire
  phi_trace=rep(0,n-burn)
 
  # initialisation de phi
  phi=0

  for (i in 1:n){
    # noyau N(phi,sigma) tronquée sur [-1,1]
    phi_proposal=rtruncnorm(mean=phi,sd=sd_phi_kernel,a=-1,b=1,n=1)
    
    # calcul P(x| mu,sigma,phi)
    p_x=log_px(x,mu,sigma,phi)
    # Calcul de P(x|mu,sigma,phi_proposal)
    p_x_proposal=log_px(x,mu,sigma,phi_proposal)
    
    #on a suppose que mu suivait une prior N(0,sigma). Calcul du rapport de proba
    r=(exp(p_x_proposal)*dtruncnorm(mean = 0,sd = sd_phi_prior, a=-1, b=1, x = phi_proposal))/
      (exp(p_x)*dtruncnorm(mean = 0,sd = sd_phi_prior, a=-1, b=1, x = phi))

    alpha=runif(n = 1,min = 0,max = 1)
    
    if (r>alpha){ #on accepte la nouvelle valeur
      if (i>burn){
        phi_trace[i-burn]=phi_proposal}
      phi=phi_proposal
    }else{ #on rejette nouvelle valeur
      if (i>burn){
        phi_trace[i-burn]=phi
      }
    }
  }
  
  plot(phi_trace)
  hist(phi_trace)
  
  sample_phi=sample(phi_trace,1)
  return(sample_phi)
}



log_px= function(x,mu,sigma,phi){
  # calcul P(x| mu,phi,sigma)
  p_x=log(dnorm(mean = 0,sd = sigma, x = phi*(x[1]-mu))) #car x[0] vaut mu par definition donc mean=x[0]-mu=0
  for (j in 2:T){
    p_x=p_x+log(dnorm(mean = phi*(x[j-1]-mu),sd = sigma, x = x[j]-mu))
  }
  return(p_x)
}

simulate_phi(x,mu,sigma)
