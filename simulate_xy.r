library(truncnorm)
T=100
y=rnorm(T)
i=1
sigma=0.1
phi=0.95
mu=0


simulate_xy <- function(mu,sigma,phi){
  x_list=rep(0,T)
  y_list=rep(0,T)
  x=mu
  for (i in 1:T){
    x=mu+rnorm(mean = phi*(x-mu), sd=sigma,n = 1)
    x_list[i]=x
    y_list[i]=rnorm(0,exp(x),n=1)
  }
  
  return(list(x_list,y_list))
}

l=simulate_xy(mu,sigma,phi)
x=l[[1]]