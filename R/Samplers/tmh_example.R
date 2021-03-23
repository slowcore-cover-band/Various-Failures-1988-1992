###TMH Example
###09/13/20

rwmh = function(x0,llh,Hfac,scale,num_draws){
  draws = matrix(0,num_draws+1,length(x0)+1)
  llh0 = -llh(x0)
  draws[1,] = c(x0,llh0)
  kept = 1
  for(i in 1:num_draws){
    xnew = x0 + (scale*Hfac) %*% rnorm(length(x0))
    llhnew = -llh(xnew)
    if((llhnew-llh0)>log(runif(1))){
      x0 = xnew
      llh0 = llhnew
      kept = .99 * kept + 1
    }else{
      kept =  .99 * kept
    }
    if(i %% 100 == 1) print(paste("kept = ", kept,"and draws = ",i))
    draws[i+1,] = c(x0,llh0)
  }
return(draws)
}

tmh = function(x0,llh,scale,num_draws){
  draws = matrix(0,num_draws+1,length(x0)+1)
  llh0 = -llh(x0)
  draws[1,] = c(x0,llh0)
  kept = 1
  for(i in 1:num_draws){
    fac = rnorm(1,0,scale^2)
    fac = fac*sample(c(-1,1),length(x0),replace = TRUE,prob=c(0.5,0.5))
    xnew = x0 + fac
    llhnew = -llh(xnew)
    if((llhnew-llh0)>log(runif(1))){
      x0 = xnew
      llh0 = llhnew
      kept = .99 * kept + 1
    }else{
      kept =  .99 * kept
    }
    if(i %% 100 == 1) print(paste("kept = ", kept,"and draws = ",i))
    draws[i+1,] = c(x0,llh0)
  }
return(draws)
}

n_rep = 50
mu_par = 4
sig_par = 6.25

fake = rnorm(n_rep,mu_par,sig_par)

##mu is pars[1] and sigma^2 = pars[2]
llh = function(pars){
  alpha = 2
  beta = 5
  llh = -0.5*length(fake)*log(2*pi)-0.5*length(fake)*log(pars[2]^2)-0.5*sum((fake-pars[1])^2)/(pars[2]^2)
  llh = llh -0.5*log(2*pi) - 0.5*log(10) -0.5*(pars[1])^2/10
  llh = llh - alpha*log(beta) - lgamma(alpha) - (alpha+1)*log(pars[2]^2)-beta/pars[2]^2
  return(-llh)
}

start = optim(c(2,1),llh,hessian=TRUE)

xd = start$par
Hfac = t(chol(start$hessian))

tmh_test = tmh(xd,llh,1,5000)
rwmh_test = rwmh(xd,llh,Hfac,1,5000)
