###Finite Difference HMC
###09/13/20

##Sims' numgrad function
numgrad <- function(fcn, x, ...) {
  ## fcn can return a vector, in which case numgrad returns a matrix.
  delta <- 1e-6
  ## delta <- 1e-8
  n <- length(x)
  ## we tolerate x's that may be n x 1, 1 x n, or R vectors (with no dim),
  ## but note that g comes out as n x k matrix regardless.
  tvec <- delta*diag(n)
  f0 <- fcn(x,...)
  k <- length(f0)
  g <- matrix(0,n,k)
  badg <- FALSE
  for (i in 1:n){
    scale <- 1
    tvecv <- tvec[,i]
    if(is.null(dim(x))){
      tvecv <- as.vector(tvecv)
    }else{
      dim(tvecv) <- dim(x)
    }
    g0 <- (fcn(x+scale*tvecv,...) - f0)/(scale*delta)
    if (max(abs(g0))< 1e15){
      g[i, ] <- as.vector(g0)
    }else{
      cat("bad gradient ------------------------\n")
      badg <- TRUE
    }
  }
  return(list(g=g,badg=badg))
}

leapfrog = function(llh,x0,p0,delT,L){
  low_del = delT/L
  for(i in 1:L){
    p0 = p0-0.5*low_del*numgrad(llh,x0)$g
    x0 = x0+low_del*p0
    p0 = p0-0.5*low_del*numgrad(llh,x0)$g
  }
  return(list(xnew=x0,pnew=p0))
}

hamiltonian = function(llh,x0,p0){
  ham = llh(x0) + 0.5*sum(p0^2)
  return(ham)
}

fd_hmc = function(x0,llh,delT,L,num_draws){
  draws = matrix(0,num_draws+1,length(x0)+1)
  llh0 = -llh(x0)
  draws[1,] = c(x0,llh0)
  kept = 1
  for(i in 1:num_draws){
    p0 = rnorm(length(x0))
    ham0 = hamiltonian(llh,x0,p0)
    new_state = leapfrog(llh,x0,p0,delT,L)
    #print(new_state$xnew)
    hamnew = hamiltonian(llh,new_state$xnew,new_state$pnew)
    if((exp(-hamnew)/exp(-ham0))>(runif(1))){
      x0 = new_state$xnew
      llh0 = -llh(x0)
      kept = 0.99*kept + 1
    }else{
      kept = 0.99*kept
    }
    if(i %% 100 == 1) print(paste("kept = ", kept,"and draws = ",i))
    draws[i+1,] = c(x0,llh0)
  }
  return(draws)
}
