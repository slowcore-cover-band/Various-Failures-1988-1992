###IV Gibbs sampling
###10/17/20

#This is basically a SUR estimation procedure, but because we have exact
#identifaction we can recover the coefficient on the endogenous variable of
#interest

#This is acceptable if you have an exactly-identified model, but it wouldn't
#work if you have an overidentified model (# instruments > # endogenous vars)

library(MCMCpack)

iv_gibbs = function(y, x, z, wrf, wi, nu, omega, b0, v0inv, num_draw = 1000, print_every = 100){
nobs = length(y)
yx = rbind(y,x)
Xrf = cbind(z,wrf)
Xi = cbind(z,wi)
Xa = cbind(Xrf,matrix(0,nobs,dim(Xi)[2]))
Xa = rbind(Xa,cbind(matrix(0,nobs,dim(Xrf)[2]),Xi))
coefs = matrix(0,dim(Xa)[2],num_draw)
sigmas = array(0,c(2,2,num_draw+1))
sigmas[,,1] = diag(2)
  for(i in 1:num_draw){
    ##Draw coefficients
    sigstack = solve(sigmas[,,i]) %x% diag(nobs)
    sighat = t(Xa) %*% sigstack %*% Xa + v0inv
    sighat = solve(sighat)
    muhat = sighat %*% (t(Xa) %*% sigstack %*% yx + v0inv %*% b0)
    coefs[,i] = muhat + chol(sighat) %*% rnorm(dim(sighat)[2])
    ##Draw covariance matrix
    errs = yx - Xa %*% coefs[,i]
    errs = cbind(errs[1:nobs],errs[-c(1:nobs)])
    S = t(errs) %*% errs
    S = S + omega
    sigmas[,,i+1] = riwish((nu+nobs),S)
    if(i %% print_every == 1){
      print(paste(i, "draw(s) done!"))
    }
  }
  sigmas = sigmas[,,-(num_draw+1)]
  rfcoef = coefs[1:dim(Xrf)[2],]
  icoef  = coefs[-c(1:dim(Xrf)[2]),]
  return(list(rfcoef = rfcoef, icoef = icoef, sigmas = sigmas))
}
