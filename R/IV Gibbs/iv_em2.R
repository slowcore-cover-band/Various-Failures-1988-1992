###IV Gibbs sampling
###10/17/20


library(MCMCpack)
library(Matrix)

##In the function below y is the dependent variable, x is the endogenous variable, wrf is a matrix
##of controls in the reduced form, wi is a matrix of controls in the first stage,
##nu and omega are parameters of the inverse wishart prior on the
##error covariance matrix for the stacked system of the reduced form and the first stage, and
##and v0scale is the variance parameter in the prior on the coefficents, which is N(0,v0scale^{-1}I).
##The difficulty with this function is that the SUR-type draws are very costly for large
##datasets because you need to compute (X'Σ^{-1}X) for each draw of Σ. Alternatively,
##you could use the inverse of the Cholesky decomposition of Σ, but I don't think that would
##be any less expensive


iv_gibbs2 = function(y, x, z, wrf, wi, nu, omega, v0scale, num_draw = 1000, print_every = 100){
  nobs = length(y)
  gamma = matrix(1,dim(z)[2],num_draw)
  beta = matrix(0,1,num_draw)
  gammaa = matrix(0,dim(wrf)[2],num_draw)
  gamma2 = matrix(0,dim(wi)[2],num_draw)
  sigmas = array(0,c(2,2,num_draw+1))
  sigmas[,,1] = diag(2)
  for(i in 1:num_draw){
    ##Draw beta, Gammas
    what = cbind((z %*% gamma[,i]),wrf)
    xm = x - z %*% gamma[,i]
    xm = rbind(y,xm)
    xa = cbind(what, matrix(0,nobs,dim(wi)[2]))
    xa = rbind(xa, cbind(matrix(0,nobs,dim(what)[2]),wi))
    sighat = t(xa) %*% (solve(sigmas[,,i]) %x% Diagonal(nobs)) %*% xa + v0scale * diag(dim(xa)[2])
    sighat = solve(sighat)
    muhat = sighat %*% (t(xa) %*% (solve(sigmas[,,i]) %x% Diagonal(nobs)) %*% xm)
    curcoef = muhat + chol(sighat) %*% rnorm(dim(sighat)[2])
    beta[i] = curcoef[1]
    gammaa[,i] = curcoef[2:(dim(wrf)[2]+1)]
    gamma2[,i] = curcoef[-(1:(dim(wrf)[2]+1))]
    ##Draw gamma
    yrf = (y - wrf %*% gammaa[,i])/beta[i]
    yrf = rbind(yrf, (x - wi %*% gamma2[,i]))
    rssig = sigmas[,,i]
    rssig[1,1] = rssig[1,1]/beta[i]^2
    rssig[2,1] = rssig[2,1]/beta[i]
    rssig[1,2] = rssig[1,2]/beta[i]
    zz = rbind(z,z)
    sighat = t(zz) %*% (solve(rssig) %x% Diagonal(nobs)) %*% zz + v0scale * diag(dim(zz)[2])
    sighat = solve(sighat)
    muhat = sighat %*% (t(zz) %*% (solve(rssig) %x% Diagonal(nobs)) %*% yrf)
    curcoef = muhat + chol(sighat) %*% rnorm(dim(sighat)[2])
    gamma[,i] = as.numeric(curcoef)
    ##Draw sigmas
    rferr = y - beta[i] * (z %*%  gamma[,i]) - wrf %*% gammaa[,i]
    ierr  = x - z %*% gamma[,i] - wi %*% gamma2[,i]
    rferr = cbind(rferr,ierr)
    rferr = t(rferr) %*% rferr + omega
    sigmas[,,i+1] = riwish((nu+nobs),rferr)
    if(i %% print_every == 1){
      print(paste(i, "draw(s) done!"))
    }
  }
  sigmas = sigmas[,,-(num_draw+1)]
  return(list(beta = beta, gamma = gamma, gammaa = gammaa, gamma2 = gamma2, sigmas = sigmas))
}


load("akdata.RData")
qobl = levels(as.factor(akdataf$qob))
yobl = levels(as.factor(akdataf$yob))
pobl = levels(as.factor(akdataf$pob))

z_dum = NULL
for(i in 2:length(qobl)){
  z_dum = cbind(z_dum,as.numeric(akdataf$qob==qobl[i]))
}
w_dum = NULL
for(i in 2:length(pobl)){
  w_dum = cbind(w_dum,as.numeric(akdataf$pob==pobl[i]))
}
for(i in 2:length(yobl)){
  w_dum = cbind(w_dum,as.numeric(akdataf$yob==yobl[i]))
}
rand = sample(1:dim(akdataf)[1],10000)
w = cbind(rep(1,dim(w_dum)[1]),w_dum)
y = as.matrix(akdataf$logwage)
x = as.matrix(akdataf$educ)
w = w[rand,]
y = as.matrix(y[rand])
x = as.matrix(x[rand])
z = as.matrix(z_dum[rand,])
test_iv2 = iv_gibbs2(y,x,z,w,w,1,diag(2),0.001)
