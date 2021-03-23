###BMA
###09/14/20

#BMA MDD
bma_mdd_em = function(g,XX,Xy,yy,n_obs){
  if(length(XX)==1){
    kj=1
    uu = Xy * (1/XX) * Xy
  }else{
    kj = dim(XX)[1]
    uu = t(Xy) %*% solve(XX) %*% Xy
  }
  uu = yy - as.numeric(uu)
  mdd = (-0.5*n_obs+0.5)*log(uu/(1+g) + g*yy/(1+g))+log(g/(1+g))*(0.5*kj)
  return(mdd)
}

#BMA Gibbs
bma_mcmc_em = function(y, x, g, n_its){
  n_obs = length(y)
  yy = sum(y^2)
  Xy = t(x) %*% y
  XX = t(x) %*% x
  beta = matrix(0,dim(x)[2],n_its)
  sigma = matrix(0.1,1,n_its)
  models = matrix(0,dim(x)[2],n_its+1)
  mdds = matrix(0,1,n_its+1)
  ##randomly draw a starting model
  curmodel = sample(1:dim(x)[2],sample(1:dim(x)[2],1))
  models[curmodel,1] = 1
  mdd0 = bma_mdd_em(g, XX[curmodel,curmodel],Xy[curmodel],yy,n_obs)
  #print(mdd0)
  mdds[1] = mdd0
  kept = 1
  print(curmodel)
  for(i in 1:n_its){
    #Draw new model
    newj = sample(1:dim(x)[2],1)
    if(length(which(curmodel==newj))==0){
      newmodel = c(curmodel,newj)
    }else{
      newmodel = curmodel[-which(curmodel==newj)]
      if(length(newmodel)==0){
        newmodel = sample(which(1:dim(x)[2]%in%curmodel),1)
      }
    }
    #print(newmodel)
    mddnew = bma_mdd_em(g, XX[newmodel,newmodel],Xy[newmodel],yy,n_obs)
    #print(mddnew)
    if((mddnew-mdd0)>log(runif(1))){
      curmodel = newmodel
      mdd0 = mddnew
      kept = 0.99*kept + 1
    }else{
      kept = 0.99*kept
    }
    models[curmodel,i+1] = 1
    mdds[i+1] = mdd0
    #Draw parameters of model
    shs = n_obs/2-0.5
    if(length(curmodel) ==1){
      ras =(g/(1+g))*Xy[curmodel] * (1/(XX[curmodel,curmodel])) * Xy[curmodel]
      ras = 0.5*yy - 0.5*(g/(1+g))*ras
      sigma[i] = 1/rgamma(1,shs,ras)
      varbet = sigma[i]^2 * (g/(1+g)) * (1/XX[curmodel,curmodel])
      beta[curmodel,i] = (g/(1+g)) * (1/XX[curmodel,curmodel]) * Xy[curmodel] + sqrt(varbet) * rnorm(length(curmodel))
    }else{
      ras = (g/(1+g))*t(Xy[curmodel])%*% solve(XX[curmodel,curmodel]) %*% Xy[curmodel]
      ras = 0.5*yy - 0.5*(g/(1+g))*ras
      sigma[i] = 1/rgamma(1,shs,ras)
      varbet = sigma[i]^2 * (g/(1+g)) * solve(XX[curmodel,curmodel])
      beta[curmodel,i] = (g/(1+g)) * solve(XX[curmodel,curmodel]) %*% Xy[curmodel] + chol(varbet) %*% rnorm(length(curmodel))
    }
    if(i%%1000 == 1){
      print(paste(i, "draw(s) and", kept))
    }
  }
  return(list(models = models, beta = beta, sigma = sigma, mdds = mdds))
}

x = matrix(rnorm(500),ncol=5)
coefs = matrix(c(1,0,2,0,-1),ncol=1)
y = x %*% coefs + rnorm(100)
y = y - mean(y)
