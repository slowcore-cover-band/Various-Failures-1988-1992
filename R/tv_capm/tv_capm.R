##Time-Varying Parameters
##09/04/2020

##I consider four ways of dealing with time-varying parameters.
##This includes rolling window estimation, Exponential weighting (two different schemes),
##and estimating a model that explicitly accounts for parameter drift

ew_ols = function(y,x,discount = 0.95,robust = TRUE){
  beta_low = matrix(0,dim(x)[2],dim(x)[1]-dim(x)[2])
  beta = matrix(0,dim(x)[2],dim(x)[1]-dim(x)[2])
  beta_high = matrix(0,dim(x)[2],dim(x)[1]-dim(x)[2])
  for(i in (dim(x)[2]+1):dim(x)[1]){
    weight_mat = diag(c(discount^(0.5*(i:1))))
    wy = weight_mat %*% y[1:i]
    wx = weight_mat %*% x[1:i,]
    cur_beta = t(wx) %*% wx
    cur_beta = solve(cur_beta)
    cur_beta = cur_beta %*% (t(wx) %*% wy)
    beta[,(i-dim(x)[2])] = cur_beta
    if(robust == TRUE){
      resids = wy-wx %*% cur_beta
      resids = resids^2
      resids = diag(c(resids))
      ##print(dim(resids))
      ##print(dim(wx))
      std_err = solve((t(wx) %*% wx)) %*% (t(wx) %*% resids %*% wx) %*% solve((t(wx) %*% wx))
      std_err = matrix(sqrt(diag(std_err)),nrow = dim(cur_beta)[1])
    }else{
      resids = wy-wx%*%cur_beta
      var_est = sum(resids^2)/(length(resids)-dim(x)[2])
      std_err = var_est*solve((t(wx) %*% wx))
      std_err = matrix(sqrt(diag(std_err)),nrow = dim(cur_beta)[1])
    }
    beta_low[,(i-dim(x)[2])] = cur_beta - 1.95*std_err
    beta_high[,(i-dim(x)[2])] = cur_beta + 1.95*std_err
  }
  return(list(beta=beta,lower = beta_low, upper = beta_high))
}

ew_ols_rw = function(y,x,discount = 0.95,robust = TRUE){
  beta_low = matrix(0,dim(x)[2],dim(x)[1]-dim(x)[2])
  beta = matrix(0,dim(x)[2],dim(x)[1]-dim(x)[2])
  beta_high = matrix(0,dim(x)[2],dim(x)[1]-dim(x)[2])
  for(i in (dim(x)[2]+1):dim(x)[1]){
    weight_mat = diag(c(discount^(0.5*abs(((1-i):(dim(x)[1]-i))))))
    #print(diag(weight_mat))
    wy = weight_mat %*% y
    wx = weight_mat %*% x
    cur_beta = t(wx) %*% wx
    cur_beta = solve(cur_beta)
    cur_beta = cur_beta %*% (t(wx) %*% wy)
    beta[,(i-dim(x)[2])] = cur_beta
    if(robust == TRUE){
      resids = wy-wx %*% cur_beta
      resids = resids^2
      resids = diag(c(resids))
      ##print(dim(resids))
      ##print(dim(wx))
      std_err = solve((t(wx) %*% wx)) %*% (t(wx) %*% resids %*% wx) %*% solve((t(wx) %*% wx))
      std_err = matrix(sqrt(diag(std_err)),nrow = dim(cur_beta)[1])
    }else{
      resids = wy-wx%*%cur_beta
      var_est = sum(resids^2)/(length(resids)-dim(x)[2])
      std_err = var_est*solve((t(wx) %*% wx))
      std_err = matrix(sqrt(diag(std_err)),nrow = dim(cur_beta)[1])
    }
    beta_low[,(i-dim(x)[2])] = cur_beta - 1.95*std_err
    beta_high[,(i-dim(x)[2])] = cur_beta + 1.95*std_err
  }
  return(list(beta=beta,lower = beta_low, upper = beta_high))
}


roll_ols = function(y, x, window = 10, robust = TRUE){
  beta_low = matrix(0,dim(x)[2],dim(x)[1] - window + 1)
  beta = matrix(0,dim(x)[2],dim(x)[1] - window + 1)
  beta_high = matrix(0,dim(x)[2],dim(x)[1] - window + 1)
  for(i in 1:(dim(x)[1]-window+1)){
    wy = y[i:(i+window-1)]
    wx = x[i:(i+window-1),]
    cur_beta = t(wx) %*% wx
    cur_beta = solve(cur_beta)
    cur_beta = cur_beta %*% (t(wx) %*% wy)
    beta[,i] = cur_beta
    if(robust == TRUE){
      resids = wy-wx %*% cur_beta
      resids = resids^2
      resids = diag(c(resids))
      ##print(dim(resids))
      ##print(dim(wx))
      std_err = solve((t(wx) %*% wx)) %*% (t(wx) %*% resids %*% wx) %*% solve((t(wx) %*% wx))
      std_err = matrix(sqrt(diag(std_err)),nrow = dim(cur_beta)[1])
    }else{
      resids = wy-wx%*%cur_beta
      var_est = sum(resids^2)/(length(resids)-dim(x)[2])
      std_err = var_est*solve((t(wx) %*% wx))
      std_err = matrix(sqrt(diag(std_err)),nrow = dim(cur_beta)[1])
    }
    beta_low[,i] = cur_beta - 1.95*std_err
    beta_high[,i] = cur_beta + 1.95*std_err
  }
  return(list(beta=beta,lower = beta_low, upper = beta_high))
}

tvp_capm = function(y,x,num_draw){
  betas = array(0,c(dim(x)[2],dim(x)[1],num_draw))
  sigmas = matrix(0,dim(x)[2]+1,num_draw+1)
  sigmas[,1] = 0.01
  sigmas[1,1] = 1000
  #print(dim(x))
  if(dim(x)[2]==1){
    value = matrix(rep(0.1,dim(x)[2]),1,dim(x)[2])
    W = diag(value)
  }else{
    W = diag(rep(0.1,dim(x)[2]))
  }
  for(i in 1:num_draw){
    sigtt = array(0,c(dim(x)[2],dim(x)[2],length(y)+1))
    stt = matrix(0,dim(x)[2],length(y)+1)
    sigtt[,,1] = W
    if(dim(x)[2] == 1){
      stt[,1] = 1
    }else{
      stt[,1] = c(0,1)
    }
    ##Run Kalman Filter
    for(j in 1:length(y)){
      sigtt1 = sigtt[,,j] + W
      psi = matrix(x[j,], nrow = 1) %*% sigtt1 %*% matrix(x[j,], ncol = 1)
      psi = as.numeric(psi) + sigmas[1,i]
      psi = 1/psi
      adjust = psi * (sigtt1 %*% matrix(x[j,], ncol = 1)) * as.numeric(y[j]-matrix(x[j,],nrow=1) %*% stt[,j])
      stt[,j+1] = stt[,j] + adjust
      adjust = matrix(x[j,],ncol=1) %*% matrix(x[j,],nrow=1)
      sigtt[,,j+1] = sigtt1 - psi * sigtt1 %*% adjust %*% sigtt1
      #print("-----------------------------------------")
    }
    ##Compute Smoother Draws
    stt[,length(y)+1] = stt[,length(y)+1] + chol(sigtt[,,length(y)+1]) %*% rnorm(dim(x)[2])
    for(k in length(y):2){
      curmean = sigtt[,,k] %*% solve(sigtt[,,k]+W) %*% stt[,k+1]
      curmean = curmean + W %*% solve(sigtt[,,k]+W) %*% stt[,k]
      curvar = sigtt[,,k] %*% solve(sigtt[,,k]+W) %*% W
      stt[,k] = curmean + chol(curvar) %*% rnorm(dim(x)[2])
    }
    betas[,,i] = stt[,-1]
    ##print(dim(predict))
    ##print(dim(x))
    predict = betas[,,i] * t(x)
    predict = apply(predict,2,sum)
    ##Draw Sigmas
    epsilons = t(y) - predict
    she = 0.1 + 0.5*length(y)
    rae = (0.1 + 0.5*sum(epsilons^2))
    factor = 10
    sigmas[1,i+1] = factor/rgamma(1,she,rae)
    uvs = stt[,-1] - stt[,-dim(stt)[2]]
    uvs = uvs^2
    shuv = 0.1 + 0.5*length(y)
    if(dim(x)[2]==1){
        raev = (0.1+0.5*sum(uvs))
        sigmas[2:(dim(x)[2]+1),i+1] = 0.001/rgamma(1,shuv,raev)
    }else{
        raev = (0.1+0.5*apply(uvs,1,sum))
        sigmas[2:(dim(x)[2]+1),i+1] = 0.001/rgamma(rep(shuv,dim(x)[2]),raev)
    }
    if(i %% 100 == 1){
      print(paste(i,"draw(s) done!"))
    }
  }
  sigmas = sigmas[,-dim(sigmas)[2]]
  return(list(betas=betas,sigmas=sigmas))
}

reformat_date = function(data,delim_in,delim_out,mdy_order){
  new_data = rep(0,length(data))
  for(i in 1:length(data)){
    pos = gregexpr(delim_in,data[i])[[1]][1:2]
    if(pos[1]==2){
      if(pos[2]==5){
        m = paste0("0",substr(data[i],1,1))
        d = substr(data[i],3,4)
        y = substr(data[i],6,9)
      }else{
        m = paste0("0",substr(data[i],1,1))
        d = paste0("0",substr(data[i],3,3))
        y = substr(data[i],5,8)
      }
    }else{
      if(pos[2]==6){
        m = substr(data[i],1,2)
        d = substr(data[i],4,5)
        y = substr(data[i],7,10)
      }else{
        m = substr(data[i],1,2)
        d = paste0("0",substr(data[i],4,4))
        y = substr(data[i],6,9)
      }
    }
    cur_date = NULL
    cur_date[1] = m
    cur_date[2] = d
    cur_date[3] = y
    cur_date[mdy_order] = cur_date[1:3]
    new_data[i] = paste0(cur_date[1],delim_out,cur_date[2],delim_out,cur_date[3])
  }
  return(new_data)
}

###Simple numerical example

rm = rnorm(1000)
ri = 1+ 1.5*rm + rnorm(1000)
x = matrix(rm,ncol=1)#cbind(rep(1,1000),rm)
y = ri

exp_test = ew_ols(y,x,discount=0.99)
exp_rw_test = ew_ols_rw(y,x,discount=0.99)
roll_test = roll_ols(y,x,window=260)
tvp_test = tvp_capm(y,x,1000)
tvp_low = apply(tvp_test$betas,2,FUN = function(x) quantile(x,0.25))
tvp_high = apply(tvp_test$betas,2,FUN = function(x) quantile(x,0.75))
tvp_med = apply(tvp_test$betas,2,median)

plot(exp_test$beta[1,258:998],ylim = c(0,3),type="l")
lines(exp_test$lower[1,258:998],col="blue")
lines(exp_test$upper[1,258:998],col="blue")

lines(exp_rw_test$beta[1,258:998],col="forestgreen")
lines(exp_rw_test$lower[1,258:998],col="goldenrod")
lines(exp_rw_test$upper[1,258:998],col="goldenrod")

lines(roll_test$beta[1,],ylim = c(1,3),col="red")
lines(roll_test$lower[1,],col="green")
lines(roll_test$upper[1,],col="green")

lines(tvp_med[-c(1:259)],col="purple")
lines(tvp_low[-c(1:259)],col="orange")
lines(tvp_high[-c(1:259)],col="orange")

###Example for Starbucks

###Cleaning the data

sbux_data = read.csv("SBUX.csv")
lr20_data = read.csv("DGS20.csv")
sp5c_data = read.csv("^GSPC.csv")

sbux_data[,1] = reformat_date(sbux_data[,1],"/","-",c(2,3,1))

sbux_data = sbux_data[which(sbux_data[,1]==lr20_data[1,1]):dim(sbux_data)[1],c(1,5)]
sp5c_data = sp5c_data[which(sp5c_data[,1]==lr20_data[1,1]):dim(sp5c_data)[1],c(1,5)]

lr20_data = lr20_data[which(c(1:dim(lr20_data)[1])%%5==1),]
sp5c_data = sp5c_data[which(c(1:dim(sp5c_data)[1])%%5==1),]
sbux_data = sbux_data[which(c(1:dim(sbux_data)[1])%%5==1),]

min_len = min(c(dim(lr20_data)[1],dim(sp5c_data)[1],dim(sbux_data)[1]))

lr20_data = lr20_data[1:min_len,]
sp5c_data = sp5c_data[1:min_len,]
sbux_data = sbux_data[1:min_len,]

lr20 = as.numeric(as.character(lr20_data[,2]))
lr20 = lr20[-1]
sp5cr = 100*sp5c_data[-1,-1]/sp5c_data[-dim(sp5c_data)[1],-1]-100
sbuxr = 100*sbux_data[-1,-1]/sbux_data[-dim(sbux_data)[1],-1]-100

sp5cr = sp5cr[-which(is.na(lr20))]
sbuxr = sbuxr[-which(is.na(lr20))]

er5c = matrix(sp5cr - lr20[-which(is.na(lr20))],ncol=1)
ersb = matrix(sbuxr - lr20[-which(is.na(lr20))],ncol=1)

###No intercept
exp_test = ew_ols(ersb,er5c,discount=0.99)
exp_rw_test = ew_ols_rw(ersb,er5c,discount=0.99)
roll_test = roll_ols(ersb,er5c,window=260)
tvp_test = tvp_capm(ersb,er5c,1000)

tvp_low = apply(tvp_test$betas,2,FUN = function(x) quantile(x,0.25))
tvp_high = apply(tvp_test$betas,2,FUN = function(x) quantile(x,0.75))
tvp_med = apply(tvp_test$betas,2,median)

##Together
plot(exp_test$beta[1,258:(length(er5c)-1)],ylim = c(0,2),type="l")
lines(roll_test$beta[1,],ylim = c(1,3),col="red")
lines(exp_test$lower[1,258:(length(er5c)-1)],col="blue")
lines(exp_test$upper[1,258:(length(er5c)-1)],col="blue")
lines(roll_test$lower[1,],col="green")
lines(roll_test$upper[1,],col="green")
lines(tvp_med[-c(1:259)],col="purple")
lines(tvp_low[-c(1:259)],col="orange")
lines(tvp_high[-c(1:259)],col="orange")
lines(exp_rw_test$beta[1,258:998],col="forestgreen")
lines(exp_rw_test$lower[1,258:998],col="goldenrod")
lines(exp_rw_test$upper[1,258:998],col="goldenrod")

##Rolling
plot(roll_test$beta[1,],ylim = c(0,2),col="red",type="l")
lines(roll_test$lower[1,],col="green")
lines(roll_test$upper[1,],col="green")

##Exponential
plot(exp_test$beta[1,258:(length(er5c)-1)],ylim = c(0,2),type="l")
lines(exp_test$lower[1,258:(length(er5c)-1)],col="blue")
lines(exp_test$upper[1,258:(length(er5c)-1)],col="blue")

##Exponential (Kernel 2)
plot(exp_rw_test$beta[1,258:(length(er5c)-1)],ylim = c(0,2),type="l",col="forestgreen")
lines(exp_rw_test$lower[1,258:(length(er5c)-1)],col="goldenrod")
lines(exp_rw_test$upper[1,258:(length(er5c)-1)],col="goldenrod")

##TVP
plot(tvp_med[-c(1:259)],col="purple",type="l",ylim=c(0,2))
lines(tvp_low[-c(1:259)],col="orange")
lines(tvp_high[-c(1:259)],col="orange")

###With intercept
er5cint = cbind(rep(1,length(er5c)),er5c)

exp_test = ew_ols(ersb,er5cint,discount=0.99)
exp_rw_test = ew_ols_rw(ersb,er5cint,discount=0.99)
roll_test = roll_ols(ersb,er5cint,window=260)
tvp_test = tvp_capm(ersb,er5cint,1000)

atvp_low = apply(tvp_test$betas[1,,],1,FUN = function(x) quantile(x,0.25))
atvp_high = apply(tvp_test$betas[1,,],1,FUN = function(x) quantile(x,0.75))
atvp_med = apply(tvp_test$betas[1,,],1,median)

btvp_low = apply(tvp_test$betas[2,,],1,FUN = function(x) quantile(x,0.25))
btvp_high = apply(tvp_test$betas[2,,],1,FUN = function(x) quantile(x,0.75))
btvp_med = apply(tvp_test$betas[2,,],1,median)


###Alpha Together

plot(exp_test$beta[1,258:(length(er5c)-2)],ylim = c(-1,4),type="l")
lines(roll_test$beta[1,],ylim = c(1,3),col="red")
lines(exp_test$lower[1,258:(length(er5c)-2)],col="blue")
lines(exp_test$upper[1,258:(length(er5c)-2)],col="blue")
lines(roll_test$lower[1,],col="green")
lines(roll_test$upper[1,],col="green")
lines(atvp_med[-c(1:259)],col="purple")
lines(atvp_high[-c(1:259)],col="orange")
lines(atvp_low[-c(1:259)],col="orange")
lines(exp_rw_test$beta[1,258:(length(er5c)-2)],col="forestgreen")
lines(exp_rw_test$lower[1,258:(length(er5c)-2)],col="goldenrod")
lines(exp_rw_test$upper[1,258:(length(er5c)-2)],col="goldenrod")

##Alpha Rolling
plot(roll_test$beta[1,],ylim = c(-2,5),type="l",col="red")
lines(roll_test$lower[1,],col="green")
lines(roll_test$upper[1,],col="green")

##Alpha Exponential
plot(exp_test$beta[1,258:(length(er5c)-2)],ylim = c(-2,5),type="l")
lines(exp_test$lower[1,258:(length(er5c)-2)],col="blue")
lines(exp_test$upper[1,258:(length(er5c)-2)],col="blue")

##Alpha Exponential (Kernel 2)
plot(exp_rw_test$beta[1,258:(length(er5c)-2)],col="forestgreen",type="l",ylim=c(-2,5))
lines(exp_rw_test$lower[1,258:(length(er5c)-2)],col="goldenrod")
lines(exp_rw_test$upper[1,258:(length(er5c)-2)],col="goldenrod")

##Alpha TVP
plot(atvp_med[-c(1:259)],type="l",ylim=c(-2,5),col="purple")
lines(atvp_high[-c(1:259)],col="orange")
lines(atvp_low[-c(1:259)],col="orange")


###Beta Together
plot(exp_test$beta[2,258:(length(er5c)-2)],ylim = c(0,3),type="l")
lines(roll_test$beta[2,],ylim = c(1,3),col="red")
lines(exp_test$lower[2,258:(length(er5c)-2)],col="blue")
lines(exp_test$upper[2,258:(length(er5c)-2)],col="blue")
lines(roll_test$lower[2,],col="green")
lines(roll_test$upper[2,],col="green")
lines(btvp_med[-c(1:259)],col="purple")
lines(btvp_high[-c(1:259)],col="orange")
lines(btvp_low[-c(1:259)],col="orange")
lines(exp_rw_test$beta[2,258:(length(er5c)-2)],col="forestgreen")
lines(exp_rw_test$lower[2,258:(length(er5c)-2)],col="goldenrod")
lines(exp_rw_test$upper[2,258:(length(er5c)-2)],col="goldenrod")

##Beta Rolling
plot(roll_test$beta[2,],ylim = c(-2,5),type="l",col="red")
lines(roll_test$lower[2,],col="green")
lines(roll_test$upper[2,],col="green")

##Alpha Exponential
plot(exp_test$beta[2,258:(length(er5c)-2)],ylim = c(-2,5),type="l")
lines(exp_test$lower[2,258:(length(er5c)-2)],col="blue")
lines(exp_test$upper[2,258:(length(er5c)-2)],col="blue")

##Beta Exponential (Kernel 2)
plot(exp_rw_test$beta[2,258:(length(er5c)-2)],col="forestgreen",type="l",ylim=c(-2,5))
lines(exp_rw_test$lower[2,258:(length(er5c)-2)],col="goldenrod")
lines(exp_rw_test$upper[2,258:(length(er5c)-2)],col="goldenrod")

##Alpha TVP
plot(btvp_med[-c(1:259)],type="l",ylim=c(-2,5),col="purple")
lines(btvp_high[-c(1:259)],col="orange")
lines(btvp_low[-c(1:259)],col="orange")
