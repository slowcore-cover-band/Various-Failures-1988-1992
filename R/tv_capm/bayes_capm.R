##Bayesian Rolling/Weighted Regressions
##09/06/20

#it's possible to use the heteroskedasticity robust model here
#but the heteroskedasticity that the model will pick up will mostly
#be due to the weights we applied. As a result, if you use that function
#use it only for rolling window regressions


bayes_ew = function(y, x, rolling = FALSE, window = 260, discount=0.95, all_obs = TRUE, num_draw = 1000){
  betas = array(0,c(dim(x)[2],dim(x)[1]-dim(x)[2],num_draw))
  sigmas = matrix(0.1,nrow=num_draw+1,ncol=dim(x)[1]-dim(x)[2])
  #print(dim(sigmas))
  if(dim(x)[2]==1){
    b = matrix(1,1,1)
  }else{
    b = matrix(1,nrow = dim(x)[2],ncol=1)
    b[1] = 0
  }
  if(rolling == TRUE){
    begin = 1
    end = dim(x)[1]-window+1
  }else{
    begin = dim(x)[2]+1
    end = dim(x)[1]
  }
  for(i in begin:end){
    if(rolling == TRUE){
      wy = y[i:(i+window-1)]
      wx = x[i:(i+window-1),]
    }else{
      if(all_obs == TRUE){
        weight_mat = diag(c(discount^(0.5*abs(((1-i):(dim(x)[1]-i))))))
        wx = weight_mat %*% x
        wy = weight_mat %*% y
      }else{
        weight_mat = diag(c(discount^(0.5*(i:1))))
        wy = weight_mat %*% y[1:i]
        wx = weight_mat %*% x[1:i,]
      }
    }
    for(j in 1:num_draw){

      if(dim(x)[2]==1){
        sighat = 1+(1/sigmas[j,i-begin+1])*sum(wx^2)
        #print(sigmas[j,i-begin+1+1])
        sighat = 1/sighat
        sighat = matrix(sighat,1,1)
      }else{
        sighat = diag(dim(x)[2])+(1/sigmas[j,i-begin+1])*(t(wx)%*%wx)
        sighat = solve(sighat)
      }

      muhat = sighat %*% b + (1/sigmas[j,i-begin+1]) * (sighat %*% t(wx) %*% wy)
      if(dim(x)[2]==1){
        sighat = sqrt(sighat)
      }else{
        sighat = chol(sighat)
      }
      betas[,i-begin+1,j] = muhat + sighat %*% rnorm(dim(x)[2])
      resids = wy-wx %*% matrix(betas[,i-begin+1,j],ncol=1)
      #print(sum(resids^2))
      shas = 0.1 + 0.5*length(wy)
      ras  = 0.1 + 0.5*sum(resids^2)
      sigmas[j+1,i-begin+1] = 1/rgamma(1,shas,ras)
    }
    if(i%%100 == 1){
      print(paste(i-begin+1,"out of", end-begin+1, "periods completed."))
    }
  }
  sigmas = sigmas[-(num_draw+1),]
  if(rolling==TRUE){
    betas = betas[,1:(dim(x)[1]-window+1),]
    sigmas = sigmas[,1:(dim(x)[1]-window+1)]
  }
  return(list(betas=betas,sigmas=sigmas))
}

bayes_ewr = function(y, x, rolling = FALSE, window = 260, discount=0.95, all_obs = TRUE, num_draw = 1000){
  betas = array(0,c(dim(x)[2],dim(x)[1]-dim(x)[2],num_draw))
  sigmas = matrix(0.1,num_draw+1,dim(x)[1]-dim(x)[2])

  if(dim(x)[2]==1){
    b = matrix(1,1,1)
  }else{
    b = matrix(1,nrow = dim(x)[2],ncol=1)
    b[1] = 0
  }
  if(rolling == TRUE){
    begin = 1
    end = dim(x)[1]-window+1
    psis = array(1,c(end,dim(x)[1]-dim(x)[2],num_draw+1))
  }else{
    begin = dim(x)[2]+1
    end = dim(x)[1]
    psis = array(1,c(dim(x)[1],dim(x)[1]-dim(x)[2],num_draw+1))
  }
  for(i in begin:end){
    if(rolling == TRUE){
      wy = y[i:(i+window-1)]
      wx = x[i:(i+window-1),]
    }else{
      if(all_obs == TRUE){
        weight_mat = diag(c(discount^(0.5*abs(((1-i):(dim(x)[1]-i))))))
        wx = weight_mat %*% x
        wy = weight_mat %*% y
      }else{
        weight_mat = diag(c(discount^(0.5*(i:1))))
        wy = weight_mat %*% y[1:i]
        wx = weight_mat %*% x[1:i,]
      }
    }

    for(j in 1:num_draw){
      if(dim(x)[2]==1){
        if(all_obs==FALSE){
          inside = sum(wx^2*psis[1:i,i-begin+1,j])
        }else{
          inside = sum(wx^2*psis[,i-begin+1,j])
        }

        inside = 1/inside
        sighat = 1+(1/sigmas[j,i-begin+1])*sum(wx^2)^2*inside
        sighat = 1/sighat
      }else{
        if(all_obs ==FALSE){
          inside = t(wx) %*% diag(psis[1:i,i-begin+1,j]) %*% wx
        }else{
          inside = t(wx) %*% diag(psis[,i-begin+1,j]) %*% wx
        }
        inside = solve(inside)
        sighat = diag(dim(x)[2])+(1/sigmas[j,i-begin+1])*(t(wx) %*% wx) %*% inside %*% (t(wx) %*% wx)
        sighat = solve(sighat)
      }
      muhat = sighat %*% b + (1/sigmas[j,i-begin+1]) * (sighat %*% t(wx) %*% wx) %*% inside %*% (t(wx) %*% wy)
      if(dim(x)[2]==1){
        sighat = sqrt(sighat)
      }else{
        sighat = chol(sighat)
      }
      betas[,i-begin+1,j] = muhat + sighat %*% rnorm(dim(x)[2])
      resids = wy-wx %*% matrix(betas[,i-begin+1,j],ncol=1)
      if(all_obs==FALSE){
        st_resids = resids/sqrt(psis[1:i,i-begin+1,j])
      }else{
        st_resids = resids/sqrt(psis[,i-begin+1,j])
      }
      shas = 0.1 + 0.5*length(wy)
      ras  = 0.1 + 0.5*sum(st_resids^2)
      sigmas[j+1,i-begin+1] = 1/rgamma(1,shas,ras)
      for(k in 1:dim(wx)[1]){
        shap = 3.5
        rap  = 2 + resids[k]^2/sigmas[j+1,i-begin+1]
        psis[k,i-begin+1,j+1] = 1/rgamma(1,shap,rap)
      }
    }
    if(i%%100 == 1){
      print(paste(i-begin+1,"out of", end-begin+1, "periods completed."))
    }
  }
  sigmas = sigmas[-(num_draw+1),]
  psis = psis[,,-(num_draw+1)]
  if(rolling==TRUE){
    betas = betas[,1:(dim(x)[1]-window+1),]
    psis = psis[,1:(dim(x)[1]-window+1),]
    sigmas = sigmas[,1:(dim(x)[1]-window+1)]
  }
  return(list(betas = betas, sigmas = sigmas,psis = psis))
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

##Exponential (Homoskedastic)
test_bew = bayes_ew(y,x,discount = 0.99,all_obs=FALSE,num_draw=100)
bew_low = apply(test_bew$betas,2,FUN = function(x) quantile(x,0.025))
bew_high = apply(test_bew$betas,2,FUN = function(x) quantile(x,0.975))
bew_med = apply(test_bew$betas,2,median)

##Exponential (Kernel 2) (Homoskedastic)
test_bew2 = bayes_ew(y,x,discount = 0.99,num_draw=100)
bew2_low = apply(test_bew2$betas,2,FUN = function(x) quantile(x,0.025))
bew2_high = apply(test_bew2$betas,2,FUN = function(x) quantile(x,0.975))
bew2_med = apply(test_bew2$betas,2,median)

##Rolling (Homoskedastic)
test_broll = bayes_ew(y,x,rolling=TRUE,num_draw=100)
broll_low = apply(test_broll$betas,1,FUN = function(x) quantile(x,0.025))
broll_high = apply(test_broll$betas,1,FUN = function(x) quantile(x,0.975))
broll_med = apply(test_broll$betas,1,median)

##Together
plot(bew_med[258:999],type="l",ylim = c(0,2))
lines(bew_high[258:999],col="blue")
lines(bew_low[258:999],col="blue")

lines(bew2_med[258:999],col="red")
lines(bew2_high[258:999],col="green")
lines(bew2_low[258:999],col="green")

lines(broll_med,col="purple")
lines(broll_high,col="orange")
lines(broll_low,col="orange")

##Apart
plot(bew_med[258:999],type="l",ylim = c(0,2))
lines(bew_high[258:999],col="blue")
lines(bew_low[258:999],col="blue")

plot(bew2_med[258:999],col="red",ylim=c(0,2),type="l")
lines(bew2_high[258:999],col="green")
lines(bew2_low[258:999],col="green")

plot(broll_med,col="purple",ylim=c(0,2),type="l")
lines(broll_high,col="orange")
lines(broll_low,col="orange")

###Starbucks Example
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

er5cint = cbind(rep(1,length(er5c)),er5c)

##No intercept
##Exponential (Homoskedastic)
test_bew = bayes_ew(ersb,er5c,discount = 0.99,all_obs=FALSE)
bew_low = apply(test_bew$betas,2,FUN = function(x) quantile(x,0.025))
bew_high = apply(test_bew$betas,2,FUN = function(x) quantile(x,0.975))
bew_med = apply(test_bew$betas,2,median)

##Exponential (Kernel 2) (Homoskedastic)
test_bew2 = bayes_ew(ersb,er5c,discount = 0.99)
bew2_low = apply(test_bew2$betas,2,FUN = function(x) quantile(x,0.025))
bew2_high = apply(test_bew2$betas,2,FUN = function(x) quantile(x,0.975))
bew2_med = apply(test_bew2$betas,2,median)

##Rolling (Homoskedastic)
test_broll = bayes_ew(ersb,er5c,rolling=TRUE)
broll_low = apply(test_broll$betas,1,FUN = function(x) quantile(x,0.025))
broll_high = apply(test_broll$betas,1,FUN = function(x) quantile(x,0.975))
broll_med = apply(test_broll$betas,1,median)


##Together
plot(bew_med[258:1273],type="l",ylim = c(0,2))
lines(bew_high[258:1273],col="blue")
lines(bew_low[258:1273],col="blue")

lines(bew2_med[258:1273],col="red")
lines(bew2_high[258:1273],col="green")
lines(bew2_low[258:1273],col="green")

lines(broll_med,col="purple")
lines(broll_high,col="orange")
lines(broll_low,col="orange")

##Apart
plot(bew_med[258:1273],type="l",ylim = c(0,2))
lines(bew_high[258:1273],col="blue")
lines(bew_low[258:1273],col="blue")

plot(bew2_med[258:1273],col="red",ylim=c(0,2),type="l")
lines(bew2_high[258:1273],col="green")
lines(bew2_low[258:1273],col="green")

plot(broll_med,col="purple",ylim=c(0,2),type="l")
lines(broll_high,col="orange")
lines(broll_low,col="orange")


##With intercept
##Exponential (Homoskedastic)
test_bew = bayes_ew(ersb,er5cint,discount = 0.99,all_obs=FALSE)
abew_low = apply(test_bew$betas[1,,],1,FUN = function(x) quantile(x,0.025))
abew_high = apply(test_bew$betas[1,,],1,FUN = function(x) quantile(x,0.975))
abew_med = apply(test_bew$betas[1,,],1,median)

bbew_low = apply(test_bew$betas[2,,],1,FUN = function(x) quantile(x,0.025))
bbew_high = apply(test_bew$betas[2,,],1,FUN = function(x) quantile(x,0.975))
bbew_med = apply(test_bew$betas[2,,],1,median)

##Exponential (Kernel 2) (Homoskedastic)
test_bew2 = bayes_ew(ersb,er5cint,discount = 0.99)
abew2_low = apply(test_bew2$betas[1,,],1,FUN = function(x) quantile(x,0.025))
abew2_high = apply(test_bew2$betas[1,,],1,FUN = function(x) quantile(x,0.975))
abew2_med = apply(test_bew2$betas[1,,],1,median)
bbew2_low = apply(test_bew2$betas[2,,],1,FUN = function(x) quantile(x,0.025))
bbew2_high = apply(test_bew2$betas[2,,],1,FUN = function(x) quantile(x,0.975))
bbew2_med = apply(test_bew2$betas[2,,],1,median)


##Rolling (Homoskedastic)
test_broll = bayes_ew(ersb,er5cint,rolling=TRUE)
abroll_low = apply(test_broll$betas[1,,],1,FUN = function(x) quantile(x,0.025))
abroll_high = apply(test_broll$betas[1,,],1,FUN = function(x) quantile(x,0.975))
abroll_med = apply(test_broll$betas[1,,],1,median)

bbroll_low = apply(test_broll$betas[2,,],1,FUN = function(x) quantile(x,0.025))
bbroll_high = apply(test_broll$betas[2,,],1,FUN = function(x) quantile(x,0.975))
bbroll_med = apply(test_broll$betas[2,,],1,median)



##Alpha Together
plot(abew_med[258:1273],type="l",ylim = c(-2,5))
lines(abew_high[258:1273],col="blue")
lines(abew_low[258:1273],col="blue")

lines(abew2_med[258:1273],col="red")
lines(abew2_high[258:1273],col="green")
lines(abew2_low[258:1273],col="green")

lines(abroll_med,col="purple")
lines(abroll_high,col="orange")
lines(abroll_low,col="orange")

##Alpha Apart
plot(abew_med[258:1273],type="l",ylim = c(-2,5))
lines(abew_high[258:1273],col="blue")
lines(abew_low[258:1273],col="blue")

plot(abew2_med[258:1273],col="red",ylim=c(-2,5),type="l")
lines(abew2_high[258:1273],col="green")
lines(abew2_low[258:1273],col="green")

plot(abroll_med,col="purple",ylim=c(-2,5),type="l")
lines(abroll_high,col="orange")
lines(abroll_low,col="orange")

##Beta Together
plot(bbew_med[258:1273],type="l",ylim = c(0,2))
lines(bbew_high[258:1273],col="blue")
lines(bbew_low[258:1273],col="blue")

lines(bbew2_med[258:1273],col="red")
lines(bbew2_high[258:1273],col="green")
lines(bbew2_low[258:1273],col="green")

lines(bbroll_med,col="purple")
lines(bbroll_high,col="orange")
lines(bbroll_low,col="orange")

##Beta Apart
plot(bbew_med[258:1273],type="l",ylim = c(0,2))
lines(bbew_high[258:1273],col="blue")
lines(bbew_low[258:1273],col="blue")

plot(bbew2_med[258:1273],col="red",ylim=c(0,2),type="l")
lines(bbew2_high[258:1273],col="green")
lines(bbew2_low[258:1273],col="green")

plot(bbroll_med,col="purple",ylim=c(0,2),type="l")
lines(bbroll_high,col="orange")
lines(bbroll_low,col="orange")
