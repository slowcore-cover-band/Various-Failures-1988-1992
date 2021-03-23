##Markov Switching CAPM
##09/08/20

##Estimates CAPM model with regime switches

library(MCMCpack)

ms_capm = function(y, x, regimes = 3, num_draw = 1000, tune_diag = 100){

  betas = array(0,c(dim(x)[2],regimes,num_draw+1))
  sigmas = matrix(0,regimes,num_draw+1)
  for(a in 1:regimes){
    for(b in 1:dim(x)[2]){
      betas[b,a,] = rep(a,num_draw+1)
    }
    sigmas[a,] = rep(a,num_draw+1)
  }
  if(dim(x)[2]==1){
    Bb = matrix(1,1,1)
  }else{
    Bb = matrix(1,dim(x)[2],1)
    Bb[1] = 0
  }
  Pis = array(1/regimes,c(regimes,regimes,num_draw+1))
  states = matrix(0,dim(x)[1],num_draw)
  p_filter = array(0,c(dim(x)[1]+1,regimes,num_draw))
  p_smooth = array(0,c(dim(x)[1],regimes,num_draw))
  for(i in 1:num_draw){
    #####Sample States
    ed_pi = eigen(Pis[,,i])
    stat_dis = ed_pi$vectors %*% diag(ed_pi$values^1000) %*% t(ed_pi$vectors) %*% matrix(0.5,regimes,1)
    ##print(stat_dis)
    stat_dis = Re(stat_dis)
    states[1,i] = sample(c(1:regimes),1,prob = stat_dis)
    p_filter[1,,i] = stat_dis
    llhe = matrix(0,regimes,dim(x)[1])
    for(j in 1:regimes){
      llhe[j,] = -0.5*(y - x %*% betas[,j,i])^2/sigmas[j,i] -0.5*log(sigmas[j,i]) - 0.5*log(2*pi)
      ##Accounts for priors
      llhe[j,] = llhe[j,] - sum(betas[,j,i]^2)/1000 - 1.1 * log(sigmas[j,i]) - 0.1/sigmas[j,i]
    }
    ##Filtering Stage
    for(k in 1:dim(x)[1]){
      pyy = 0
      for(l in 1:regimes){
        for(m in 1:regimes){
          pyy = pyy + p_filter[k,m,i]*exp(llhe[l,k])*Pis[m,l,i]
        }
      }
      for(n in 1:regimes){
        pmst = 0
        for(p in 1:regimes){
           pmst = pmst + p_filter[k,p,i]*exp(llhe[n,k])*Pis[p,n,i]
        }
        p_filter[k+1,n,i] = pmst/pyy
      }
    }
    ##Sampling Stage
    states[dim(x)[1],i] = sample(c(1:regimes),1,prob = p_filter[dim(x)[1]+1,,i])
    p_smooth[dim(x)[1],,i] = p_filter[dim(x)[1]+1,,i]
    for(q in (dim(x)[1]-1):1){
      ##Compute Smoothed Probabilities
      num = rep(0,regimes)
      for(cspa in 1:regimes){
        num[cspa] = Pis[cspa,states[q+1,i],i]*p_filter[q+1,cspa,i]
      }
      num = num/sum(num)
      states[q,i] = sample(c(1:regimes),1,prob = num)
      p_smooth[q,,i] = num
    }
    ##For data
    for(sbg in 1:regimes){
      ##Sample Betas
      wx = x[which(states[,i]==sbg),]
      wy = y[which(states[,i]==sbg)]
      if(dim(x)[2]==1){
        sighat = 0.002+(1/sigmas[sbg,i])*sum(wx^2)
        #print(sigmas[j,i-begin+1+1])
        sighat = 1/sighat
        sighat = matrix(sighat,1,1)
      }else{
        sighat = 0.002*diag(dim(x)[2])+(1/sigmas[sbg,i])*(t(wx)%*%wx)
        sighat = solve(sighat)
      }

      muhat = sighat %*% Bb + (1/sigmas[sbg,i]) * (sighat %*% t(wx) %*% wy)
      if(dim(x)[2]==1){
        sighat = sqrt(sighat)
      }else{
        sighat = chol(sighat)
      }
      betas[,sbg,i+1] = muhat + sighat %*% rnorm(dim(x)[2])

      ##Sample Sigmas
      resids = wy-wx %*% matrix(betas[,sbg,i+1],ncol=1)
      #print(sum(resids^2))
      shas = 0.1 + 0.5*length(wy)
      ras  = 0.1 + 0.5*sum(resids^2)
      sigmas[sbg,i+1] = 1/rgamma(1,shas,ras)
    }

    ##Sample Transition Matrix
    for(stm in 1:regimes){
      curval = which(states[,i]==stm)+1
      curval = curval[-which(curval==(dim(x)[1]+1))]
      stp1 = states[,i]
      counts = matrix(0,regimes,regimes)
      for(count_reg in 1:regimes){
        counts[count_reg,stm]=length(which(stp1 == count_reg))+1+as.numeric(count_reg==stm)*tune_diag
      }
      Pis[stm,,i+1] = rdirichlet(1,counts[,stm])
      #Compute likelihood
      ed_pi1 = eigen(Pis[,,i+1])
      pbar = ed_pi1$vectors %*% diag(ed_pi1$values^1000) %*% t(ed_pi1$vectors) %*% matrix(0.5,regimes,1)
      pbar = Re(pbar)
      llhn = sum(counts * log(Pis[,,i+1])) + log(pbar[states[1,i]])
      llho = sum(counts * log(Pis[,,i+1])) + log(stat_dis[states[1,i]])
      ##Rejection step to account for probability of first state
      if((llhn-llho)<log(runif(1))){
        Pis[,,i+1] = Pis[,,i]
      }
    }
    if(i%%100==0){
        print(paste(i,"draw(s) done!"))
    }

  }
  p_filter=p_filter[-1,,]
  ##Betas,sigmas,pis
  betas = betas[,,-(num_draw+1)]
  Pis = Pis[,,-(num_draw+1)]
  sigmas = sigmas[,-(num_draw+1)]
  return(list(Pi = Pis, beta = betas, sigma = sigmas, pf = p_filter, ps = p_smooth, st = states))
}
reorder_draws = function(out_draws,num_draw=1000){
  new_list = out_draws
  for(i in 1:num_draw){
    if(length(dim(out_draws$beta))==2){
      reorderit = sort(out_draws$beta[,i],index.return=TRUE)$ix
    }else{
      reorderit = sort(out_draws$beta[1,,i],index.return=TRUE)$ix
    }
    new_list$Pi[,,i] = out_draws$Pi[reorderit,reorderit,i]
    new_list$sigma[,i] = out_draws$sigma[reorderit,i]
    new_list$pf[,,i] = out_draws$pf[,reorderit,i]
    new_list$ps[,,i] = out_draws$ps[,reorderit,i]
    for(j in 1:dim(out_draws$st)[1]){
        new_list$st[,i] = reorderit[out_draws$st[j,i]]
    }
    if(length(dim(out_draws$beta))==2){
      new_list$beta[,i] = out_draws$beta[reorderit,i]
    }else{
      new_list$beta[,,i] = out_draws$beta[,reorderit,i]
    }
  }
  return(new_list)
}

weight_beta = function(coef,probs){
  weighteds = matrix(0,dim(probs)[1],dim(probs)[3])
  for(i in 1:dim(probs)[1]){
    weighteds[i,] = apply(probs[i,,] * coef,2,sum)
  }
  return(weighteds)
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


x = rnorm(100)
y = 1.5*x + 2*as.numeric(1:100>20)*x+ 3*as.numeric(1:100>50)*x + rnorm(100)

y = matrix(y,ncol=1)
x = matrix(x,ncol=1)

ms_test = ms_capm(y,x)
ms_test_ro = reorder_draws(ms_test)


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

ms_test = ms_capm(ersb,er5c)
ms_test_ro = reorder_draws(ms_test)

ms_test_int = ms_capm(ersb,er5cint)
ms_test_ro_int = reorder_draws(ms_test_int)
