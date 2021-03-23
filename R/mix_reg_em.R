###Mixtures of Regressions
###09/20/2020

library(MCMCpack)

mix_reg_em = function(y, x, n_mix = 5, n_its = 5000){
  betas = array(0,c(dim(x)[2],n_mix,n_its+1))
  sigma = matrix(0.1,n_mix,n_its+1)
  for(a in 1:n_mix){
    for(b in 1:dim(x)[2]){
      betas[b,a,] = rep(a,n_its+1)
    }
    sigma[a,] = rep(a,n_its+1)
  }
  p_assign = matrix(1/n_mix,n_mix,n_its+1)
  assignment = matrix(0,dim(x)[1],n_its)
  for(i in 1:n_its){
    ##Sample assignment
    ##First compute likelihood elements
    llhe = matrix(0,n_mix,dim(x)[1])
    for(j in 1:n_mix){
      llhe[j,] = -0.5*(y - x %*% betas[,j,i])^2/sigma[j,i] -0.5*log(sigma[j,i]) - 0.5*log(2*pi)
      ##Accounts for priors
      llhe[j,] = llhe[j,] - sum(betas[,j,i]^2)/1000 - 1.5 * log(sigma[j,i]) - 0.5/sigma[j,i]
    }
    for(dimex  in 1:dim(x)[1]){
      llhe[,dimex] = exp(llhe[,dimex])*p_assign[,i]/sum(exp(llhe[,dimex])*p_assign[,i])
      assignment[dimex,i] = sample(1:n_mix,1,prob = llhe[,dimex])
    }
    #print(assignment[,i])
    ##Sample coefficients and errors
    for(sbg in 1:n_mix){
      ##Sample betass
      wx = x[which(assignment[,i]==sbg),]
      wy = y[which(assignment[,i]==sbg)]
      #print(dim(wx))
      if(!is.null(dim(wx)[1])){
        if(dim(wx)[1]>0){
          if(dim(x)[2]==1){
            sighat = 0.002+(1/sigma[sbg,i])*sum(wx^2)
            #print(sigmas[j,i-begin+1+1])
            sighat = 1/sighat
            sighat = matrix(sighat,1,1)
          }else{
            #print(dim(wx))
            sighat = 0.002*diag(dim(x)[2])+(1/sigma[sbg,i])*(t(wx)%*%wx)
            sighat = solve(sighat)
          }

          muhat = (1/sigma[sbg,i]) * (sighat %*% t(wx) %*% wy)
          if(dim(x)[2]==1){
            sighat = sqrt(sighat)
          }else{
            sighat = chol(sighat)
          }
          betas[,sbg,i+1] = muhat + sighat %*% rnorm(dim(x)[2])

          ##Sample Sigmas
          resids = wy-wx %*% matrix(betas[,sbg,i+1],ncol=1)
          #print(sum(resids^2))
          shas = 1 + 0.5*length(wy)
          ras  = 1 + 0.5*sum(resids^2)
          sigma[sbg,i+1] = 1/rgamma(1,shas,ras)
        }else{
          betas[,sbg,i+1] = sqrt(500) * rnorm(dim(x)[2])
          sigma[sbg,i+1] = 1/rgamma(1,0.5,0.5)
        }
        #print(sigma[sbg,i+1])
      }else{
        betas[,sbg,i+1] = sqrt(500) * rnorm(dim(x)[2])
        sigma[sbg,i+1] = 1/rgamma(1,0.5,0.5)
      }
    }

    ##Sample probabilities of assignment
    counts = rep(0,n_mix)
    for(ll in 1:n_mix){
      counts[ll] = length(which(assignment[,i]==ll))+1
    }
    p_assign[,i+1] = rdirichlet(1,counts)
    if(i%%100==0){
        print(paste(i,"draw(s) done!"))
    }

  }
  betas = betas[,,-(n_its+1)]
  p_assign = p_assign[,-(n_its+1)]
  sigma = sigma[,-(n_its+1)]

  betas = betas[,,-1]
  p_assign = p_assign[,-1]
  sigma = sigma[,-1]
  assignment = assignment[,-1]

return(list(betas = betas, sigma = sigma, Pi = p_assign, st = assignment))
}

reorder_draws_mr = function(out_draws,n_its=1000){
  new_list = out_draws
  for(i in 1:n_its){
    if(length(dim(out_draws$betas))==2){
      reorderit = sort(out_draws$betas[,i],index.return=TRUE)$ix
    }else{
      reorderit = sort(out_draws$betas[1,,i],index.return=TRUE)$ix
    }
    new_list$Pi[,i] = out_draws$Pi[reorderit,i]
    new_list$sigma[,i] = out_draws$sigma[reorderit,i]

    for(j in 1:dim(out_draws$st)[1]){
        new_list$st[,i] = reorderit[out_draws$st[j,i]]
    }
    if(length(dim(out_draws$betas))==2){
      new_list$betas[,i] = out_draws$betas[reorderit,i]
    }else{
      new_list$betas[,,i] = out_draws$betas[,reorderit,i]
    }
  }
  return(new_list)
}

##Function to compute the conditional mean

compute_cm = function(x,out_draws,n_its){
  ppd_y = matrix(0,dim(x)[1],n_its)
  for(i in 1:n_its){
    for(j in 1:dim(x)[1]){
      cur_st = out_draws$st[j,i]
      ppd_y[j,i] = sum(x[j,]*out_draws$betas[,cur_st,i])
    }
  }
  return(ppd_y)
}

x = matrix(rnorm(500),ncol=5)
coefs = matrix(c(1,0,2,0,-1),ncol=1)
coefs2 = matrix(c(0,2,0,-3,0),ncol=1)
y1 = x %*% coefs + rnorm(100)
y2 = x %*% coefs2 + rnorm(100)
y = rbind(y1,y2)
y = as.matrix(y)
x = rbind(x,x)
x = as.matrix(x)


gg_dat = read.csv("data_gdp_gr.csv")
gg_dat = gg_dat[,-1]
gg_ord = sort(gg_dat[,1],index.return=TRUE)$ix
gr = as.matrix(gg_dat[gg_ord,2])/sd(gg_dat[,2])
gdp_c = apply(as.matrix(cbind(rep(1,dim(gg_dat)[1]),gg_dat[gg_ord,1],gg_dat[gg_ord,1]^2)),2,function(x) x/sd(x))
gdp_c[,1] =1


xs = rnorm(100)
xs = (xs-min(xs))
xs = xs/max(xs)
xs = xs[sort(xs,index.return=TRUE)$ix]
mu_ys = 3*sin(4*pi*xs)
ys = mu_ys + 0.1*rnorm(100)

xs = cbind(rep(1,length(xs)),xs,xs^2)
xs = as.matrix(xs)
ys = as.matrix(ys)
