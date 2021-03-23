###Particle Slice Sampling
###09/13/20

particle_slice = function(x0,llh,scale=5,particles,num_draws){
  draws = matrix(0,num_draws+1,length(x0)+1)
  draws[1,] = c(x0,-llh(x0))
  for(i in 1:num_draws){
    for(j in 1:length(x0)){
      x_prop = x0[j] + rnorm(particles)*scale
      llh_part = rep(0,particles)
      for(k in 1:particles){
        if(j==1){
          curpar = c(x_prop[k],x0[-j])
        }else if(j==particles){
          curpar = c(x0[-j],x_prop[k])
        }else{
          curpar = c(x0[1:(j-1)],x_prop[k],x0[(j+1):particles])
        }
        llh_part[k] = -llh(curpar)
      }
      #print(x_prop)
      #print(llh_part)
      llh_part = exp(llh_part)/sum(exp(llh_part))
      x0[j] = sample(x_prop,1,prob=llh_part)
    }
  draws[i+1,] = c(x0,-llh(x0))
  }
  return(draws)
}

rosenbrock_llh = function(pars){
  llh = log((1-pars[1])^2 + 100*(pars[2]-pars[1]^2)^2+1)
  return(llh)
}
