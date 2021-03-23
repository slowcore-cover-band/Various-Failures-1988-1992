###Predictor Corrector Terminal Velocity
###02/27/2021




pred_cor = function(vel0,Bcoef,h_val,num_steps){
  velos = rep(0,num_steps+1)
  times = c(0:num_steps)*h_val
  grav = 9.8

  v0 = vel0
  velos[1] = v0

  for(i in 1:num_steps){
    dvdt0 = grav - Bcoef*v0^2
    v1 = v0+h_val*dvdt0
    v0 = v0+grav*h_val-0.5*h_val*Bcoef*(v0^2+v1^2)
    velos[i+1] = v0
  }
  return(list(velos=velos,times = times))
}

Cd = 1.2
A = 0.25
rho = 1.225
mass = 91

Binnit = rho*Cd*A/(2*mass)

terminal_velo = sqrt(9.8/Binnit)

vpath1 = pred_cor(0,Binnit,0.01,100)
