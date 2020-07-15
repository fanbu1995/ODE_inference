# 07/15/2020
# try inference framework on Lotka-Volterra model

library(deSolve)
library(tidyverse)
library(minpack.lm)

# 1. ODE definition of LV model
LV_ODE <- function(t, y, parms){
  x1 = y[1]
  x2 = y[2]
  with(as.list(parms),{
    dx1 = alpha * x1 - beta * x1 * x2
    dx2 = delta * x1 * x2 - x2
    list(c(dx1, dx2))
  })
}

# try it out (and obtain numerical solution over time)
times = seq(from=0, to=100, by=1)
y0 = c(50,20)
p = c(alpha=2, beta=1, delta=0.5)

res = ode(y0, times, LV_ODE, p)
head(res)

matplot(res[,1], res[,2:3], type="l", 
        xlab="time", ylab="population",
        main = "LV dynamics: ODE")

## try out a "stiff ODE solver"
res2 = ode(y0, times, LV_ODE, p, method="ode23")
head(res2)
### hmmm.. it stops integrating after a few time steps
### AND it's really slow!

res3 = ode(y0, times, LV_ODE, p, method="ode45")
head(res3)
### it also stops early... maybe because the prey population can't drop below zero?


# 2. CTMC simulator
sim_LV <- function(y0, rates, tmax = 100){
  # rates: should be a named vector
  alpha = rates['alpha']
  beta = rates['beta']
  delta = rates['delta']
  
  values = y0 
  times = 0
  
  x1 = y0[1]; x2 = y0[2]
  t = 0
  
  while(t < tmax){
    ## calculate rates
    ## (a little adjustment, in case x1 or x2 < 0)
    cat("x1=",x1,"x2=",x2,"\t")
    if(x1 <= 0){
      if(x2 <= 0){
        R = c(alpha, 1e-3, 1e-3)
      }else{
        R = c(alpha, 1e-3, x2)
      }
    }else if(x2 <= 0){
      R = c(alpha * x1, beta * x1, 1e-3)
    }else{
      R = c(alpha * x1, beta * x1 * x2, x2)
    }
    totalR = sum(R)
    
    cat('Rates=',R,"\n")
    
    ## next event time
    delta.t = rexp(1, rate = totalR)
    next.t = t + delta.t
    times = c(times, next.t)
    
    ## next event type
    z = sample(3, 1, replace = FALSE, prob = R)
    if(z == 1){
      ## x1 growth
      x1 = x1 + 1
    }else if(z == 2){
      ## x1 decrease and x2 increase
      x1 = x1 - 1
      x2 = x2 + delta/beta #(predator grows proportionally)
    }else{
      ## x2 death
      x2 = x2 - 1
    }
    values = rbind(values, c(x1,x2))
    
    ## forward clock
    t = next.t
  }
  
  return(list(values = values, times = times, 
              params = rates, tmax = tmax))
}

## try it out
CTMC_res = sim_LV(y0, p)

matplot(CTMC_res$times, CTMC_res$values, type = "l",
        xlab="time", ylab="population",
        main = "LV dynamics: CTMC")


# 3. get discrete observations
get_discrete_obs <- function(LV_list, tau){
  
  times = LV_list$times
  values = LV_list$values
  tmax = LV_list$tmax
  
  obs = values[1,] %>% as.numeric()
  ts = seq(from=0, to=tmax, by=tau)
  
  for(t in ts){
    if(t!=0){
      ## find the event time t_i such that t_i <= t < t_{i+1}
      ind = max(which(times <= t))
      y_t = values[ind,] %>% as.numeric()
      obs = rbind(obs, y_t)
    }
  }
  
  return(list(times=ts, y=obs))
}

## try it out
observed = get_discrete_obs(CTMC_res, 1)
matplot(observed$times, observed$y, type = "l")

# dat = data.frame(t = observed$times, y1=observed$y[,1], y2=observed$y[,2])
# 
# nls.fit = nls( ~ ode(c(50,20), t, LV_ODE, params=c(alpha,beta,delta)),
#               data = dat,
#               start=c(alpha=5,beta=1,delta=1))


# 4. define residual function used for LM NLS optimization function
residFun <- function(params, observed, y0, weights=NULL){
  # weights: a n*2 matrix of inverse weights = gamma * x_i(t)
  alpha = params[1]
  beta = params[2]
  delta = params[3]
  
  times = observed$times
  values = observed$y
  
  pred = ode(y0, times, LV_ODE, list(alpha=alpha, beta=beta, delta=delta))
  
  # error = pred - values
  # if(is.null(weights)){
  #   Res = sqrt(error[,2]^2 + error[,3]^2)
  # }else{
  #   Res = sqrt(error[,2]^2/weights[,1] + error[,3]^2/weights[,2])
  # }
  
  cat("nrow of pred:", nrow(pred), "; nrow of values", nrow(values), "\n")
  error = pred[,3] - values[which(times %in% pred[,1]),2]
  if(is.null(weights)){
    Res = error
  }else{
    Res = error[,2] * sqrt(1/weights[,2])
  }
  
  return(Res)
}

## try it out
LM_res = nls.lm(par=c(alpha=1, beta=2, delta=0.5), fn = residFun, 
                observed = observed, y0 = c(50,20))
summary(LM_res)
# it doesn't really work!!
# the estimates don't move far from the starting points
