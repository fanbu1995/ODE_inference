# 07/15/2020
# try inference framework on Lotka-Volterra model

library(deSolve)
library(tidyverse)
library(minpack.lm)

# 1. ODE definition of LV model
# 07/17/2020: change to beta=delta case
LV_ODE <- function(t, y, parms){
  x1 = y[1]
  x2 = y[2]
  with(as.list(parms),{
    dx1 = alpha * x1 - beta * x1 * x2
    #dx2 = delta * x1 * x2 - x2
    dx2 = beta * x1 * x2 - x2
    list(c(dx1, dx2))
  })
}

# try it out (and obtain numerical solution over time)
times = seq(from=0, to=20, by=1)
y0 = c(35,50)
#p = c(alpha=2, beta=1, delta=0.5)
p = c(alpha=0.5, beta=0.02)

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
### after re-parameterization, it kind of works!
matplot(res2[,1], res2[,2:3], type="l", 
        xlab="time", ylab="population",
        main = "LV dynamics: ODE")

res3 = ode(y0, times, LV_ODE, p, method="ode45")
head(res3)
### it also stops early... maybe because the prey population can't drop below zero?
### after re-parameterization, it kind of works!
matplot(res3[,1], res3[,2:3], type="l", 
        xlab="time", ylab="population",
        main = "LV dynamics: ODE")

# 2. CTMC simulator
# 07/17/2020: outdated... right now try GillespieSSA package
# sim_LV <- function(y0, rates, tmax = 100){
#   # rates: should be a named vector
#   alpha = rates['alpha']
#   beta = rates['beta']
#   delta = rates['delta']
#   
#   values = y0 
#   times = 0
#   
#   x1 = y0[1]; x2 = y0[2]
#   t = 0
#   
#   while(t < tmax){
#     ## calculate rates
#     ## (a little adjustment, in case x1 or x2 < 0)
#     cat("x1=",x1,"x2=",x2,"\t")
#     if(x1 <= 0){
#       if(x2 <= 0){
#         R = c(alpha, 1e-3, 1e-3)
#       }else{
#         R = c(alpha, 1e-3, x2)
#       }
#     }else if(x2 <= 0){
#       R = c(alpha * x1, beta * x1, 1e-3)
#     }else{
#       R = c(alpha * x1, beta * x1 * x2, x2)
#     }
#     totalR = sum(R)
#     
#     cat('Rates=',R,"\n")
#     
#     ## next event time
#     delta.t = rexp(1, rate = totalR)
#     next.t = t + delta.t
#     times = c(times, next.t)
#     
#     ## next event type
#     z = sample(3, 1, replace = FALSE, prob = R)
#     if(z == 1){
#       ## x1 growth
#       x1 = x1 + 1
#     }else if(z == 2){
#       ## x1 decrease and x2 increase
#       x1 = x1 - 1
#       x2 = x2 + delta/beta #(predator grows proportionally)
#     }else{
#       ## x2 death
#       x2 = x2 - 1
#     }
#     values = rbind(values, c(x1,x2))
#     
#     ## forward clock
#     t = next.t
#   }
#   
#   return(list(values = values, times = times, 
#               params = rates, tmax = tmax))
# }
# 
# ## try it out
# CTMC_res = sim_LV(y0, p)
# 
# matplot(CTMC_res$times, CTMC_res$values, type = "l",
#         xlab="time", ylab="population",
#         main = "LV dynamics: CTMC")

# 2.5: use GillespieSSA package to simulate
library(GillespieSSA)

parms <- c(alpha = 0.5, beta = 0.02, gamma = 1)
tf <- 20                                     # Final time
simName <- "Lotka-Volterra CTMC"              # Name

x0 <- c(Y1=35, Y2=50)

nu <- matrix(c(+1, -1, 0, 0, 1, -1), nrow = 2, byrow = TRUE)

a  <- c("alpha*Y1", "beta*Y1*Y2","gamma*Y2")

set.seed(1)
out <- ssa(
  x0 = x0,
  a = a,
  nu = nu,
  parms = parms,
  tf = tf,
  method = ssa.d(),
  simName = simName,
  verbose = FALSE,
  consoleInterval = 1
) 
ssa.plot(out, show.title = TRUE, show.legend = TRUE)

# 2.5.x get discrete observations from SSA output
get_discrete_SSA <- function(SSA_out, tau){
  times = SSA_out$data[,"t"]
  values = SSA_out$data[, 2:3]
  tmax = SSA_out$args$tf
  
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
observed2 = get_discrete_SSA(out, 1)
matplot(observed2$times, observed2$y, type = "l")



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
# 07/17/2020: change to beta=delta, still gamma=1
residFun <- function(params, observed, y0, weights=NULL){
  # weights: a n*2 matrix of inverse weights = gamma * x_i(t)
  alpha = params[1]
  beta = params[2]
  #delta = params[3]
  
  times = observed$times
  values = observed$y
  
  pred = ode(y0, times, LV_ODE, list(alpha=alpha, beta=beta))
  
  error = pred[,2:3] - values
  if(is.null(weights)){
    Res = sqrt(error[,1]^2 + error[,2]^2)
  }else{
    Res = sqrt(error[,1]^2/weights[,1] + error[,2]^2/weights[,2])
  }
  
  # cat("nrow of pred:", nrow(pred), "; nrow of values", nrow(values), "\n")
  # error = pred[,2] - values[which(times %in% pred[,1]),1]
  # if(is.null(weights)){
  #   Res = error
  # }else{
  #   Res = error * sqrt(1/weights[,2])
  # }
  
  return(Res)
}

## try it out
LM_res = nls.lm(par=c(alpha=1, beta=2, delta=0.5), fn = residFun, 
                observed = observed, y0 = c(50,20))
summary(LM_res)
# it doesn't really work!!
# the estimates don't move far from the starting points

# 07/17/2020
## try it with the modified setting
## (also play around a bit with the weights to see if it works)
LM_res = nls.lm(par=c(alpha=0.3, beta=0.01), fn = residFun, 
                observed = observed2, y0 = c(35,50), 
                weights = matrix(runif(42), nrow = 21, ncol = 2))
summary(LM_res)
# estimation isn't very accurate but at least the algorithm is doing sth.


# 5. a function that runs the iGLS procedure
iGLS_LV <- function(observed, init.params, tol=1e-4, maxIter=1000){
  # init.params = list(alpha=..., beta=...)
  times = observed$times
  y = observed$y
  y0 = y[1,]
  N = length(times)
  
  alphaA = NULL
  betaA = NULL
  gammas = numeric(2)
  
  # 0) initial estimates
  alpha = init.params$alpha
  beta = init.params$beta
  
  ODE_solutions = ode(y0, times, LV_ODE, 
                      c(alpha=alpha, beta=beta))
  xvalues = ODE_solutions[,2:3]
  ## a bit adjust to avoid 1/0 situations
  xvalues[xvalues < 0.1] = 0.1
  ww = 1/xvalues
  init.fit = nls.lm(par=c(alpha=init.params$alpha, beta=init.params$beta), 
                    fn = residFun, 
                    observed = observed, y0 = y0, weights=ww)
  
  alpha.old = alpha; beta.old = beta
  params = coef(init.fit) %>% as.numeric()
  alpha = params[1]; beta = params[2]
  
  iter = 0
  while((abs(alpha - alpha.old) > tol | abs(beta - beta.old) > tol) & (iter < maxIter)){
    iter = iter + 1
    # 1) estimate gamma
    ODE_solutions = ode(y0, times, LV_ODE, 
                        c(alpha=alpha, beta=beta))
    xvalues = ODE_solutions[,2:3]
    gammas[1] = sum((y[,1] - xvalues[,1])^2/xvalues[,1])/(N-2)
    gammas[2] = sum((y[,2] - xvalues[,2])^2/xvalues[,2])/(N-2)
    
    #sigma2 = sum((observed$y - xvalues)^2/xvalues)/(N-1)
    #sigma2A = c(sigma2A, sigma2)
    
    # 2) estimate lambda
    alpha.old = alpha; beta.old = beta
    ## a bit adjust to avoid 1/0 situations
    xvalues[xvalues < 0.1] = 0.1
    ww = 1/t(gammas * t(xvalues))
    
    alpha.old = alpha; beta.old = beta
    this.fit = nls.lm(par=c(alpha=alpha.old, beta=beta.old), 
                      fn = residFun, 
                      observed = observed, y0 = y0, weights=ww)
    
    params = coef(init.fit) %>% as.numeric()
    alpha = params[1]; beta = params[2]
    
    alphaA = c(alphaA, alpha)
    betaA = c(betaA, beta)
    # 3) print info
    cat("Iteration", iter, ", alpha =", alpha, "beta =", beta, "\n")
  }
  
  return(list(data = data.frame(t=times, y1 = y[,1], y2 = y[,2]),
              #alpha=alpha, beta=beta, 
              est=c(alpha=alpha, beta=beta),
              alphaList=alphaA, betaList = betaA))
}

# try it out
try1 = iGLS_LV(observed2, list(alpha=0.5, beta=0.02))

# 6. a function for the entire pipeline
test_iGLS_VL <- function(y0=c(35,50), alpha=0.5, beta=0.02, tmax=20, 
                          plot=TRUE, tau = 1, init.params, tol=1e-4, maxIter=1000){
  # 1) simulate                                   
  # set.seed(1)
  events <- ssa(
    x0 = c(Y1=y0[1], Y2=y0[2]),
    a = c("alpha*Y1", "beta*Y1*Y2","gamma*Y2"),
    nu = matrix(c(+1, -1, 0, 0, 1, -1), nrow = 2, byrow = TRUE),
    parms = c(alpha = alpha, beta = beta, gamma = 1),
    tf = tmax,
    method = ssa.d(),
    simName = "Lotka-Volterra CTMC",
    verbose = FALSE,
    consoleInterval = 1
  ) 
  # if(plot){
  #   ssa.plot(events, show.title = TRUE, show.legend = TRUE)
  # }
  cat("simulation done!\n")
  # 2) discretize
  obs = get_discrete_SSA(events, tau)
  cat("Discrete data acquired with stepsize =",tau,"\n")
  if(plot){
    matplot(obs$times, obs$y, type="l", main="discrete observations")
  }
  # # 2.5) if simulated result is a flat line, then repeat
  # while(all(obs$y == 1)){
  #   # 1) simulate
  #   events = sim_branching(lambda, tmax, y0)
  #   cat("simulation done!\n")
  #   # 2) discretize
  #   obs = get_discrete_obs(events, tau)
  #   cat("Discrete data acquired with stepsize =",tau,"\n")
  #   plot(obs$times, obs$y)
  # }
  # 3) estimate
  res = iGLS_LV(obs, init.params, tol, maxIter)
  # 4) compare with NLS
  NLS.fit = nls.lm(par=c(alpha=init.params$alpha, beta=init.params$beta), 
                   fn = residFun, 
                   observed = obs, y0 = y0)
  
  NLS.est = coef(NLS.fit)
  cat("NLS estimate: ", NLS.est,"\n")
  
  # 5) and true lambda value and NLS estimate
  res$paramTrue = c(alpha=alpha, beta=beta)
  res$estNLS = NLS.est
  
  return(res)
}

test_res = test_iGLS_VL(alpha=1, beta=0.02, 
                        init.params = list(alpha=1.1, beta=0.015))


# 7. function for repeated test
Alphas= c(0.5, 1)
Betas = c(0.02)
Taus = c(0.2, 1)

repeat_test <- function(AlphaList, BetaList, TauList, R=20, ...){
  resTable = NULL
  for(alpha in AlphaList){
    for(beta in BetaList){
      for(tau in TauList){
        cat("Experiment with alpha =", alpha, "beta =", beta,
            "and tau =", tau, "...\n")
        for(r in 1:R){
          this.test = test_iGLS_VL(alpha=alpha, beta=beta, tau=tau, 
                                   init.params = list(alpha=alpha+0.1, beta=beta+0.01),
                                   ...)
          est.error = this.test$est - this.test$paramTrue
          NLS.error = this.test$estNLS - this.test$paramTrue
          resTable = rbind(resTable, c(alpha, beta, tau, est.error))
          resTable = rbind(resTable, c(alpha, beta, tau, NLS.error))
          cat(r,"\t")
        }
        cat("Experiment done and results recorded!\n")
      }
    }
  }
  resTable = as.data.frame(resTable)
  names(resTable) = c("alpha","beta","tau", "alpha.error", "beta.error")
  n = length(AlphaList) * length(BetaList) * length(TauList) * R
  resTable$method = rep(c("Approx","NLS"), n)

  return(resTable)
}

repRes1 = repeat_test(Alphas, Betas, Taus)

# fix a bug of true beta value manually...
repRes1$beta.error = repRes1$beta.error + repRes1$alpha - repRes1$beta

saveRDS(repRes1, "Lotka_Volterra_erros.rds")

ggplot(data=repRes1, aes(x=method,y=alpha.error)) +
  geom_hline(yintercept = 0, size=1, color="gray") +
  geom_boxplot() +
  #geom_violin() +
  theme_bw(base_size = 14)+
  facet_grid(tau~alpha)

ggplot(data=repRes1, aes(x=method,y=beta.error)) +
  geom_hline(yintercept = 0, size=1, color="gray") +
  geom_boxplot() +
  #geom_violin() +
  theme_bw(base_size = 14)+
  facet_grid(tau~alpha)
