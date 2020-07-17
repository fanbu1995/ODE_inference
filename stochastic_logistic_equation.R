# 07/12/2020

# try my ODE inference framework on 
# logistic equation <--> stochastic logistic equation

library(tidyverse)

# 1. a function to simulate a SLE
# 07/17/2020: change the diffusion term : linear to Yt
sim_SLE <- function(y0=100, alpha=0.2, K=1000, sigma=1, timestep=0.1, tmax=50, plot=FALSE){
  times = seq(from=0, to=tmax, by=timestep)
  values = numeric(length=length(times))
  values[1] = y0
  
  Steps = length(times)-1
  y = y0
  SD = sigma*sqrt(timestep)
  for(s in 1:Steps){
    #delta.y = y*(1-y/K) * (alpha*timestep + rnorm(1, mean = 0, sd = SD))
    # 07/17/2020 change below
    delta.y = y*(1-y/K) * alpha*timestep + y * rnorm(1, mean = 0, sd = SD)
    y = y + delta.y
    values[s+1] = y
  }
  
  if(plot){
    plot(values~times, type="l")
  }
  
  return(list(times = times, values=values, tmax=tmax, timestep=timestep))
}

SLE1 = sim_SLE(plot=TRUE)

SLE2 = sim_SLE(alpha=0.5, sigma=0.2, plot = TRUE)

SLE3 = sim_SLE(alpha=0.5, sigma = 1.1, plot = TRUE)

# 2. a function to obtain discrete observations
discretize_SLE <- function(SLE_list, tau){
  timestep = SLE_list$timestep
  times = SLE_list$times
  values = SLE_list$values
  tmax = SLE_list$tmax
  
  ts = seq(from=0, to=tmax, by=tau)
  
  if(tau %% timestep){
    every = tau %/% timestep
    inds = seq(from=0, to=length(times), by=every)
    y = values[inds]
    return(list(times=ts, y=y))
  }else{
    obs = values[1]
    for(t in ts){
      if(t!=0){
        ## find the event time t_i such that t_i <= t < t_{i+1}
        ind = max(which(times <= t))
        y_t = values[ind]
        obs = c(obs, y_t)
      }
    }
    return(list(times=ts, y=obs))
  }
}

observed = discretize_SLE(SLE2, tau=1)
plot(observed$times, observed$y, type="l")


# 3. a function for iterative GLS procedure for logistic equation X
# (assume population capacity K is known...)
# 07/17/2020: LSE diffusion term linear to yt, 
# AND use closed form solution of SLE
iGLS_SLE <- function(observed, init.params, K=1000, tol=1e-4, maxIter=1000){
  # init.params = list(alpha=..., sigma2=...)
  times = observed$times
  y = observed$y
  y0 = x0 = y[1]
  dat = data.frame(times = times, y = y)
  N = nrow(dat)
  
  alphaA = NULL
  sigma2A = NULL
  
  #xformula = "K*x0/(x0+(K-x0)*exp(-alpha * times))"
  # 07/17/2020 change below
  xformula = "(alpha - sigma2/2)/(((alpha - sigma2/2)/x0 - alpha/K) * exp(-(alpha - sigma2/2) * times) + alpha/K)"
  
  # 0) initial estimates
  alpha = init.params$alpha
  sigma2 = init.params$sigma2
  #xvalues = K*x0/(x0+(K-x0)*exp(-alpha * times))
  # 07/17/2020 change below
  xvalues = (alpha - sigma2/2)/(((alpha - sigma2/2)/x0 - alpha/K) * exp(-(alpha - sigma2/2) * times) + alpha/K)
  ww = 1/xvalues
  #xformula = paste("exp(", as.character(init.param), " * times)", sep="") 
  init.fit = nls(as.formula(paste("y",xformula,sep=" ~ ")), data = dat,
                 start = list(alpha = alpha, sigma2 = sigma2), weights = ww,
                 control = nls.control(warnOnly = T))
  alpha.old = alpha
  #alpha = coef(init.fit) %>% as.numeric()
  parms = coef(init.fit)
  #cat(parms)
  alpha = parms[1]; sigma2 = parms[2]
  
  iter = 0
  while((abs(alpha - alpha.old) > tol) & (iter < maxIter)){
    iter = iter + 1
    # 1) estimate gamma
    xvalues = K*x0/(x0+(K-x0)*exp(-alpha * times))
    # 07/17/2020 change below
    # xvalues = (alpha - sigma2/2)/(((alpha - sigma2/2)/x0 - alpha/K) * exp(-(alpha - sigma2/2) * times) + alpha/K)
    gamma = sum((observed$y - xvalues)^2/xvalues)/(N-2)
    #sigma2 = sum((observed$y - xvalues)^2/xvalues)/(N-1)
    #sigma2A = c(sigma2A, sigma2)
    
    # 2) estimate lambda
    alpha.old = alpha
    ww = 1/(gamma * xvalues)
    this.fit = nls(as.formula(paste("y",xformula,sep=" ~ ")), data = dat,
                   start = list(alpha = alpha.old, sigma2 = sigma2), weights = ww,
                   control = nls.control(warnOnly = T))
    #alpha = coef(this.fit) %>% as.numeric()
    parms = coef(init.fit)
    #print(summary(this.fit))
    alpha = parms[1]; sigma2 = parms[2]
    alphaA = c(alphaA, alpha)
    sigma2A = c(sigma2A, sigma2)
    # 3) print info
    cat("Iteration", iter, ", alpha =", alpha, "sigma2 =", sigma2, "\n")
  }
  
  return(list(data = dat,
              alpha=alpha, sigma2=sigma2, 
              alphaList=alphaA, sigma2List = sigma2A))
}

#try1 = iGLS_SLE(observed, list(alpha=0.2, sigma2=1))

observed2 = discretize_SLE(SLE2,0.5)
try2 = iGLS_SLE(observed2, list(alpha=0.6, sigma2=0.05))


# 4. a function that tests the whole process
test_iGLS_SLE <- function(y0=100, alpha=0.5, K=1000, sigma=0.5, timestep=0.1, tmax=50, 
                          plot=TRUE, tau = 1, init.params, tol=1e-4, maxIter=1000){
  # 1) simulate
  events = sim_SLE(y0, alpha, K, sigma, timestep, tmax, plot)
  cat("simulation done!\n")
  # 2) discretize
  obs = discretize_SLE(events, tau)
  cat("Discrete data acquired with stepsize =",tau,"\n")
  if(plot){
    plot(obs$times, obs$y, type="l", main="discrete observations")
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
  res = iGLS_SLE(obs, init.params, K, tol, maxIter)
  # 4) compare with NLS
  x0 = obs$y[1]
  NLS.fit = nls(y ~ (alpha - sigma2/2)/(((alpha - sigma2/2)/x0 - alpha/K) * exp(-(alpha - sigma2/2) * times) + alpha/K), 
                res$data, 
                start = list(alpha = init.params$alpha, sigma2 = init.params$sigma2),
                control = nls.control(warnOnly = T))
  NLS.est = coef(NLS.fit)[1] %>% as.numeric()
  cat("NLS estimate for alpha: ", NLS.est,"\n")
  
  # 5) and true lambda value and NLS estimate
  res$alphaTrue = alpha
  res$alphaNLS = NLS.est
  
  return(res)
}

res1 = test_iGLS_SLE(init.params = list(alpha=0.3, sigma2=0.3))
res2 = test_iGLS_SLE(alpha=1, sigma=0.5, tau=0.2, init.params = list(alpha=0.8, sigma2=0.2))


# 5. Conduct a series of experiments with
# alpha = 0.5
# sigma = 0.2, 0.5, 1
# tau = 0.2, 1
# each repeat 20 times?

alpha=1
Sigmas = c(0.2, 0.5, 0.8)
#Sigmas = c(0.5, 1, 1.1)
Taus = c(0.2, 1)

repeat_test <- function(SigmaList, TauList, R=20, showCurves=TRUE, ...){
  resTable = NULL
  for(sigma in SigmaList){
    for(tau in TauList){
      cat("Experiment with sigma =", sigma, "and tau =", tau,"...\n")
      if(showCurves){
        allData = NULL
      }
      for(r in 1:R){
        this.test = test_iGLS_SLE(sigma=sigma, tau=tau, ...)
        if(showCurves){
          this.data = this.test$data
          this.data$rep = r
          allData = rbind(allData, this.data)
        }
        est.error = this.test$alpha - this.test$alphaTrue
        NLS.error = this.test$alphaNLS - this.test$alphaTrue
        resTable = rbind(resTable, c(sigma, tau, est.error))
        resTable = rbind(resTable, c(sigma, tau, NLS.error))
        cat(r,"\t")
      }
      cat("Experiment done and results recorded!\n")
      if(showCurves){
        print(
          ggplot(allData, aes(x=times, y=y, group = rep)) +
            geom_line(size=0.2) +
            labs(x="time", y="y", title="discrete observations",
                 caption = paste("alpha=", alpha, "sigma=", sigma,"tau=",tau)) +
            theme_bw(base_size = 14)
        )
      }
    }
  }
  resTable = as.data.frame(resTable)
  names(resTable) = c("sigma","tau", "error")
  n = length(SigmaList) * length(TauList) * R
  resTable$method = rep(c("Approx","NLS"), n)
  
  return(resTable)
}


repTest1 = repeat_test(Sigmas, Taus, alpha=alpha, plot=FALSE, 
                       init.params = list(alpha=0.6, sigma2=0.6))

repTest2 = repeat_test(Sigmas, Taus, alpha=alpha, plot=FALSE, 
                       init.params = list(alpha=0.5, sigma2=1))

## Findings:
## 1. For some sequences, the algorithm doesn't convergence: 
##    two (or more) solutions fit the data equally well, so the algorithm oscillates between them...
## 2. My framework doesn't really work that well, because
##   a) for this SDE, the ODE isn't the mean function; bistability of SDE isn't present in ODE system
##   b) it works comparably as NLS, so it doesn't solve the real issue; 
##      numerical problems encountered when noise level > growth rate

#saveRDS(repTest1, "SLE_simulation_res.rds")
saveRDS(repTest1, "SLE_simulation_res_updated.rds")

ggplot(data=repTest1, aes(x=method,y=error)) +
  geom_hline(yintercept = 0, size=1, color="gray") +
  geom_boxplot() +
  #geom_violin() +
  theme_bw(base_size = 14)+
  facet_grid(tau~sigma)

errors_summary = 
repTest1 %>% group_by(sigma, tau, method) %>% 
  summarise(mean_error = mean(error), 
            error_Q1 = quantile(error, 0.25),
            error_Q3 = quantile(error, 0.75))

library(xtable)

xtable(errors_summary,digits = 3)
