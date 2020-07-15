# 07/15/2020
# try inference framework on Lotka-Volterra model

library(deSolve)
library(tidyverse)

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

matplot(res[,1], res[,2:3], type="l")

