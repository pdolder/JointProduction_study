library(TMB)
library(VAST)

load('Save.RData')

dyn.load(dynlib('VAST_v2_1_0'))

## Extract the hessian
Hess <- optimHess(Save$Obj$par, Save$Obj$fn, Save$Obj$gn)

# Invert  numerically
library(Corpcor)

p <- pseudoinverse(Hess)

save(Hess, p, file = 'InvHess.RData')


