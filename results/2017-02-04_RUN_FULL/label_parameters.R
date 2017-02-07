#############################################################
### Simple script to label up the parameters from the run  ##
#############################################################
load('parameter_estimates.RData')

library(dplyr)

params <- parameter_estimates$diagnostics

class(params)

n_params <- print(group_by(params, Param) %>% summarise(n()))
print(sum(n_params[,2])) 

## So....
## 1522 parameters total
##
# 468 beta1_ct      = 18 spp * 26 years - year effect encounter per species
# 468 beta2_ct      = 18 spp * 26 years - year effect +ve catch per species
# 6   gamma1_j      = Static covariates:  6 habitat factors (-1) + depth (cont) - encounter. 
# 6   gamma2_j      = Static covariates: 6 habitat factors (-1) + depth (cont)  - +ve catch
# 6   lambda1_k     = catchability variables : 7 gears - 1 - encounter
# 6   lambda2_k     = catchability variables: 7 gears -1 - +ve catch
# 135 L_epsilon1_z  = Spatio-temporal loading encounter - lower diag of loading matrix!!!
# 135 L_epsilon2_z  = Spatio-temporal loading +ve catch - lower diag of loading matrix!!!
# 2   ln_H_input    = Anisotropic values (encounter / +ve catch)
# 1   logkappa1     = Range for correlation distance matern encounter
# 1   logkappa2     = Range for correlation distance matern +ve catch 
# 18  logSigmaM     = 18 (species) - variance for species +ve catch ??
# 135 L_omega1_z    = Spatial loading encounter - lower diag of loading matrix!!!
# 135 L_omega2_z    = Spatial loading +ve catch - lower diag of loading matrix!!!


