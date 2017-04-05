
load('Save.RData')

Save$Opt$SD$par.fixed[names(Save$Opt$SD$par.fixed) %in% c('lambda1_k','lambda2_k')]

Save$Opt$SD$sd[names(Save$Opt$SD$sd) %in% c('lambda1_k','lambda2_k')]


sd.sum <- summary(Save$Opt$SD, 'fixed', p.value = TRUE)

lambda1_k  <- sd.sum[rownames(sd.sum) %in% c('lambda1_k'),]

plot(lambda1_k[,'Estimate'] ~ c(1:6), ylim = c(-1,1))
points(lambda1_k[,'Estimate'] + 1.96 * lambda1_k[,'Std. Error'], col  = 'red')
points(lambda1_k[,'Estimate'] - 1.96 * lambda1_k[,'Std. Error'], col  = 'red')
