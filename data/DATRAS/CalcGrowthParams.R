library(ggplot2)

## http://derekogle.com/fishR/examples/oldFishRVignettes/LengthWeight.pdf
DF <- read.csv('SMALK_EVHOE.csv')
DF <- DF[!is.na(DF$IndWgt),]

spp <- c('Gadus morhua','Merlangius merlangus','Melanogrammus aeglefinus')

DF <- DF[DF$Species %in% spp,]

ggplot(DF, aes(x = LngtClass, y = IndWgt)) + geom_point(aes(colour = factor(Year))) + facet_wrap(~Species, scale = 'free')

# W = a L^b
# log(Wt) = log(a) + b*log(L) + E
# weight is in g, length in mm

DF$lWt <- log(DF$IndWgt)
DF$lL  <- log(DF$LngtClass)

ggplot(DF, aes(x = lL, y = lWt)) + geom_point(aes(colour = factor(Year))) + facet_wrap(~Species, scale = 'free' )

lm1 <- glm(lWt ~ lL + Year + Species, data = DF)

summary(lm1)

ggplot(DF, aes(x = lL, y = lWt)) + geom_point(aes(colour = factor(Year))) + facet_wrap(~Species, scale = 'free' ) +
	geom_smooth(method = 'lm')

predDF <- expand.grid(lL = log(seq(min(DF$LngtClass), max(DF$LngtClass), 10)), Year = seq(min(DF$Year), max(DF$Year), 1), Species  = spp)

predDF$lWt <- predict(lm1, newdata = predDF)

predDF$L <- exp(predDF$lL)
predDF$Wt <- exp(predDF$lWt)

ggplot(DF, aes(x = LngtClass, y = IndWgt)) + geom_point() +
	facet_wrap(Species ~ Year) + geom_line(data = predDF, aes(x = L, y = Wt))

drop1(lm1)

#######################################################
## Year not massively important, given gaps drop it...
#######################################################
lm2 <- glm(lWt ~ lL + Species, data = DF)

predDF <- data.frame(lL = log(c(seq(min(DF$LngtClass[DF$Species == spp[1]]), max(DF$LngtClass[DF$Species == spp[1]]), l = 80),
			       seq(min(DF$LngtClass[DF$Species == spp[2]]), max(DF$LngtClass[DF$Species == spp[2]]), l = 80),
			       seq(min(DF$LngtClass[DF$Species == spp[3]]), max(DF$LngtClass[DF$Species == spp[3]]), l = 80))),
                     Species = rep(spp, each = 80))


predDF$lWt <- predict(lm2, newdata = predDF)

predDF$L <- exp(predDF$lL)
predDF$Wt <- exp(predDF$lWt)

## Now we need to bias correct due to the fact that the mean on the logscale is
## the geometric mean...

corr.fact <- exp(sigma(lm2)^2/2)
predDF$WtCorr <- predDF$Wt * corr.fact

ggplot(DF, aes(x = LngtClass, y = IndWgt)) + geom_point() +
	facet_wrap(~ Species, scale = 'free') + geom_line(data = predDF, aes(x = L, y = Wt), col = 'red') +
        geom_line(data = predDF, aes(x = L, y = WtCorr), col = 'blue')

## Save the fit and the correction factor
save(lm2, corr.fact, file = 'LengthWeightPredictGadoids.RData')

###################################################
## Liklihood approach to estimate a and b params ##
###################################################

# data
Wt = DF$IndWgt[DF$Species == 'Gadus morhua']
L =  DF$LngtClass[DF$Species == 'Gadus morhua'] 

reg.ll <- function (theta) {
	mean.vec <- theta[1] * (L^theta[2])
	ll.vec <- dnorm(Wt, mean = mean.vec, sd = theta[3], log = TRUE)
	nll <- sum(-ll.vec)
	return(nll)
}

## starting params
a <- 0.0002 
b <- 3
sigma <-  2 
par <- c(a, b, sigma)

fit <- optim(par, reg.ll, control = list(trace = 3))

plot(x = DF$LngtClass[DF$Species == 'Gadus morhua'], y = DF$IndWgt[DF$Species == 'Gadus morhua'])
x  <-  seq(min(DF$LngtClass[DF$Species == 'Gadus morhua']), max(DF$LngtClass[DF$Species == 'Gadus morhua']), 10)
lines(x = x,y = fit$par[1] * (x^fit$par[2]), col = 'red')
lines(x = predDF$L[predDF$Species== 'Gadus morhua'], y = predDF$Wt[predDF$Species == 'Gadus morhua'], col = 'blue')


## Fit doesn't 'look' as good...
