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

predDF <- expand.grid(lL = log(seq(min(DF$LngtClass), max(DF$LngtClass), 10)), Species  = spp)

predDF$lWt <- predict(lm2, newdata = predDF)

predDF$L <- exp(predDF$lL)
predDF$Wt <- exp(predDF$lWt)

ggplot(DF, aes(x = LngtClass, y = IndWgt)) + geom_point() +
	facet_wrap(~ Species ) + geom_line(data = predDF, aes(x = L, y = Wt), col = 'red')

save(lm2, file = 'LengthWeightPredictGadoids.RData')
