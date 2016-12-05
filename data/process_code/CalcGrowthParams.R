
## http://derekogle.com/fishR/examples/oldFishRVignettes/LengthWeight.pdf
DF <- read.csv('SMALK    _2016-11-29 14_42_18.csv')

DF <- DF[DF$Species == 'Gadus morhua',]
plot(DF$IndWgt ~ DF$LngtClass, col = DF$Year)

# W = a L^b

# log(Wt) = log(a) + b*log(L) + E

DF$lWt <- log(DF$IndWgt)
DF$lL  <- log(DF$LngtClass)

plot(DF$lWt ~ DF$lL, col = DF$Year)

lm1 <- glm(lWt ~ lL + Year, data = DF)

summary(lm1)

plot(DF$lWt ~ DF$lL)
abline(lm1, col = 'red')

predDF <- expand.grid(lL = log(seq(min(DF$LngtClass), max(DF$LngtClass), 10)), Year = seq(min(DF$Year), max(DF$Year), 1))

predDF$lWt <- predict(lm1, newdata = predDF)

predDF$L <- exp(predDF$lL)
predDF$Wt <- exp(predDF$lWt)

plot(DF$IndWgt ~ DF$LngtClass, col = DF$Year)
lines(predDF$Wt ~ predDF$L, col = predDF$Year)

library(ggplot2)

ggplot(DF, aes(x = LngtClass, y = IndWgt)) + geom_point() +
	facet_wrap(~ Year) + geom_line(data = predDF, aes(x = L, y = Wt))

