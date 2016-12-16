#################################################################
## Calculations to estimate relative survey catch rates at length 
## based on raw numbers caught per survey
## From Clark 2013 Direct calculation of relative fishery and 
## survey selectivities
#################################################################
## House cleaning
rm(list = ls())
library(dplyr)
library(ggplot2)
#################################################################
## Data
# Load in data
load(file.path('..','Cleaned', 'DATRAS_No_At_Length.RData')) # Datras data by length 
load(file.path('..','Cleaned','Cefas_No_At_Length.RData')) # Cefas data by length 

Dln <- DF; Cln <- FSS; rm(DF, FSS)

##################################################

## Combine the datasets
Ln <- data.frame(Survey    = c(Dln$Survey      , as.character(Cln$fldSeriesName)),
		 Year      = c(Dln$Year        , Cln$Year),
		 Month     = c(Dln$Month       , Cln$Month),
		 HaulNo    = c(Dln$HaulNo      , Cln$fldCruiseStationNumber),
		 Lon       = c(Dln$HaulLonMid  , Cln$HaulLonMid),
		 Lat       = c(Dln$HaulLatMid  , Cln$HaulLatMid),
		 HaulDur   = c(Dln$HaulDur     , Cln$fldTowDuration),
		 SweptArea = c(Dln$SweptAreaAdj, Cln$SweptAreaAdj),
		 Species   = c(as.character(Dln$SpeciesName) , Cln$Species),
		 Length    = c(Dln$LngtClass   , Cln$Length),
		 Number    = c(Dln$Number      , Cln$Numbers),
		 Kg        = c(Dln$Kg          , Cln$Kg))

rm(Dln, Cln)
gc()

###################################################

## First the IE-IFGS relative to the EVHOE only
## Choose comparable years, 2003 - 2015
DF <- filter(Ln, Survey %in% c('EVHOE', 'IE-IGFS'), Year %in% c(2003:2015))

## Sum across all survey locations and reformat the data

DF2 <- group_by(DF, Survey, Species, Length) %>% 
	summarise(N = sum(Number)) %>% 
	reshape2::dcast(Species + Length ~ Survey, value.var = 'N')

colnames(DF2)[4] <- 'IE_IGFS'

## Plot the survey catch at length
ggplot(DF2, aes(x = Length)) + geom_line(aes(y = EVHOE), col = 'red') + geom_line(aes(y = IE_IGFS), col = 'green') +
	facet_wrap(~ Species, ncol = 1, scale = 'free_y')

## Proportion at length in the survey

TotNo <- group_by(DF2, Species) %>% summarise(EVHOE = sum(EVHOE, na.rm = T), IE_IGFS = sum(IE_IGFS, na.rm = T))
DF3 <- group_by(DF2, Species, Length) %>% summarise(EVHOE = EVHOE/sum(EVHOE, na.rm = T))

DF2$TotalEvo <- TotNo$EVHOE[match(DF2$Species, TotNo$Species)]
DF2$TotalIE <- TotNo$IE_IGFS[match(DF2$Species, TotNo$Species)]

DF2$PropEVHOE <- DF2$EVHOE / DF2$TotalEvo
DF2$PropIE <- DF2$IE_IGFS / DF2$TotalIE

ggplot(DF2, aes(x = Length)) + geom_line(aes(y = PropEVHOE), col = 'red') +
	geom_line(aes(y = PropIE), col = 'green') +
	facet_wrap(~Species, ncol = 1)

######################
### Calculation is ###
# Ps1,a / Ps2,a = S1,a/S1,2 * Sum(Ssb * Pb)/Sum(Scb * Pb)
#####################

## For the RHS we simply divide the total biomasses by each other
B <- group_by(DF, Survey, Species) %>%
	summarise(B = sum(Kg)) %>%
	reshape2::dcast(Species ~ Survey, value.var = 'B')

# Ratio of IE-IGFS to EVHOE
B_ratio <- data.frame(Species = B[,1],ratio = B[,3]/B[,2])

## Now calculate the ratios

DF2$Ratios <- DF2[,'IE_IGFS'] / DF2[,'EVHOE']


Cod <- filter(DF2, Species == 'Gadus morhua')

Cod$Ratio * B_ratio[B_ratio$Species == 'Gadus morhua','ratio']
