###################################################
# Some plots of the survey data for diagnostics ###
###################################################

rm(list = ls())

library(dplyr)
library(ggplot2)

##################################################

# Load in data
load(file.path('CelticSurveyFormattedSize.RData')) # Datras data by weight
load(file.path('CelticSurvey2FormattedSize.RData')) # Cefas data by weight

DWt <- DF; CWt <- FSS; rm(DF, FSS) # Rename to avoid confusion

load(file.path('DATRAS_No_At_Length.RData')) # Datras data by length 
load(file.path('Cefas_No_At_Length.RData')) # Cefas data by length 

Dln <- DF; Cln <- FSS; rm(DF, FSS)

##################################################

## Combine the datasets

Wt <- data.frame(Survey    = c(DWt$Survey      , CWt$fldSeriesName),
		 Year      = c(DWt$Year        , CWt$Year),
		 Month     = c(DWt$Month       , CWt$Month),
		 HaulNo    = c(DWt$HaulNo      , CWt$fldCruiseStationNumber),
		 Lon       = c(DWt$HaulLonMid  , CWt$HaulLonMid),
		 Lat       = c(DWt$HaulLatMid  , CWt$HaulLatMid),
		 HaulDur   = c(DWt$HaulDur     , CWt$fldTowDuration),
		 SweptArea = c(DWt$SweptAreaAdj, CWt$SweptAreaAdj),
		 Species   = c(DWt$Species,      CWt$Species),
		 Kg        = c(DWt$Kg          , CWt$Kg))
rm(DWt, CWt)

Ln <- data.frame(Survey    = c(Dln$Survey      , Cln$fldSeriesName),
		 Year      = c(Dln$Year        , Cln$Year),
		 Month     = c(Dln$Month       , Cln$Month),
		 HaulNo    = c(Dln$HaulNo      , Cln$fldCruiseStationNumber),
		 Lon       = c(Dln$HaulLonMid  , Cln$HaulLonMid),
		 Lat       = c(Dln$HaulLatMid  , Cln$HaulLatMid),
		 HaulDur   = c(Dln$HaulDur     , Cln$fldTowDuration),
		 SweptArea = c(Dln$SweptAreaAdj, Cln$SweptAreaAdj),
		 Length    = c(Dln$LngtClass   , Cln$fldLengthGroup),
		 Number    = c(Dln$Number      , Cln$Number),
		 Kg        = c(Dln$Kg          , Cln$Kg))


plotDF <- filter(Wt, Year == 2011, Species == 'Gadus morhua_Adu')

library(maps)

map(regions = c('ireland', 'uk', 'france'), xlim = c(-12, 2), ylim = c(48, 52))
points(plotDF$Lat[plotDF$Kg != 0] ~ plotDF$Lon[plotDF$Kg != 0], col = 'blue')
points(plotDF$Lat[plotDF$Kg == 0] ~ plotDF$Lon[plotDF$Kg == 0], pch = '+', col = 'red')

