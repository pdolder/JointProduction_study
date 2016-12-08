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

Wt <- data.frame(Survey = c(DWt$Survey, CWt$fldSeriesName),
		 Year   = c(DWt$Year  , CWt$Year),
		 Month  = c(DWt$Month , CWt$Month),
		 HaulNo = c(DWt$HaulNo, CWt$fldCruiseStationNumber),

