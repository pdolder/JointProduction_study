#########################################
##### Explore and format the CEFAS ######
############ survey data ################

library(ggmap); library(ggplot2); library(dplyr)

rm(list=ls())

FSS <- read.csv(file = 'WesternSurveys_V20160905.dat')
nrow(FSS)
summary(FSS)

######################################
# Plot shoot / haul by survey

Stations <- group_by(FSS, fldSeriesName, fldCruiseName, Year, Month, Day, Time, fldValidityCode, fldTowDuration) %>% 
	    summarise(ShootLat = mean(fldShotLatDecimalDegrees),
		  ShootLon = mean(fldShotLonDecimalDegrees),
		  HaulLat  = mean(fldHaulLatDecimalDegrees),
		  HaulLon  = mean(fldHaulLonDecimalDegrees)) %>% as.data.frame()
######################################

# Keep only valid hauls
table(Stations$fldValidityCode)
Stations <- filter(Stations, fldValidityCode == 'V')

##################################################################
## There are some tows which are clearly too large a distance...
## so to clean the data
# but note tow distance varies greatly over time - so standardise
##################################################################

## https://www.r-bloggers.com/great-circle-distance-calculations-in-r/i

# Convert degrees to radians
deg2rad <- function(deg) return(deg*pi/180)

# Calculates the geodesic distance between two points specified by radian
# latitude/longitude using the
# Haversine formula (hf)
gcd.hf <- function(long1, lat1, long2, lat2) {
	  R <- 6371 # Earth mean radius [km]
 	  delta.long <- (long2 - long1)
   	  delta.lat <- (lat2 - lat1)
    	  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
      	  c <- 2 * asin(min(1,sqrt(a)))
      	  d = R * c
          return(d) # Distance in km
}
############################3

Stations$Dist <- mapply(gcd.hf, long1 = deg2rad(Stations$ShootLon),
       	         lat1  = deg2rad(Stations$ShootLat),
	         long2 = deg2rad(Stations$HaulLon),
	         lat2  = deg2rad(Stations$HaulLat))

summary(Stations$Dist)

Stations$DistStand <- (Stations$Dist / Stations$fldTowDuration) * 60

ggplot(Stations, aes(DistStand)) + geom_histogram() +
	facet_wrap(~fldSeriesName, scale = 'free') 

# Remove any over or under the SE of median distance for the survey (robust
# detection of outliers
# https://www.r-bloggers.com/absolute-deviation-around-the-median/)

StationsClean    <- group_by(Stations, fldSeriesName) %>% 
	            summarise(median = median(DistStand),
			      mean   = mean(DistStand),
		    MAD = mad(DistStand, center = median(DistStand))) %>%
		    as.data.frame()

StationsClean$Up <- StationsClean$median + 5 * StationsClean$MAD
StationsClean$Lo <- StationsClean$median - 5 * StationsClean$MAD

## Now add the upper and lower thresholds to the stations
Stations$LoThres <- StationsClean$Lo[match(Stations$fldSeriesName,
					   StationsClean$fldSeriesName)]
Stations$UpThres <- StationsClean$Up[match(Stations$fldSeriesName,
					   StationsClean$fldSeriesName)]

Stations$InTol   <- ifelse(Stations$DistStand >= Stations$LoThres & Stations$DistStand <= Stations$UpThres, 'KEEP', 'LOSE')

table(Stations$InTol)

ggplot(filter(Stations, InTol == 'KEEP'), aes(DistStand)) + geom_histogram() +
	facet_wrap(~fldSeriesName, scale = 'free') 


########################################################################

map <- get_map("Celtic Sea", zoom = 4)

pdf(file = file.path("CEFASSurveyStations.pdf"), paper = 'a4r')
print(ggmap(map) + geom_segment(data = Stations,aes(x = ShootLon, xend = HaulLon, 
					            y = ShootLat, yend = HaulLat, colour = fldSeriesName)) +
	scale_colour_manual(values = c("red","darkgreen",'darkblue','orange','purple')) + 
	facet_wrap(~InTol))
dev.off()

#######

## Now match and raise the data

FSS$InTol <- Stations$InTol[match(paste(FSS$fldSeriesName, FSS$Year, FSS$Month, FSS$Day, FSS$Time),
				  paste(Stations$fldSeriesName, Stations$Year, Stations$Month, Stations$Day, Stations$Tim))]

nrow(FSS)
FSS <- filter(FSS, InTol == 'KEEP')
nrow(FSS)

# Convert all lengths to cm and round to 5cm size class
FSS$fldLengthGroup <- (FSS$fldLengthGroup / 10) + 0.5

# load a/b parameters
## Add a and b parameters for length-weight

# Only gadoids
FSS <- filter(FSS, fldScientificName %in% c('GADUS MORHUA', 'MELANOGRAMMUS AEGLEFINUS', 'MERLANGIUS MERLANGUS'))

## Load the modelled length weight relationships....
load(file.path('..', 'DATRAS', 'LengthWeightPredictGadoids.RData'))

# Add log length
FSS$lL <- log(FSS$fldLengthGroup * 10)

# Scientific names to small case except first letter
FSS$Species <- paste(toupper(substring(FSS$fldScientificName, 1, 1)), tolower(substring(FSS$fldScientificName, 2, 1000)), sep = '')

# Predict log weight
FSS$LogWtLength <- predict(lm2, newdata = FSS)

FSS$WtLength <- exp(FSS$LogWtLength) # convert back to weight in grams

FSS$Wt <- FSS$WtLength * FSS$Numbers # Total weight in g
FSS$Wt <- FSS$Wt / 1000 # Weight in Kg
FSS$Wt <- FSS$Wt * corr.fact  # bias correct

FSS <- FSS[!is.na(FSS$Wt),]  ## Lack length measurements

COD <- filter(FSS, fldScientificName == 'GADUS MORHUA')

boxplot(COD$Wt / COD$Numbers ~ COD$fldLengthGroup) # Looks OK

# Now merge in the station details: lat, lon etc..
# midpoint of haul locations - small enough distances to not worry about
# spherical distances
an <- as.numeric
FSS$HaulLatMid <- (an(FSS$fldShotLatDecimalDegrees) + an(FSS$fldHaulLatDecimalDegrees)) / 2 
FSS$HaulLonMid <- (an(FSS$fldShotLonDecimalDegrees) + an(FSS$fldHaulLonDecimalDegrees )) / 2

## Calculate the swept area per gear
## for beam trawls its easy, for otter trawls need to include the doorspread
## for effective swept area

## Add distances to the FSS file

FSS$Dist <- mapply(gcd.hf, long1 = deg2rad(FSS$fldShotLonDecimalDegrees),
       	         lat1  = deg2rad(FSS$fldShotLatDecimalDegrees),
	         long2 = deg2rad(FSS$fldHaulLonDecimalDegrees),
	         lat2  = deg2rad(FSS$fldHaulLatDecimalDegrees))

Surveys <- unique(FSS$fldGearDescription)

## No details for otter trawl deployment, so use the standard 87m doorspread
FSS$GearWidth <- ifelse(FSS$fldGearDescription == Surveys[1], 4,  ifelse(
			FSS$fldGearDescription == Surveys[2], 4,  ifelse(
			FSS$fldGearDescription == Surveys[3], 4,  ifelse(
			FSS$fldGearDescription == Surveys[4], 87, ifelse(
			FSS$fldGearDescription == Surveys[5], 4,  ifelse(
			FSS$fldGearDescription == Surveys[6], 87, ifelse(
			FSS$fldGearDescription == Surveys[7], 4,  ifelse(
			FSS$fldGearDescription == Surveys[8], 4,  ifelse(
			FSS$fldGearDescription == Surveys[9], 4,  ifelse(
			FSS$fldGearDescription == Surveys[10], 4, ifelse(
			FSS$fldGearDescription == Surveys[11], 87,ifelse(
			FSS$fldGearDescription == Surveys[12], 87,ifelse(
			FSS$fldGearDescription == Surveys[13], 87,ifelse(
			FSS$fldGearDescription == Surveys[14], 87,ifelse(
			FSS$fldGearDescription == Surveys[15], 87,ifelse(
			FSS$fldGearDescription == Surveys[16], 87,ifelse(
			FSS$fldGearDescription == Surveys[17], 3, ifelse(
			FSS$fldGearDescription == Surveys[18], 87,ifelse(
			FSS$fldGearDescription == Surveys[19], 87,ifelse(
			FSS$fldGearDescription == Surveys[20], 87,ifelse(
			FSS$fldGearDescription == Surveys[21], 87, ifelse(
			FSS$fldGearDescription == Surveys[22], 2, ifelse(
			FSS$fldGearDescription == Surveys[23], NA,ifelse(
			FSS$fldGearDescription == Surveys[24], NA,ifelse(
			FSS$fldGearDescription == Surveys[25], NA,ifelse(
			FSS$fldGearDescription == Surveys[26], 2.5, ifelse(
			FSS$fldGearDescription == Surveys[27], 87,  ifelse(
			FSS$fldGearDescription == Surveys[28], 87,  ifelse(
			FSS$fldGearDescription == Surveys[29], 3,   ifelse(
			FSS$fldGearDescription == Surveys[30], 3,NA))))))))))))))))))))))))))))))

FSS$SweptArea <- FSS$Dist * (FSS$GearWidth / 1000)

plot(FSS$SweptArea ~ FSS$fldSeriesName) 

## Adjust swept area for gear efficiencies, after Piet et al for roundfish:

# BT: 0.19
# OT: 0.22 - 0.54 (Juv, ad). 0.38

FSS <- FSS[!is.na(FSS$GearWidth),]

FSS$SweptAreaAdjFac <- sapply(FSS$GearWidth, function(x) {
	       if(x %in% c(2.5, 3, 4))  return(0.19) 
	       if(x == 87)  return(0.38)    
	       else return(NA)
		 })

FSS$SweptAreaAdj <- FSS$SweptArea * FSS$SweptAreaAdjFac

plot(FSS$SweptAreaAdj ~ FSS$fldSeriesName)

# Aggregate

# split into Ju and Ad

FSS$Species <- ifelse(FSS$Species == 'Gadus morhua' & FSS$fldLengthGroup <  34.5, paste(FSS$Species,'Juv', sep = '_'), 
ifelse(FSS$Species == 'Gadus morhua' & FSS$fldLengthGroup >= 34.5, paste(FSS$Species,'Adu', sep = '_'),
ifelse(FSS$Species == 'Melanogrammus aeglefinus' & FSS$fldLengthGroup <  29.5, paste(FSS$Species,'Juv', sep = '_'),
ifelse(FSS$Species == 'Melanogrammus aeglefinus' & FSS$fldLengthGroup >= 29.5, paste(FSS$Species,'Adu', sep = '_'),
ifelse(FSS$Species == 'Merlangius merlangus' & FSS$fldLengthGroup <  26.5, paste(FSS$Species,'Juv', sep = '_'),
ifelse(FSS$Species == 'Merlangius merlangus' & FSS$fldLengthGroup >= 26.5, paste(FSS$Species,'Adu', sep = '_'),
ifelse(FSS$Species == 'Merluccius merluccius' & FSS$fldLengthGroup <  26.5, paste(FSS$Species,'Juv', sep = '_'),
ifelse(FSS$Species == 'Merluccius merluccius' & FSS$fldLengthGroup >= 26.5, paste(FSS$Species,'Adu', sep = '_'),
ifelse(FSS$Species == 'Pollachius virens' & FSS$fldLengthGroup <  34.5, paste(FSS$Species,'Juv', sep = '_'),
ifelse(FSS$Species == 'Pollachius virens' & FSS$fldLengthGroup >= 34.5, paste(FSS$Species,'Adu', sep = '_'),
ifelse(FSS$Species == 'Pleuronectes platessa' & FSS$fldLengthGroup < 26.5, paste(FSS$Species,'Juv', sep = '_'),
ifelse(FSS$Species == 'Pleuronectes platessa' & FSS$fldLengthGroup >= 26.5, paste(FSS$Species,'Adu', sep = '_'),
ifelse(FSS$Species == 'Solea solea' & FSS$fldLengthGroup < 23.5, paste(FSS$Species,'Juv', sep = '_'),
ifelse(FSS$Species == 'Solea solea' & FSS$fldLengthGroup >= 23.5, paste(FSS$Species,'Adu', sep = '_'),
ifelse(FSS$Species == 'Lepidorhombus whiffiagonis' & FSS$fldLengthGroup >= 19.5, paste(FSS$Species,'Adu', sep = '_'),
ifelse(FSS$Species == 'Lepidorhombus whiffiagonis' & FSS$fldLengthGroup <  19.5, paste(FSS$Species,'Juv', sep = '_'),
ifelse(FSS$Species == 'Dicentratchus labrax' & FSS$fldLengthGroup <  35.5, paste(FSS$Species,'Juv', sep = '_'),
ifelse(FSS$Species == 'Dicentrarchus labrax' & FSS$fldLengthGroup >= 35.5, paste(FSS$Species,'Adu', sep = '_'),
       paste(FSS$Species,'All', sep ='_')))))))))))))))))))


FSS <- group_by(FSS, fldSeriesName, Year, Month, fldPrimeStation, fldCruiseStationNumber, HaulLatMid, HaulLonMid, fldTowDuration, SweptArea, SweptAreaAdj, Species) %>% summarise(Kg = sum(Wt)) %>% as.data.frame()

FSS <- FSS[c('fldSeriesName','Year', 'Month','fldPrimeStation','fldCruiseStationNumber','HaulLatMid','HaulLonMid','fldTowDuration', 'SweptArea', 'SweptAreaAdj' ,'Species','Kg')]

## Trim to only keep data within core Celtic Sea

FSS <- filter(FSS, HaulLonMid > -12 &  HaulLonMid < -2) # remove extreme Lons
FSS <- filter(FSS, HaulLatMid >  48 &  HaulLatMid < 52) # remove extreme Lats

plot(FSS$SweptAreaAdj ~ FSS$fldSeriesName)
boxplot(FSS$SweptAreaAdj ~ FSS$Year)

table(FSS$Month, FSS$Year, FSS$fldSeriesName)

# save(FSS, file = file.path('..','CelticSurvey2Formatted.RData'))
save(FSS, file = file.path('..', 'Cleaned','CelticSurvey2FormattedSize.RData'))

