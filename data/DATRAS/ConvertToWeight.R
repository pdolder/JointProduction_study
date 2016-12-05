###############################################################
### Convert the at length data from the survey to a biomass ###
###############################################################

# packages
library(dplyr) 

load('CelticSurveyData.RData') # pre-downloaded data, HH is station, HL is catch

## Some initial cleaning
HH <- filter(HH, HaulVal == 'V') # only valid hauls

##################################################################
# Want to calculate the effective swept area..
# From the distance * Doorspread
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
an <- as.numeric
HH$Dist <- mapply(gcd.hf, long1 = deg2rad(an(HH$ShootLong)),
       	         lat1  = deg2rad(an(HH$ShootLat)),
	         long2 = deg2rad(an(HH$HaulLong)),
	         lat2  = deg2rad(an(HH$HaulLat)))

plot(an(HH$Distance[(HH$Distance != -9)])/1000 ~ HH$Dist[(HH$Distance != -9)])
## Looks good - use the calculated estimates ##

# Now look at door spread
plot(HH$DoorSpread[an(HH$DoorSpread) != -9]) # Differences by survey?

par(mfrow = c(1,2))
plot(HH$DoorSpread[an(HH$DoorSpread) != -9 & HH$Survey == 'EVHOE'], ylim = c(0, 150), col = HH$Year[an(HH$DoorSpread) != -9 & HH$Survey == 'EVHOE'])
plot(HH$DoorSpread[an(HH$DoorSpread) != -9 & HH$Survey == 'IE-IGFS'], ylim = c(0, 150), col = HH$Year[an(HH$DoorSpread) != -9 & HH$Survey == 'IE-IGFS'])
# Yes, but also a split within the IE-IGFS, not by year though

# Seems to be a relationship between depth and doorspread
plot(HH$DoorSpread[an(HH$DoorSpread) !=-9 & an(HH$Depth) != -9] ~
     HH$Depth[an(HH$DoorSpread) !=-9 & an(HH$Depth) != -9])

table(an(HH$Depth[an(HH$DoorSpread) == -9]))
#  As only 5 records have both no depth and no doorspread, lets fill the door
#  spread based on this relationship

Spread <- an(HH$DoorSpread[an(HH$DoorSpread) !=-9 & an(HH$Depth) != -9])
Depth  <- an(HH$Depth[an(HH$DoorSpread) !=-9 & an(HH$Depth) != -9])

model <- lm(Spread ~ poly(Depth,4))
summary(model)

plot(fitted(model), residuals(model))

pred <- predict(model, newdata = data.frame(Depth = seq(min(Depth), max(Depth), l = 1000), 
					    Spread = seq(min(Spread), max(Spread), l = 1000)))

plot(Spread ~ Depth)
lines(pred, col = 'red')

## looks OK, now predict the missing door spreads..
HH$PredSpread <- predict(model, newdata = data.frame(Depth = an(HH$Depth), Spread = an(HH$DoorSpread)))

## Take the actuals or predicted spreads to calculate the swept area...

HH$SweptArea <- HH$Dist * HH$PredSpread/1000

plot(HH$SweptArea)

## OK 

## Now for adjustment factor from Piet et al

HH$SweptAreaAdjFac <- 0.38

HH$SweptAreaAdj <- HH$SweptArea * HH$SweptAreaAdjFac

################################################################
## Add species names
load('DatrasSpeciesCodes.RData')
HL$SpeciesName <- DatrasSpeciesCodes$scientific.name[match(HL$SpecCode, DatrasSpeciesCodes$code_number)]

## Add a and b parameters for length-weight
load('length.weight.RData')

HL$SpeciesA <- length.weight$a[match(HL$SpeciesName, length.weight$scientific.name)]
HL$SpeciesB <- length.weight$b[match(HL$SpeciesName, length.weight$scientific.name)]

# Weight at length class

# need as numeric
an <- as.numeric
HL$LngtClass  <- an(HL$LngtClass)
HL$HLNoAtLngt <- an(HL$HLNoAtLngt)
HL$SubFactor  <- an(HL$SubFactor)


# Deal with different length codes - standarise to cm
HL$LngtClass[HL$LngtCode == ". "] <- HL$LngtClass[HL$LngtCode == '. ']/10
HL$LngtClass[HL$LngtCode == 0] <- HL$LngtClass[HL$LngtCode == 0]/10

# Round down length classes & add 0.5
HL$LngtClass[HL$LngtCode != "5"] <- round(HL$LngtClass[HL$LngtCode != "5"])
HL$LngtClass[HL$LngtCode != "5"] <- HL$LngtClass[HL$LngtCode != "5"]+0.5

# boxplot(HL$LngtClass ~ HL$SpeciesName)

# Now raise the length data to weights
# a*L^b

HL$Wt <- (HL$SpeciesA * HL$LngtClass)^HL$SpeciesB * (HL$HLNoAtLngt * HL$SubFactor) * 100 # in kg

# COD <- filter(HL, SpeciesName == 'Gadus morhua')
# boxplot((COD$Wt/COD$HLNoAtLngt) ~ COD$LngtClass) 

## And aggregate across lengths

# split into Ju and Ad

HL$SpeciesName <- toupper(HL$SpeciesName)

HL$SpeciesName <- ifelse(HL$SpeciesName == 'GADUS MORHUA' & HL$LngtClass <  34.5, paste(HL$SpeciesName,'Juv', sep = '_'), 
ifelse(HL$SpeciesName == 'GADUS MORHUA' & HL$LngtClass >= 34.5, paste(HL$SpeciesName,'Adu', sep = '_'),
ifelse(HL$SpeciesName == 'MELANOGRAMMUS AEGLEFINUS' & HL$LngtClass <  29.5, paste(HL$SpeciesName,'Juv', sep = '_'),
ifelse(HL$SpeciesName == 'MELANOGRAMMUS AEGLEFINUS' & HL$LngtClass >= 29.5, paste(HL$SpeciesName,'Adu', sep = '_'),
ifelse(HL$SpeciesName == 'MERLANGIUS MERLANGUS' & HL$LngtClass <  26.5, paste(HL$SpeciesName,'Juv', sep = '_'),
ifelse(HL$SpeciesName == 'MERLANGIUS MERLANGUS' & HL$LngtClass >= 26.5, paste(HL$SpeciesName,'Adu', sep = '_'),
ifelse(HL$SpeciesName == 'MERLUCCIUS MERLUCCIUS' & HL$LngtClass <  26.5, paste(HL$SpeciesName,'Juv', sep = '_'),
ifelse(HL$SpeciesName == 'MERLUCCIUS MERLUCCIUS' & HL$LngtClass >= 26.5, paste(HL$SpeciesName,'Adu', sep = '_'),
ifelse(HL$SpeciesName == 'POLLACHIUS VIRENS' & HL$LngtClass <  34.5, paste(HL$SpeciesName,'Juv', sep = '_'),
ifelse(HL$SpeciesName == 'POLLACHIUS VIRENS' & HL$LngtClass >= 34.5, paste(HL$SpeciesName,'Adu', sep = '_'),
ifelse(HL$SpeciesName == 'PLEURONECTES PLATESSA' & HL$LngtClass < 26.5, paste(HL$SpeciesName,'Juv', sep = '_'),
ifelse(HL$SpeciesName == 'PLEURONECTES PLATESSA' & HL$LngtClass >= 26.5, paste(HL$SpeciesName,'Adu', sep = '_'),
ifelse(HL$SpeciesName == 'SOLEA SOLEA' & HL$LngtClass < 23.5, paste(HL$SpeciesName,'Juv', sep = '_'),
ifelse(HL$SpeciesName == 'SOLEA SOLEA' & HL$LngtClass >= 23.5, paste(HL$SpeciesName,'Adu', sep = '_'),
ifelse(HL$SpeciesName == 'LEPIDORHOMBUS WHIFFIAGONIS' & HL$LngtClass >= 19.5, paste(HL$SpeciesName,'Adu', sep = '_'),
ifelse(HL$SpeciesName == 'LEPIDORHOMBUS WHIFFIAGONIS' & HL$LngtClass <  19.5, paste(HL$SpeciesName,'Juv', sep = '_'),
ifelse(HL$SpeciesName == 'DICENTRATCHUS LABRAX' & HL$LngtClass <  35.5, paste(HL$SpeciesName,'Juv', sep = '_'),
ifelse(HL$SpeciesName == 'DICENTRARCHUS LABRAX' & HL$LngtClass >= 35.5, paste(HL$SpeciesName,'Adu', sep = '_'),
       paste(HL$SpeciesName,'All', sep ='_')))))))))))))))))))


DF <- HL[!is.na(HL$Wt),]
DF <- DF %>% group_by(Survey, Quarter, Country, Ship, Gear, StNo, HaulNo,Year, SpeciesName) %>%
	summarise(Kg = sum(Wt)) %>% as.data.frame()

# Now merge in the station details: lat, lon etc..
# midpoint of haul locations - small enough distances to not worry about
# spherical distances
HH$HaulLatMid <- (an(HH$ShootLat) + an(HH$HaulLat)) / 2 
HH$HaulLonMid <- (an(HH$ShootLon) + an(HH$HaulLon )) / 2

# Fix blank spaces in variables...
DF$Survey <- gsub(' ','',DF$Survey)
DF$Gear   <- gsub(' ','',DF$Gear)
DF$Ship   <- gsub(' ','',DF$Ship)
DF$StNo   <- gsub(' ','',DF$StNo)

DF2 <- left_join(x = DF, y = HH)

# Standardise hauls per 30 mins - no, use as effort measure
#DF2$HaulDur <- an(DF2$HaulDur)
#DF2$Kg <- an(as.character(DF2$Kg))
#DF2$Kg <- (DF2$Kg / DF2$HaulDur) * 30

# Subset to variables of interest
DF <- DF2[c('Survey','Ship','StNo','HaulNo','Year','SpeciesName','HaulLatMid','HaulLonMid','HaulDur', 'SweptArea', 'SweptAreaAdj','Kg')]


# Remove marginal areas
DF <- filter(DF, HaulLonMid < -2 & HaulLonMid > -12)
DF <- filter(DF, HaulLatMid >  48 & HaulLatMid < 52)

# Save

# save(DF, file = file.path('..','CelticSurveyFormatted.RData'))
save(DF, file = file.path('..','CelticSurveyFormattedSize.RData'))


