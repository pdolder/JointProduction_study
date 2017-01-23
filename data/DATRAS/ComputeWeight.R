###############################################################
### Convert the at length data from the survey to a biomass ###
###############################################################

# packages
library(dplyr) 
library(ggplot2)
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

###########################
# Now look at door spread #
###########################

# Covert numeric variables so we can explore the covariates
HH$SweepLngt   <- as.numeric(HH$SweepLngt) ;  HH$HaulDur     <- as.numeric(HH$HaulDur)
HH$DoorSpread  <- as.numeric(HH$DoorSpread) ; HH$Depth       <- as.numeric(HH$Depth)
HH$Netopening  <- as.numeric(HH$Netopening) ; HH$Warplngt    <- as.numeric(HH$Warplngt)
HH$Warpdia     <- as.numeric(HH$Warpdia) ;    HH$DoorSurface <- as.numeric(HH$DoorSurface)
HH$DoorWgt     <- as.numeric(HH$DoorWgt) ;    HH$WingSpread  <- as.numeric(HH$WingSpread)
HH$KiteDim     <- as.numeric(HH$KiteDim);     HH$TowDir      <- as.numeric(HH$TowDir)
HH$GroundSpeed <- as.numeric(HH$GroundSpeed); HH$SpeedWater  <- as.numeric(HH$SpeedWater)
HH$SurCurDir   <- as.numeric(HH$SurCurDir) ;  HH$SurCurSpeed <- as.numeric(HH$SurCurSpeed)
HH$BotCurDir   <- as.numeric(HH$BotCurDir);   HH$BotCurSpeed <- as.numeric(HH$BotCurSpeed)
HH$WindDir     <- as.numeric(HH$WindDir) ;    HH$WindSpeed   <- as.numeric(HH$WindSpeed)
HH$SwellDir    <- as.numeric(HH$SwellDir) ;   HH$SwellHeight <- as.numeric(HH$SwellHeight)
HH$SurTemp     <- as.numeric(HH$SurTemp) ;    HH$BotTemp     <- as.numeric(HH$BotTemp)
HH$SurSal      <- as.numeric(HH$SurSal) ;     HH$BotSal      <- as.numeric(HH$BotSal)

HH[HH == -9] <- NA

ggplot(HH, aes(x = Depth, y = DoorSpread)) + geom_point() + 
	theme_bw() + ggtitle('Relationship between depth of gear and
     door spread')

library(gridExtra)

p1 <- ggplot(HH, aes(x = Depth, y = DoorSpread)) + geom_point(aes(colour = factor(DoorWgt))) +
theme_bw() + ggtitle('..with door weight') + theme(legend.position = 'top')

p2 <- ggplot(HH, aes(x = Depth, y = DoorSpread)) + geom_point(aes(colour = Warplngt)) + theme_bw() + 
	ggtitle('..with warp length') + theme(legend.position = 'top')

p3 <- ggplot(HH, aes(x = Depth, y = DoorSpread)) + geom_point(aes(colour = factor(Warpdia))) + theme_bw() + 
	ggtitle('..with warp diameter') + theme(legend.position = 'top')

p4 <- ggplot(HH, aes(x = Depth, y = DoorSpread)) + geom_point(aes(colour = factor(SweepLngt) )) + theme_bw() + 
	ggtitle('..with warp sweep length') + theme(legend.position = 'top')

grid.arrange(p1, p2, p3, p4, ncol = 2)

library(mgcv)

# Without covariate
m1 <- gam(DoorSpread ~ Depth, data = HH)
#summary(m1)

# With all covariate, no interactions
m2 <- gam(DoorSpread ~ Depth + factor(DoorWgt) + Warplngt + factor(Warpdia) +
	  factor(SweepLngt), data =   HH)
#summary(m2)

### full interactions
m3 <- gam(DoorSpread ~ Depth * factor(DoorWgt) * Warplngt * factor(Warpdia) *
	  factor(SweepLngt), data =   HH)
#summary(m3)

HHresid <- filter(HH, !is.na(DoorWgt), !is.na(Warplngt), 
		  !is.na(Warpdia),  !is.na(SweepLngt), 
		  !is.na(Depth), !is.na(DoorSpread))

HHresid$residm3   <- resid(m3)
HHresid$predictm3 <- fitted(m3)

ggplot(HHresid, aes(x = predictm3, y = residm3)) + geom_point() +
	geom_smooth(method = 'loess', col = 'red') +
	theme_bw() + ggtitle('fitted values against residuals') +
	geom_hline(yintercept = 0)

p1 <- ggplot(HHresid, aes(x = factor(DoorWgt), y = residm3)) + geom_boxplot() +
	theme_bw()

p2 <- ggplot(HHresid, aes(x = Warplngt, y = residm3)) + geom_point() +
	geom_smooth(method = 'loess', colour = 'red') + theme_bw() +
	geom_hline(yintercept = 0)

p3 <- ggplot(HHresid, aes(x = factor(Warpdia), y = residm3)) + geom_boxplot() +
	theme_bw()

p4 <- ggplot(HHresid, aes(x = factor(SweepLngt), y = residm3)) + geom_boxplot() +
	theme_bw()

grid.arrange(p1, p2, p3, p4, ncol = 2)

HH$PredSpread <- predict(m3, newdata = HH)

ggplot(HH, aes(x = DoorSpread, y = PredSpread )) + geom_point(colour = 'grey') +
	geom_abline(slope = 1, intercept = 0, col = 'red') + theme_bw() +
	ylab('Predicted door spread') + xlab('Measured door spred') +
	ggtitle('Door spread predictions against measurements')

## Take the or predicted spreads to calculate the swept area...

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

# need as numeric
an <- as.numeric
HL$LngtClass  <- an(HL$LngtClass)
HL$HLNoAtLngt <- an(HL$HLNoAtLngt)
HL$SubFactor  <- an(HL$SubFactor)

# Deal with different length codes - standarise to cm
HL$LngtClass[(HL$LngtClass == 2460 & HL$SpeciesName == 'Merlangius merlangus')] <-  HL$LngtClass[(HL$LngtClass == 2460 & HL$SpeciesName == 'Merlangius merlangus')] / 10 ## Dodgy datapoint!
HL$LngtClass[HL$LngtCode == ". "] <- HL$LngtClass[HL$LngtCode == '. ']/10
HL$LngtClass[HL$LngtCode == 0] <- HL$LngtClass[HL$LngtCode == 0]/10

# Round down length classes & add 0.5
HL$LngtClass[HL$LngtCode != "5"] <- round(HL$LngtClass[HL$LngtCode != "5"])
HL$LngtClass[HL$LngtCode != "5"] <- HL$LngtClass[HL$LngtCode != "5"]+0.5

# boxplot(HL$LngtClass ~ HL$SpeciesName)

########################################
## Now raise with the Model predictions
load('LengthWeightPredictGadoids.RData')

# Filter to the 3 gadoids
HL$Species <- HL$SpeciesName
spp <- c('Gadus morhua','Merlangius merlangus','Melanogrammus aeglefinus')
HL <- filter(HL, Species %in% spp)

# Add log(length)
HL$lL  <- log(HL$LngtClass * 10)

# Predict log weight
HL$LogWtLength <- predict(lm2, newdata = HL)
HL$WtLength <- exp(HL$LogWtLength) # convert back to weight in grams

HL$Wt <- HL$WtLength * HL$HLNoAtLngt * HL$SubFactor # Total weight in g
HL$Wt <- HL$Wt / 1000 # Weight in Kg
HL$Wt <- HL$Wt * corr.fact  # bias correct


#########################################

## And aggregate across lengths
# split into Ju and Ad

HL$SpeciesName <- ifelse(HL$SpeciesName == 'Gadus morhua' & HL$LngtClass <  34.5, paste(HL$SpeciesName,'Juv', sep = '_'), 
ifelse(HL$SpeciesName == 'Gadus morhua' & HL$LngtClass >= 34.5, paste(HL$SpeciesName,'Adu', sep = '_'),
ifelse(HL$SpeciesName == 'Melanogrammus aeglefinus' & HL$LngtClass <  29.5, paste(HL$SpeciesName,'Juv', sep = '_'),
ifelse(HL$SpeciesName == 'Melanogrammus aeglefinus' & HL$LngtClass >= 29.5, paste(HL$SpeciesName,'Adu', sep = '_'),
ifelse(HL$SpeciesName == 'Merlangius merlangus' & HL$LngtClass <  26.5, paste(HL$SpeciesName,'Juv', sep = '_'),
ifelse(HL$SpeciesName == 'Merlangius merlangus' & HL$LngtClass >= 26.5, paste(HL$SpeciesName,'Adu', sep = '_'),
ifelse(HL$SpeciesName == 'Merluccius merluccius' & HL$LngtClass <  26.5, paste(HL$SpeciesName,'Juv', sep = '_'),
ifelse(HL$SpeciesName == 'Merluccius merluccius' & HL$LngtClass >= 26.5, paste(HL$SpeciesName,'Adu', sep = '_'),
ifelse(HL$SpeciesName == 'Pollachius virens' & HL$LngtClass <  34.5, paste(HL$SpeciesName,'Juv', sep = '_'),
ifelse(HL$SpeciesName == 'Pollachius virens' & HL$LngtClass >= 34.5, paste(HL$SpeciesName,'Adu', sep = '_'),
ifelse(HL$SpeciesName == 'Pleuronectes platessa' & HL$LngtClass < 26.5, paste(HL$SpeciesName,'Juv', sep = '_'),
ifelse(HL$SpeciesName == 'Pleuronectes platessa' & HL$LngtClass >= 26.5, paste(HL$SpeciesName,'Adu', sep = '_'),
ifelse(HL$SpeciesName == 'Solea solea' & HL$LngtClass < 23.5, paste(HL$SpeciesName,'Juv', sep = '_'),
ifelse(HL$SpeciesName == 'Solea solea' & HL$LngtClass >= 23.5, paste(HL$SpeciesName,'Adu', sep = '_'),
ifelse(HL$SpeciesName == 'Lepidorhombus whiffiagonis' & HL$LngtClass >= 19.5, paste(HL$SpeciesName,'Adu', sep = '_'),
ifelse(HL$SpeciesName == 'Lepidorhombus whiffiagonis' & HL$LngtClass <  19.5, paste(HL$SpeciesName,'Juv', sep = '_'),
ifelse(HL$SpeciesName == 'dicentratchus labrax' & HL$LngtClass <  35.5, paste(HL$SpeciesName,'Juv', sep = '_'),
ifelse(HL$SpeciesName == 'Dicentrarchus labrax' & HL$LngtClass >= 35.5, paste(HL$SpeciesName,'Adu', sep = '_'),
       paste(HL$SpeciesName,'All', sep ='_')))))))))))))))))))


DF <- HL[!is.na(HL$Wt),]
DF <- DF %>% group_by(Survey, Quarter, Country, Ship, Gear, StNo, HaulNo, Year, SpeciesName) %>%
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

## Create a haul record for each species
HH <- merge(x = HH, y = data.frame(SpeciesName = unique(DF$SpeciesName)))

# Join on the catch data
DF2 <- full_join(x = HH, y = DF)
DF2$Kg[is.na(DF2$Kg)] <- 0  #NAs are zero catches of the species
# Standardise hauls per 30 mins - no, use as effort measure
#DF2$HaulDur <- an(DF2$HaulDur)
#DF2$Kg <- an(as.character(DF2$Kg))
#DF2$Kg <- (DF2$Kg / DF2$HaulDur) * 30

# Subset to variables of interest
DF <- DF2[c('Survey','Ship','StNo','HaulNo','Year', 'Month','SpeciesName','HaulLatMid','HaulLonMid','HaulDur', 'SweptArea', 'SweptAreaAdj','Kg')]


# Remove marginal areas
DF <- filter(DF, HaulLonMid < -2 & HaulLonMid > -12)
DF <- filter(DF, HaulLatMid >  48 & HaulLatMid < 52)


ggplot(DF, aes(x = HaulLonMid, y = HaulLatMid)) + geom_point(aes(size = log(Kg))) + facet_wrap(~SpeciesName, ncol = 2)

# Save

# save(DF, file = file.path('..','CelticSurveyFormatted.RData'))
save(DF, file = file.path('..', 'Cleaned','CelticSurveyFormattedSize.RData'))


