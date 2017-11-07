############################################################
## Downloading the NOAA 0.25 * 0.25 high res SST dataset ###
## ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/
############################################################

library(ncdf4)
library(ncdf4.helpers)
library(PCICt)

###################################
## Data frame to store the temps ##
###################################

## Read in Spatial List
load(file.path('..', '..','results', 'CovariatesAtKnot.RData'))

MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

HabDF <- data.frame(x = Spatial_List$Kmeans$centers[,1], 
		    y= Spatial_List$Kmeans$centers[,2],
		    knot = 1:nrow(Spatial_List$Kmeans$centers),
		    Habitat = Hab$Habitat, Depth = Depths)

DF <- data.frame(X = Centers[,'E_km'], Y = Centers[,'N_km'])
attr(DF, 'projection') = 'UTM'
attr(DF, "zone") <- 29 

LLs <- PBSmapping::convUL(DF)

HabDF$Lat <- LLs$Y
HabDF$Lon <- LLs$X

## Data frame to store the temperatures
Temps <- data.frame(knot = HabDF$knot, Lon = HabDF$Lon, Lat = HabDF$Lat, 
		    Year = rep(1990:2015, each = 250), TempMean = NA, TempCumSum = NA)

# Generate set of grid cell latitudes (center of cell) from south to north
	lats = seq(-89.875,89.875,0.25)
	# Generate set of grid cell longitudes (center of cell)
	lons = seq(0.125,359.875,0.25) 

yrs <- 1990:2015

for (y in yrs)  {
	print(y)

if(y == 1990) { tmpFile <- nc_open(file.path('NOAA_TempData', '1990_1995.nc'))}
if(y == 1996) { tmpFile <- nc_open(file.path('NOAA_TempData', '1996_2001.nc'))}
if(y == 2002) { tmpFile <- nc_open(file.path('NOAA_TempData', '2002_2008.nc'))}
if(y == 2009) { tmpFile <- nc_open(file.path('NOAA_TempData', '2009_2015.nc'))}

if(y %in% c(1990, 1996, 2002, 2009)) {
# Extract the data
sst <- ncvar_get(tmpFile, "sst")
lat <- ncvar_get(tmpFile, "lat")
lon <- ncvar_get(tmpFile, "lon")
ncdates <- tmpFile$dim$time$vals
ncdates <- as.Date(ncdates, origin = '1800-1-1')

dimnames(sst) <- list(lon = lon, lat = lat, time= as.character(ncdates))
sst = ifelse(sst== 32767, NA, sst)
ssts <- as.data.frame.table(sst)
ssts$lat <- as.numeric(as.character(ssts$lat))
ssts$lon <- as.numeric(as.character(ssts$lon))

}

# extract the data for each knot 
for (i in 1:nrow(Temps[Temps$Year == y,])) {
	print(i)

	# knot nearest temp locations
	kLon <- Temps[Temps$Year == y & Temps$knot == i, 'Lon'] +360 
	kLon <- lons[which.min(abs(kLon - lons))]
	kLat <- Temps[Temps$Year == y & Temps$knot == i, 'Lat']
	kLat <- lats[which.min(abs(kLat - lats))]

	
# subset the ssts based on these lat & lons & keep only dates
	# in the May - Nov range
ssts_k <- ssts[ssts$lat == kLat & ssts$lon == kLon & substr(ssts$time,1,4) == y  &
	       substr(ssts$time,6,7) %in% c("05","06","07","08","09","10","11"),'Freq']

# save the mean
Temps[Temps$Year == y & Temps$knot == i,'TempMean']  <-  mean(ssts_k, na.rm = T)
Temps[Temps$Year == y & Temps$knot == i,'TempCumSum']  <-  sum(ssts_k, na.rm = T)

} # end of knots loop

if(y %in% c(1995, 2001, 2008, 2015)) nc_close(tmpFile)

} # end of year loop

save(Temps, file = 'YearlyMeanandCumSumSSTatKnot.RData')

library(ggplot2)

ggplot(Temps, aes(x = Year, y = TempMean)) + geom_line(aes(colour = factor(knot)))
