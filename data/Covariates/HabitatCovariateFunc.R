######################################################################
## Generic function for assigning the habitat classification to the ##
#################### knot generated within VAST ######################

HabAssignFunc <- function(KmeansCenters = NULL, zone = 29, locationHabMap = NULL, nameHabMap = NULL) {
library(rgdal)
library(VAST)
# Create a dataframe of the knots
DF <- data.frame(X = KmeansCenters[,'E_km'], Y = KmeansCenters[,'N_km'])
attr(DF, 'projection') = 'UTM'
attr(DF, "zone") <- zone

LLs <- PBSmapping::convUL(DF)

HabMap <- readOGR(dsn = file.path(locationHabMap), layer = nameHabMap)

# joint the spatial points..
LLs <- SpatialPoints(LLs)
proj4string(LLs) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')

join <- over(x = LLs, y = HabMap)

LLs <- SpatialPointsDataFrame(LLs, join)
KmeanHab <- data.frame(Habitat = LLs$Substrate)

return(KmeanHab)

}
