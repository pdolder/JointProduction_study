######################################################################
######### Generic function for assigning a depth to the ##############
#################### knot generated within VAST ######################
######################################################################

DepthAssignFunc <- function(KmeansCenters = NULL, zone = 29, locationDepths = NULL) {
library(rgdal)
library(VAST)
# Create a dataframe of the knots
DF <- data.frame(X = KmeansCenters[,'E_km'], Y = KmeansCenters[,'N_km'])
attr(DF, 'projection') = 'UTM'
attr(DF, "zone") <- zone

# Get the Lon Lats in decimal degrees
LLs <- PBSmapping::convUL(DF)

load(locationDepths)

# Find the nearest depth measurement
Depthfun <- function(x,y) {
BD$Depth[which(abs(BD$Lon - x) == min(abs(BD$Lon - x)) &
	 abs(BD$Lat - y) == min(abs(BD$Lat - y)))]
}

LLs$Depth <- mapply(x = LLs$X, y = LLs$Y, FUN = Depthfun)

return(LLs$Depth)

}
