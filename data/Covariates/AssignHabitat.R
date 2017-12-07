
## This can be generalised, but basically takes the Kmeans knots and assigns a
## habitat based on joining the EMOD data to produce a
#### dataframe of habitats associated with the kmean centers for a covariate

library(VAST)
library(rgdal)

# Load the kmeans
load(file.path('Kmeans-100.RData'))

# Create dataframe of kmeans and convert back to LLs
DF <- data.frame(X = Kmeans$centers[,'E_km'], Y = Kmeans$centers[,'N_km'])
attr(DF, 'projection') = 'UTM'
attr(DF, "zone") <-  29 
LLs <- PBSmapping::convUL(DF)

# Distances between sample points
hist(dist(coordinates(LLs)))

## http://rstudio-pubs-static.s3.amazonaws.com/7993_6b081819ba184047802a508a7f3187cb.html

# Read in frame
HabMap <- readOGR(file.path('201208_EUSeaMap_Atlantic_Habitats'),'201208_EUSeaMap_Atlantic_Habitats')

## Join on the spatial points

LLs <- SpatialPoints(LLs)
proj4string(LLs) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')

join <- over(x = LLs, y = HabMap)

LLs <- SpatialPointsDataFrame(LLs, join)
head(LLs)

table(LLs$substrate)

KmeanHab <- data.frame(Habitat = LLs$substrate)

print(KmeanHab)

## 61 is missing...

## Some plotting 
coordsDF <- data.frame(Lon = coordinates(LLs)[,1], Lat = coordinates(LLs)[,2])

library(maps)
library(mapdata)
map('worldHires',c('UK','Ireland','France'), xlim = c(-12,2), ylim = c(46,54))
points(coordsDF$Lon, coordsDF$Lat, col = KmeanHab$Habitat)
#points(coordsDF$Lon[61], coordsDF$Lat[61], col = 'red', pch = 'x')
legend(x = -2, y = 48, legend = levels(KmeanHab$Habitat), col = 1:length(unique(KmeanHab$Habitat)), pch = 'o', cex = 0.8)
savePlot(file = 'HabitatAssignment.png')
## Its very close to land, assign the nearest point
KmeanHab$Habitat[61] <- 'Rock or other hard substrata'

# plot(HabMap, col = 'grey', border = NULL, axes = T, xlim = c(-12,2), ylim = c(46,54))


## Check the colours are right...
ggDF <- cbind(coordsDF,KmeanHab)

ggplot(ggDF, aes(x = Lon, y = Lat)) + geom_point(aes(colour = Habitat))

# Save for input
save(KmeanHab, file = 'KmeansHab.RData')


