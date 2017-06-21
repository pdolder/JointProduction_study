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
Temps <- data.frame(knot = HabDF$knot, Lon = HabDF$Lon, Lat = HabDF$Lat, Year = rep(1990:2015, each = 250), Temp = NA)

## http://lukemiller.org/index.php/2014/11/extracting-noaa-sea-surface-temperatures-with-ncdf4/
source('NOAA_OISST_ncdf4.R')

# The land mask file
LM <- 'ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/lsmask.oisst.v2.nc'
tmpFileLM <- tempfile()
download.file(LM, tmpFileLM)


yrs <- 1990:2015

for (y in yrs)  {

rm(tmpFile, file); gc() # remove files and clear memory

file <- paste('ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.mean.',y,'.v2.nc', sep = "")
tmpFile <- tempfile()
download.file(file, tmpFile)

# extract the data for each knot 
 
for (i in 1:nrow(Temps[Temps$Year == y,])) {

	tmpTemps <- Temps[Temps$Year == y,]

ssts = extractOISSTdaily(tmpFile,
			 tmpFileLM,
			 lonW=tmpTemps[tmpTemps$knot == i,'Lon']+180,
			 lonE=tmpTemps[tmpTemps$knot == i,'Lon']+180,
			 latS=tmpTemps[tmpTemps$knot == i,'Lat'],
			 latN=tmpTemps[tmpTemps$knot == i,'Lat'],
			 date1 = paste(y,"-10-01",sep=""),
			 date2 = paste(y,"-12-31",sep=""))
Temps[Temps$Year == y & Temps$knot == i,'Temp']  <-  mean(ssts)
} # end of knots loop


} # end of year loop

save(Temps, file = 'YearlyMeanSSTatKnot.RData')

