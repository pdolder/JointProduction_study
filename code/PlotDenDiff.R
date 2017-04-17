#####################################################################
### Code for plotting differences in density for cod, had and whg ###
#####################################################################

rm(list = ls())

load(file.path('..', 'results', 'CovariatesAtKnot.RData'))
load(file.path('..', 'results', '2017-04-08_M2', 'Save.RData'))

MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

cod <- Save$Report$Index_xcyl[,1,,]
had <- Save$Report$Index_xcyl[,9,,]
whg <- Save$Report$Index_xcyl[,11,,]

# Check overall index
par(mfrow = c(1,3))
plot(colSums(cod), type = 'b')
plot(colSums(had), type = 'b') 
plot(colSums(whg), type = 'b')

# Now standardise the 2015 values
cod <- cod[,26] / sum(cod[,26]) * 100
had <- had[,26] / sum(had[,26]) * 100 
whg <- whg[,26] / sum(whg[,26]) * 100

# Make the plot dataframe

plotDF <- data.frame(x = Spatial_List$Kmeans$centers[,1], 
		    y= Spatial_List$Kmeans$centers[,2],
		    knot = 1:nrow(Spatial_List$Kmeans$centers))

DF <- data.frame(X = Centers[,'E_km'], Y = Centers[,'N_km'])
attr(DF, 'projection') = 'UTM'
attr(DF, "zone") <- 29 

LLs <- PBSmapping::convUL(DF)

plotDF$Lon <- LLs$Y
plotDF$Lat <- LLs$X

## Add the densities
plotDF$cod <- cod
plotDF$had <- had
plotDF$whg <- whg

## Area densities relative to... 
plotDF$cod_had  <- plotDF$cod - plotDF$had
plotDF$cod_whg  <- plotDF$cod - plotDF$whg
plotDF$had_whg  <- plotDF$had - plotDF$whg

#############
## the Map ##
#############

library(ggplot2); library(cowplot)
library(mapplots)
library(rworldmap)
library(rworldxtra)
library(data.table)
library(broom)
library(RColorBrewer)

mapdata =ggplot2::fortify(rworldmap::getMap(resolution = 'high'))

xlim = c(-12, -2)
ylim = c(48, 52)
 # if you want to fast subset by lat/long you can do the following:
 # subset to just regions in xlim and ylim, see
# http://stackoverflow.com/a/16574176
 mapdata = data.table(mapdata)
  mapdata = mapdata[mapdata[,.I[any(long %between% xlim) & any(lat %between%
 ylim)], by = list(group)]$V1]

  # make geoms
  coast.poly = geom_polygon(data=mapdata, aes(x=long, y=lat, group=group), colour="#999999", fill="#999999", lwd=0.2)
  coast.outline = geom_path(data=mapdata, aes(x=long, y=lat, group=group), colour="#000000", lwd=0.2)


## The map for all areas

DF2 <- MapDetails_List[['PlotDF']]

DF2$cod_had <- plotDF$cod_had[match(DF2$x2i, plotDF$knot)]
DF2$cod_whg <- plotDF$cod_whg[match(DF2$x2i, plotDF$knot)]
DF2$had_whg <- plotDF$had_whg[match(DF2$x2i, plotDF$knot)]

## The plots
p1 <- ggplot() + geom_point(data = DF2, aes(x = Lon, y =  Lat, colour = cod_had), size = 0.02) + geom_point(data = DF2, aes(x = Lat, y = Lon), size = 0.5, shape = 3) + 
	coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = 'bottom', legend.text = element_text(angle = -90), legend.title = element_blank()) + 
	scale_colour_gradient2(low = 'red', mid = 'white', high = 'blue', midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar")

p2 <- ggplot() + geom_point(data = DF2, aes(x = Lon, y = Lat, colour = cod_whg), size = 0.02) + geom_point(data = DF2, aes(x = Lat, y = Lon), size = 0.5, shape = 3) +
	coast.poly + coast.outline + coord_quickmap(xlim, ylim) + theme(legend.position = 'bottom', legend.text = element_text(angle = -90), legend.title = element_blank()) + 
	scale_colour_gradient2(low = 'red', mid = 'white', high = 'blue', midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar")

p3 <- ggplot() + geom_point(data = DF2, aes(x = Lon, y = Lat, colour = had_whg), size = 0.02) + geom_point(data = DF2, aes(x = Lat, y = Lon), size = 0.5, shape = 3) +
	coast.poly + coast.outline + coord_quickmap(xlim, ylim) + theme(legend.position = 'bottom', legend.text = element_text(angle = -90), legend.title = element_blank()) + 
	scale_colour_gradient2(low = 'red', mid = 'white', high = 'blue', midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar")

#plot_grid(p1,p2, p3, labels = c('(a) cod:haddock', '(b) cod:whiting', '(c) haddock:whiting'), vjust = c(1.5, 1.5,1.5), ncol = 3)

DF3 <- reshape2::melt(DF2, id = c('Lat', 'Lon', 'x2i', 'Include'))

ggplot() + geom_point(data = DF3, aes(x = Lon, y = Lat, colour = value), size = 0.02) + geom_point(data = DF3, aes(x = Lat, y = Lon), size = 0.5, shape = 3) +
	coast.poly + coast.outline + coord_quickmap(xlim, ylim) + theme(legend.position = 'bottom', legend.text = element_text(angle = -90), legend.title = element_blank()) + 
	scale_colour_gradient2(low = 'red', mid = 'white', high = 'blue', midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar") + facet_wrap(~variable)

ggsave(file.path('..', 'plots', 'Density_Differences.png'), width = 12, height = 4)

